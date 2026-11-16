#!/usr/bin/env python3
"""Create a GSP-truncated, exon-focused GTF from GENCODE or RefSeq GTF.

The target gene is replaced by gene/transcript/exon records. Other records
are preserved verbatim. The primer interval is retained, not removed.

Run using:

From blastn against GSP: GGCGGTTTTCCATGTTGAGATCCCAGACTTTG
GSP_START=146600308
GSP_END=146600339

python3 cap_gtf.py \
  /path/to/gencode.annotation.gtf.gz \
  /path/to/refseq.annotation.gtf.gz \
  --primer-mode retained \
  --outdir /path/to/capped_references
"""
import argparse, collections, csv, gzip, pathlib, re

PRIMER = 'GGCGGTTTTCCATGTTGAGATCCCAGACTTTG'
def attrs(s):
    return dict(re.findall(r'(\S+)\s+"?([^";]+)"?\s*;', s))
def op(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)
def target(d):
    return (any(d.get(k)=='PPP2R2B' for k in ('gene_name','gene','gene_id'))
            or d.get('gene_id','').split('.')[0]=='ENSG00000156475')

def read_rows(path):
    with op(path) as f:
        for n,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'):
                yield line, None, None
                continue
            r=line.rstrip('\n').split('\t')
            if len(r)!=9: raise ValueError(f'{path}:{n}: expected 9 fields')
            yield line, r, attrs(r[8])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('gtf', nargs='+', type=pathlib.Path)
    p.add_argument('--primer-start', type=int, default=146600308, help='BLAST genomic start, 1-based inclusive (default: 146600308)')
    p.add_argument('--primer-end', type=int, default=146600339, help='BLAST genomic end, 1-based inclusive (default: 146600339)')
    p.add_argument('--primer-mode', choices=['retained','removed'], default='retained', help='Whether aligned reads retain genomic primer-binding bases')
    p.add_argument('--outdir', type=pathlib.Path, default=pathlib.Path('capped_gtf'))
    a=p.parse_args()
    start,end=a.primer_start,a.primer_end
    if start<1 or end-start+1!=len(PRIMER): p.error('Coordinates must span the 32-base primer')
    cutoff=start if a.primer_mode=='retained' else end+1
    # Preflight all output names before writing anything.
    names=[path.name.removesuffix('.gz').removesuffix('.gtf') for path in a.gtf]
    if len(set(names))!=len(names): p.error('Input basenames must be distinct')
    for path,name in zip(a.gtf,names):
        for suffix in ('.primer_capped.gtf','.audit.tsv'):
            out=a.outdir/(name+suffix)
            if out.exists() or out.resolve()==path.resolve(): p.error(f'Refusing to overwrite {out}')
    a.outdir.mkdir(parents=True,exist_ok=True)
    for path,name in zip(a.gtf,names):
        # Pass 1 stores only PPP2R2B; memory does not grow with the entire GTF.
        groups=collections.defaultdict(list); loci=set()
        for line,r,d in read_rows(path):
            if r is None or not target(d): continue
            loci.add((r[0],r[6]))
            if d.get('transcript_id'): groups[d['transcript_id']].append((r,d,True))
        if not groups: raise ValueError(f'{path}: no PPP2R2B transcripts found')
        if len(loci)!=1: raise ValueError(f'{path}: multiple PPP2R2B loci; select primary assembly annotation')
        chrom,strand=next(iter(loci))
        if strand!='-' or chrom not in ('chr5','5','NC_000005.10'):
            raise ValueError(f'{path}: expected human GRCh38 chromosome 5, minus strand; found {chrom}, {strand}')
        primer=PRIMER
        replacements={}; report=[]
        for tid,items in groups.items():
            exons=[r for r,d,t in items if r[2]=='exon']
            compatible=any(int(r[3])<=start and int(r[4])>=end for r in exons)
            kept=[]
            if compatible:
                for r in exons:
                    z=r.copy(); z[3]=str(max(int(z[3]),cutoff))
                    if int(z[3])<=int(z[4]): kept.append(z)
            kept.sort(key=lambda r:int(r[4]),reverse=True)
            report.append((tid,'kept' if kept else 'excluded_no_primer_spanning_exon',len(exons),len(kept)))
            if not kept: replacements[tid]=[]; continue
            template=next((r for r,d,t in items if r[2]=='transcript'),items[0][0]).copy()
            template[2]='transcript'; template[3]=str(min(int(r[3]) for r in kept)); template[4]=str(max(int(r[4]) for r in kept)); template[7]='.'
            template[8]=re.sub(r'\bexon_number\s+[^;]+;\s*','',template[8])
            template[8]+=' race_capped "true";'
            for i,r in enumerate(kept,1):
                r[7]='.'
                r[8]=re.sub(r'\bexon_number\s+[^;]+;\s*','',r[8])+f' exon_number "{i}";'
                # Reference exon IDs cease to describe the original interval after clipping.
                if int(r[3])==cutoff:
                    r[8]=re.sub(r'\bexon_id\s+[^;]+;\s*','',r[8])
            replacements[tid]=[template]+kept
        out=a.outdir/(name+'.primer_capped.gtf')
        if out.resolve()==path.resolve() or out.exists(): raise ValueError(f'Refusing to overwrite {out}')
        emitted=set()
        with open(out,'w') as f:
            f.write(f'# PPP2R2B assay reference; {chrom}:{start}-{end}; strand=-; primer={primer}; mode={a.primer_mode}; cutoff={cutoff}; target features=transcript,exon\n')
            # Pass 2 copies every other gene's original text verbatim.
            for line,r,d in read_rows(path):
                if r is None or not target(d):
                    f.write(line)
                    continue
                tid=d.get('transcript_id')
                if tid and tid not in emitted:
                    for z in replacements[tid]: f.write('\t'.join(z)+'\n')
                    emitted.add(tid)
        with open(a.outdir/(name+'.audit.tsv'),'w') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(['transcript_id','status','original_exons','retained_exons']); w.writerows(report)
        print(f'{out}: {sum(bool(v) for v in replacements.values())}/{len(groups)} transcripts retained; primer={chrom}:{start}-{end}; cutoff={cutoff}')
if __name__=='__main__': main()
