DIR="/mnt/faruq2/lab_users/ankit/nanopore_analysis/MF_5RACE_PPP2R2B_29082025/analysis_long_read"

SJ_DIR="$DIR/short_read_junctions"

for file in "$SJ_DIR"/*.SJ.out.tab; do

    echo "======================================"
    echo "Checking: $(basename "$file")"
    echo "======================================"

    awk '
    BEGIN { FS="\t" }

    NF != 9 {
        malformed++
    }

    NF == 9 && $4 !~ /^[012]$/ {
        invalid_strand++
        if (invalid_strand <= 5)
            print "Invalid strand:", NR, $0
    }

    END {
        print "Malformed rows:", malformed+0
        print "Invalid strand rows:", invalid_strand+0
    }
    ' "$file"

done
