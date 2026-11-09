import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_clean_data(file_path=None):
    """
    Loads and standardizes the 'Clean' sheet of 358_DATA.xlsx.
    """
    if file_path is None:
        file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'metadata', '358_DATA.xlsx')
    
    # Load sheet
    df = pd.read_excel(file_path, sheet_name='Clean')
    
    # Clean string columns: strip spaces and standardize categories
    df['FIRST_SYMPTOM'] = df['FIRST_SYMPTOM'].astype(str).str.strip().str.title()
    
    # MRI Categories
    df['MRI_CATEGORIES'] = df['MRI_CATEGORIES'].astype(str).str.strip().str.upper()
    mri_map = {
        'UNKNOWN': 'Unknown',
        'CEREBELLAR ATROPHY': 'Cerebellar Atrophy',
        'CEREBRAL AND CEREBELLAR ATROPHY': 'Cerebral & Cerebellar Atrophy',
        'CEREBRAL ATROPHY': 'Cerebral Atrophy',
        'NORMAL': 'Normal',
        'OTHERS': 'Others'
    }
    df['MRI_CATEGORIES'] = df['MRI_CATEGORIES'].map(mri_map).fillna('Unknown')
    
    # Parasympathetic Reactivity
    df['PARASYMPATHETIC_ACTIVITY_CATEGORIES'] = df['PARASYMPATHETIC_ACTIVITY_CATEGORIES'].astype(str).str.strip().str.upper()
    para_map = {
        'NORMAL': 'Normal',
        'MILD LOSS': 'Mild Loss',
        'MODERATE LOSS': 'Moderate Loss',
        'SEVERE LOSS': 'Severe Loss',
        'MILD TO MODERATE LOSS': 'Moderate Loss',
        'UNKNOWN': 'Unknown'
    }
    df['PARASYMPATHETIC_ACTIVITY_CATEGORIES'] = df['PARASYMPATHETIC_ACTIVITY_CATEGORIES'].map(para_map).fillna('Unknown')
    
    # Sympathetic Reactivity
    df['SYMPATHETIC_ACTIVITY_CATEGORIES'] = df['SYMPATHETIC_ACTIVITY_CATEGORIES'].astype(str).str.strip().str.upper()
    symp_map = {
        'NORMAL': 'Normal',
        'MILD LOSS': 'Mild Loss',
        'MODERATE LOSS': 'Moderate Loss',
        'SEVERE LOSS': 'Severe Loss',
        'MODERATE TO SEVERE LOSS': 'Moderate Loss',
        'UNKNOWN': 'Unknown'
    }
    df['SYMPATHETIC_ACTIVITY_CATEGORIES'] = df['SYMPATHETIC_ACTIVITY_CATEGORIES'].map(symp_map).fillna('Unknown')
    
    # Bradykinesia
    df['BRADYKINESIA_CATEGORIES'] = df['BRADYKINESIA_CATEGORIES'].astype(str).str.strip().str.upper()
    ep_map = {'0': 'Absent', '1': 'Present', 'UNKNOWN': 'Unknown'}
    df['BRADYKINESIA_CATEGORIES'] = df['BRADYKINESIA_CATEGORIES'].map(ep_map).fillna('Unknown')
    
    # Rigidity
    df['RIGIDITY_CATEGORIES'] = df['RIGIDITY_CATEGORIES'].astype(str).str.strip().str.upper()
    df['RIGIDITY_CATEGORIES'] = df['RIGIDITY_CATEGORIES'].map(ep_map).fillna('Unknown')
    
    # Plantar Response
    df['PLANTAR_CATEGORIES'] = df['PLANTAR_CATEGORIES'].astype(str).str.strip().str.upper()
    plantar_map = {
        'FLEXOR': 'Flexor (Normal)',
        'EXTENSOR': 'Extensor (Abnormal)',
        'EQUIVOCAL': 'Equivocal',
        'UNKNOWN': 'Unknown'
    }
    df['PLANTAR_CATEGORIES'] = df['PLANTAR_CATEGORIES'].map(plantar_map).fillna('Unknown')
    
    return df

def plot_first_symptom(df, out_dir):
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')
    
    counts = df['FIRST_SYMPTOM'].value_counts()
    percentages = (counts / len(df)) * 100
    
    bars = ax.barh(counts.index, counts.values, color='#2b5c8f', edgecolor='none', height=0.6)
    ax.set_title('First Clinical Symptoms (n = 356)', fontsize=14, fontweight='bold', pad=15, color='#2c3e50')
    ax.set_xlabel('Number of Patients', fontsize=11, fontweight='semibold')
    ax.set_ylabel('Symptom', fontsize=11, fontweight='semibold')
    ax.set_xlim(0, max(counts.values) * 1.30)
    
    for bar, pct in zip(bars, percentages):
        width = bar.get_width()
        ax.text(width + 3, bar.get_y() + bar.get_height()/2, f'{int(width)} ({pct:.1f}%)', 
                ha='left', va='center', fontsize=10, fontweight='bold', color='#333333')
        
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '08a_first_symptom.png'), dpi=600, bbox_inches='tight')
    plt.close()

def plot_mri_categories(df, out_dir):
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')
    
    # Order categories logically
    order = ['Normal', 'Cerebellar Atrophy', 'Cerebral Atrophy', 'Cerebral & Cerebellar Atrophy', 'Others', 'Unknown']
    counts = df['MRI_CATEGORIES'].value_counts().reindex(order).fillna(0)
    percentages = (counts / len(df)) * 100
    
    colors = ['#2ca02c', '#d62728', '#ff7f0e', '#9467bd', '#bcbd22', '#7f7f7f']
    bars = ax.bar(counts.index, counts.values, color=colors, edgecolor='none', width=0.5)
    
    ax.set_title('Brain MRI Classifications (n = 356)', fontsize=14, fontweight='bold', pad=15, color='#2c3e50')
    ax.set_ylabel('Number of Patients', fontsize=11, fontweight='semibold')
    ax.set_xticks(range(len(counts.index)))
    ax.set_xticklabels(counts.index, rotation=30, ha='right', fontsize=10, fontweight='semibold')
    ax.set_ylim(0, max(counts.values) * 1.15)
    
    for bar, pct in zip(bars, percentages):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height + 3, f'{int(height)}\n({pct:.1f}%)', 
                ha='center', va='bottom', fontsize=9, fontweight='bold', color='#333333')
        
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '08b_mri_categories.png'), dpi=600, bbox_inches='tight')
    plt.close()

def plot_autonomic_activity(df, out_dir):
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')
    
    categories = ['Normal', 'Mild Loss', 'Moderate Loss', 'Severe Loss', 'Unknown']
    para_counts = df['PARASYMPATHETIC_ACTIVITY_CATEGORIES'].value_counts().reindex(categories).fillna(0)
    symp_counts = df['SYMPATHETIC_ACTIVITY_CATEGORIES'].value_counts().reindex(categories).fillna(0)
    
    para_pct = (para_counts / len(df)) * 100
    symp_pct = (symp_counts / len(df)) * 100
    
    x = np.arange(len(categories))
    width = 0.35
    
    rects1 = ax.bar(x - width/2, para_counts.values, width, label='Parasympathetic Activity', color='#2b5c8f')
    rects2 = ax.bar(x + width/2, symp_counts.values, width, label='Sympathetic Activity', color='#d95f02')
    
    ax.set_title('Autonomic Dysfunction Profiles (n = 356)', fontsize=14, fontweight='bold', pad=15, color='#2c3e50')
    ax.set_ylabel('Number of Patients', fontsize=11, fontweight='semibold')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=10, fontweight='semibold')
    ax.legend(frameon=True, facecolor='#ffffff')
    ax.set_ylim(0, max(max(para_counts.values), max(symp_counts.values)) * 1.15)
    
    # Add values on top of bars
    for i, (r1, r2) in enumerate(zip(rects1, rects2)):
        h1 = r1.get_height()
        h2 = r2.get_height()
        ax.text(r1.get_x() + r1.get_width()/2, h1 + 2, f'{int(h1)}\n({para_pct.iloc[i]:.1f}%)', ha='center', va='bottom', fontsize=8, color='#333333', fontweight='semibold')
        ax.text(r2.get_x() + r2.get_width()/2, h2 + 2, f'{int(h2)}\n({symp_pct.iloc[i]:.1f}%)', ha='center', va='bottom', fontsize=8, color='#333333', fontweight='semibold')
        
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '08c_autonomic_activity.png'), dpi=600, bbox_inches='tight')
    plt.close()

def plot_extrapyramidal(df, out_dir):
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')
    
    categories = ['Absent', 'Present', 'Unknown']
    brady_counts = df['BRADYKINESIA_CATEGORIES'].value_counts().reindex(categories).fillna(0)
    rigid_counts = df['RIGIDITY_CATEGORIES'].value_counts().reindex(categories).fillna(0)
    
    brady_pct = (brady_counts / len(df)) * 100
    rigid_pct = (rigid_counts / len(df)) * 100
    
    x = np.arange(len(categories))
    width = 0.35
    
    rects1 = ax.bar(x - width/2, brady_counts.values, width, label='Bradykinesia (Extrapyramidal)', color='#4682b4')
    rects2 = ax.bar(x + width/2, rigid_counts.values, width, label='Rigidity (Extrapyramidal)', color='#e07a5f')
    
    ax.set_title('Extrapyramidal Sign Prevalence (n = 356)', fontsize=14, fontweight='bold', pad=15, color='#2c3e50')
    ax.set_ylabel('Number of Patients', fontsize=11, fontweight='semibold')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=10, fontweight='semibold')
    ax.legend(frameon=True, facecolor='#ffffff')
    ax.set_ylim(0, max(max(brady_counts.values), max(rigid_counts.values)) * 1.15)
    
    for i, (r1, r2) in enumerate(zip(rects1, rects2)):
        h1 = r1.get_height()
        h2 = r2.get_height()
        ax.text(r1.get_x() + r1.get_width()/2, h1 + 2, f'{int(h1)}\n({brady_pct.iloc[i]:.1f}%)', ha='center', va='bottom', fontsize=8, color='#333333', fontweight='semibold')
        ax.text(r2.get_x() + r2.get_width()/2, h2 + 2, f'{int(h2)}\n({rigid_pct.iloc[i]:.1f}%)', ha='center', va='bottom', fontsize=8, color='#333333', fontweight='semibold')
        
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '08d_extrapyramidal.png'), dpi=600, bbox_inches='tight')
    plt.close()

def plot_plantar(df, out_dir):
    fig, ax = plt.subplots(figsize=(6, 6))
    fig.patch.set_facecolor('#fafafa')
    
    categories = ['Flexor (Normal)', 'Extensor (Abnormal)', 'Equivocal', 'Unknown']
    counts = df['PLANTAR_CATEGORIES'].value_counts().reindex(categories).fillna(0)
    
    colors = ['#2ca02c', '#d62728', '#ff7f0e', '#7f7f7f']
    
    wedges, texts = ax.pie(
        counts.values, 
        labels=[f'{cat}\n(n={int(val)}, {val/sum(counts)*100:.1f}%)' for cat, val in zip(categories, counts.values)], 
        startangle=140, 
        colors=colors, 
        wedgeprops=dict(width=0.4, edgecolor='w', linewidth=2),
        textprops=dict(fontsize=10, color='#2c3e50', fontweight='semibold')
    )
    
    ax.text(0, 0, f'Plantar\nResponse\n(n={len(df)})', ha='center', va='center', fontsize=12, fontweight='bold', color='#2c3e50')
    ax.set_title('Pyramidal Reflex: Plantar Response (n = 356)', fontsize=14, fontweight='bold', pad=15, color='#2c3e50')
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '08e_plantar.png'), dpi=600, bbox_inches='tight')
    plt.close()

def main():
    print("Generating Clinical Phenotypic Proportion Plots...")
    df = load_clean_data()
    
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'plots')
    os.makedirs(out_dir, exist_ok=True)
    
    plot_first_symptom(df, out_dir)
    plot_mri_categories(df, out_dir)
    plot_autonomic_activity(df, out_dir)
    plot_extrapyramidal(df, out_dir)
    plot_plantar(df, out_dir)
    
    print("All proportion plots successfully generated and saved to plots/ directory.")

if __name__ == '__main__':
    main()
