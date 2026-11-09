import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
from data_loader import load_and_clean_data

def plot_clinical_progression(output_path=None):
    """
    Generates the clinical progression visualization (Disease duration vs. SARA Score).
    """
    if output_path is None:
        output_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plots', '06_clinical_progression.png')
    
    df = load_and_clean_data()
    # Filter rows with valid SARA Score and Disease duration
    df_analysis = df.dropna(subset=['SARA Score', 'Disease duration']).copy()
    
    sns.set_theme(style='whitegrid')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11
    
    fig, ax = plt.subplots(figsize=(9, 7))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')
    
    # Calculate correlation and regression parameters
    x = df_analysis['Disease duration']
    y = df_analysis['SARA Score']
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_sq = r_value ** 2
    
    # Filter rows that also have UPPER ALLELE for color coding
    df_scatter = df_analysis.dropna(subset=['UPPER ALLELE'])
    
    scatter = ax.scatter(
        df_scatter['Disease duration'], 
        df_scatter['SARA Score'], 
        c=df_scatter['UPPER ALLELE'], 
        cmap='coolwarm', 
        s=80, 
        alpha=0.85, 
        edgecolors='w', 
        linewidths=0.5,
        label='Patients'
    )
    
    # For points missing UPPER ALLELE but having valid SARA/duration
    missing_genetics = df_analysis[df_analysis['UPPER ALLELE'].isna()]
    if not missing_genetics.empty:
        ax.scatter(
            missing_genetics['Disease duration'], 
            missing_genetics['SARA Score'], 
            c='#7f8c8d', 
            s=80, 
            alpha=0.6, 
            edgecolors='w', 
            linewidths=0.5,
            label='Patients (missing genetics)'
        )
        
    # Fit regression line
    sns.regplot(
        x='Disease duration', 
        y='SARA Score', 
        data=df_analysis, 
        scatter=False, 
        color='#1b3a5b', 
        line_kws={'linewidth': 2.5, 'label': 'Regression Line'}, 
        ax=ax
    )
    
    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('Expanded CAG Repeat Size (UPPER ALLELE)', rotation=270, labelpad=15, fontweight='semibold')
    
    # Annotation text
    annotation_str = (
        f"Pearson r = {r_value:.3f}\n"
        f"R² = {r_sq:.3f}\n"
        f"p-value = {p_value:.1e}\n"
        f"n = {len(df_analysis)} patients\n"
        f"Progression Rate: {slope:.3f} SARA points/year\n"
        f"Equation: SARA = {slope:.2f}*(Duration) + {intercept:.2f}"
    )
    ax.text(
        0.05, 0.65, 
        annotation_str, 
        transform=ax.transAxes, 
        fontsize=11, 
        fontweight='semibold', 
        color='#2c3e50',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#ffffff', edgecolor='#bdc3c7', alpha=0.9)
    )
    
    ax.set_title('SARA Score vs. Disease Duration', fontsize=15, fontweight='bold', pad=15, color='#2c3e50')
    ax.set_xlabel('Disease Duration (years)', fontsize=12, fontweight='semibold')
    ax.set_ylabel('Ataxia Severity (SARA Score)', fontsize=12, fontweight='semibold')
    ax.set_xlim(-1, 35)
    ax.set_ylim(-1, 35)
    
    plt.suptitle('Predicting Disease Severity over Time in SCA12 Cohort', fontsize=16, fontweight='bold', color='#1a252f', y=0.98)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Progression plot saved to: {output_path}")

if __name__ == '__main__':
    plot_clinical_progression()
