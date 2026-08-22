import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
from data_loader import load_and_clean_data

def plot_genetic_anticipation(output_path=None):
    """
    Generates the genetic anticipation visualization (UPPER ALLELE size vs. Age of Onset).
    """
    if output_path is None:
        output_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plots', '05_genetic_anticipation.png')
    
    df = load_and_clean_data()
    # Filter rows with both valid UPPER ALLELE and AO
    df_analysis = df.dropna(subset=['UPPER ALLELE', 'AO']).copy()
    
    sns.set_theme(style='whitegrid')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11
    
    fig, ax = plt.subplots(figsize=(9, 7))
    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')
    
    # Calculate correlation and regression parameters
    x = df_analysis['UPPER ALLELE']
    y = df_analysis['AO']
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_sq = r_value ** 2
    
    # Scatter plot colored by SARA Score to show severity
    # SARA Score might be missing for some, we can dropna for SARA as well or color by it if available
    df_scatter = df_analysis.dropna(subset=['SARA Score'])
    
    scatter = ax.scatter(
        df_scatter['UPPER ALLELE'], 
        df_scatter['AO'], 
        c=df_scatter['SARA Score'], 
        cmap='viridis_r', 
        s=80, 
        alpha=0.85, 
        edgecolors='w', 
        linewidths=0.5,
        label='Patients'
    )
    
    # For points that have missing SARA but valid AO/allele, plot in grey
    missing_sara = df_analysis[df_analysis['SARA Score'].isna()]
    if not missing_sara.empty:
        ax.scatter(
            missing_sara['UPPER ALLELE'], 
            missing_sara['AO'], 
            c='#7f8c8d', 
            s=80, 
            alpha=0.6, 
            edgecolors='w', 
            linewidths=0.5,
            label='Patients (missing SARA)'
        )
    
    # Fit line using seaborn regplot or custom line to match colors
    sns.regplot(
        x='UPPER ALLELE', 
        y='AO', 
        data=df_analysis, 
        scatter=False, 
        color='#e74c3c', 
        line_kws={'linewidth': 2.5, 'label': 'Regression Line'}, 
        ax=ax
    )
    
    # Colorbar
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('Ataxia Severity (SARA Score)', rotation=270, labelpad=15, fontweight='semibold')
    
    # Annotation text
    annotation_str = (
        f"Pearson r = {r_value:.3f}\n"
        f"R² = {r_sq:.3f}\n"
        f"p-value = {p_value:.1e}\n"
        f"n = {len(df_analysis)} patients\n"
        f"Equation: AO = {slope:.2f}*(CAG) + {intercept:.2f}"
    )
    ax.text(
        0.05, 0.05, 
        annotation_str, 
        transform=ax.transAxes, 
        fontsize=11, 
        fontweight='semibold', 
        color='#2c3e50',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#ffffff', edgecolor='#bdc3c7', alpha=0.9)
    )
    
    ax.set_title('CAG Repeat Size vs. Age of Onset', fontsize=15, fontweight='bold', pad=15, color='#2c3e50')
    ax.set_xlabel('Expanded CAG Repeat Size (UPPER ALLELE)', fontsize=12, fontweight='semibold')
    ax.set_ylabel('Age of Onset (AO, years)', fontsize=12, fontweight='semibold')
    ax.set_xlim(44, 72)
    ax.set_ylim(20, 70)
    
    plt.suptitle('Linear Relationship and Impact of Repeat Expansion on Onset Age', fontsize=16, fontweight='bold', color='#1a252f', y=0.98)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Anticipation plot saved to: {output_path}")

if __name__ == '__main__':
    plot_genetic_anticipation()
