import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from data_loader import load_and_clean_data, get_clinical_data_for_analysis

def plot_symptom_concordance(output_path=None):
    """
    Generates the clinical scale correlation heatmap.
    """
    if output_path is None:
        output_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plots', '07_symptom_concordance.png')
    
    df = load_and_clean_data()
    # Filter rows with valid values for the core clinical scales
    scales = {
        'SARA Score': 'Ataxia Severity (SARA)',
        'INAS Score': 'Non-Ataxic Signs (INAS)',
        'FTMTRS Grand Total': 'Tremor Severity (FTMTRS)',
        'TETRAS Grand Total': 'Tremor Severity (TETRAS)',
        'ADL Total': 'Daily Living Impairment (ADL)'
    }
    
    df_scales = df[list(scales.keys())].dropna()
    # Rename columns for the plot
    df_scales = df_scales.rename(columns=scales)
    
    # Calculate Pearson correlation matrix
    corr_matrix = df_scales.corr()
    
    sns.set_theme(style='white')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11
    
    fig, ax = plt.subplots(figsize=(9, 8))
    fig.patch.set_facecolor('#fafafa')
    
    # Draw heatmap
    # We use a beautiful custom colormap (warm diverging or sequential)
    sns.heatmap(
        corr_matrix, 
        annot=True, 
        fmt=".3f", 
        cmap='Blues', 
        vmin=0.6, 
        vmax=1.0, 
        square=True, 
        linewidths=2, 
        cbar_kws={"shrink": 0.8, "label": "Pearson Correlation Coefficient (r)"},
        annot_kws={"size": 12, "weight": "bold"},
        ax=ax
    )

    for text in ax.texts:
        value = float(text.get_text())
        
        if value >= 0.85:
            text.set_color("white")
        else:
            text.set_color("#1c2833")
    
    # Rotate tick labels for readability
    plt.xticks(rotation=45, ha='right', fontsize=11, fontweight='semibold', color='#2c3e50')
    plt.yticks(rotation=0, fontsize=11, fontweight='semibold', color='#2c3e50')
    
    ax.set_title('Concordance of Clinical Scales (n = 74 patients)', fontsize=15, fontweight='bold', pad=20, color='#2c3e50')
    
    plt.suptitle('Cross-Scale Correlation Heatmap of Motor and Daily Living Metrics', fontsize=16, fontweight='bold', color='#1a252f', y=0.98)
    plt.tight_layout()
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Symptom concordance plot saved to: {output_path}")

if __name__ == '__main__':
    plot_symptom_concordance()
