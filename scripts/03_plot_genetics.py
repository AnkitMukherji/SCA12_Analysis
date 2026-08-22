import os
import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import load_and_clean_data


def plot_allele_distribution(output_path=None):
    """
    Plot the distributions of LOWER and UPPER CAG repeat alleles.

    LOWER ALLELE: smaller allele
    UPPER ALLELE: larger allele

    Highlights the pathogenic expansion threshold at >=43 repeats.
    """

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'plots',
            '03_allele_distribution.png'
        )

    # ------------------------------------------------------------------
    # Load and prepare data
    # ------------------------------------------------------------------

    df = load_and_clean_data()

    df_gen = df.dropna(
        subset=['LOWER ALLELE', 'UPPER ALLELE']
    ).copy()

    # ------------------------------------------------------------------
    # Global style
    # ------------------------------------------------------------------

    sns.set_theme(style='whitegrid')

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.titleweight'] = 'bold'

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # ------------------------------------------------------------------
    # KDE distributions
    # ------------------------------------------------------------------

    sns.kdeplot(
        data=df_gen,
        x='LOWER ALLELE',
        fill=True,
        color='#2b5c8f',
        alpha=0.35,
        linewidth=2,
        label='Lower allele',
        ax=ax
    )

    sns.kdeplot(
        data=df_gen,
        x='UPPER ALLELE',
        fill=True,
        color='#9467bd',
        alpha=0.35,
        linewidth=2,
        label='Upper allele',
        ax=ax
    )

    # ------------------------------------------------------------------
    # Clinically relevant regions
    # ------------------------------------------------------------------

    # Borderline / intermediate range
    ax.axvspan(
        28,
        43,
        color='grey',
        alpha=0.10,
        zorder=0
    )

    # Pathogenic threshold
    ax.axvline(
        43,
        color='#b22222',
        linestyle='--',
        linewidth=1.8,
        label='Expansion threshold (≥43)'
    )

    # ------------------------------------------------------------------
    # Threshold annotation
    # ------------------------------------------------------------------

    ax.text(
        43.8,
        ax.get_ylim()[1] * 0.90,
        '≥43 repeats',
        color='#b22222',
        fontsize=10,
        fontweight='semibold',
        va='top'
    )

    ax.text(
        35.5,
        ax.get_ylim()[1] * 0.08,
        'Intermediate range',
        color='#666666',
        fontsize=9,
        ha='center',
        style='italic'
    )

    # ------------------------------------------------------------------
    # Labels
    # ------------------------------------------------------------------

    ax.set_title(
        'CAG Repeat-Size Distribution',
        fontsize=15,
        fontweight='bold',
        pad=12
    )

    ax.set_xlabel(
        'CAG repeat size (PPP2R2B)',
        fontsize=12,
        labelpad=8, 
        fontweight = 'bold'
    )

    ax.set_ylabel(
        'Density',
        fontsize=12,
        labelpad=8,
        fontweight = 'bold'
    )

    ax.set_xlim(0, 80)

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------

    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        frameon=False,
        fontsize=10
    )

    # ------------------------------------------------------------------
    # Clean appearance
    # ------------------------------------------------------------------

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(
        axis='y',
        alpha=0.25
    )

    ax.grid(
        axis='x',
        visible=False
    )

    plt.tight_layout()
    plt.subplots_adjust(right=0.78)

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------

    os.makedirs(
        os.path.dirname(output_path),
        exist_ok=True
    )

    plt.savefig(
        output_path,
        dpi=600,
        bbox_inches='tight'
    )

    plt.close()

    print(
        f"Allele distribution plot saved to: {output_path}"
    )


if __name__ == '__main__':
    plot_allele_distribution()

import os
import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import load_and_clean_data


def plot_allele_scatter(output_path=None):
    """
    Plot LOWER vs UPPER CAG repeat alleles.

    Individuals with LOWER ALLELE >= 43 are classified as
    double-expansion / homozygous expansion cases.
    """

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'plots',
            '04_allele_scatter.png'
        )

    # ------------------------------------------------------------------
    # Load and prepare data
    # ------------------------------------------------------------------

    df = load_and_clean_data()

    df_gen = df.dropna(
        subset=['LOWER ALLELE', 'UPPER ALLELE']
    ).copy()

    # ------------------------------------------------------------------
    # Classification
    # ------------------------------------------------------------------

    df_gen['Class'] = 'Heterozygous expansion'

    double_expansion = (
        df_gen['LOWER ALLELE'] >= 43
    )

    df_gen.loc[
        double_expansion,
        'Class'
    ] = 'Biallelic expansion'

    # ------------------------------------------------------------------
    # Style
    # ------------------------------------------------------------------

    sns.set_theme(style='whitegrid')

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(8, 7)
    )

    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # ------------------------------------------------------------------
    # Scatter plot
    # ------------------------------------------------------------------

    sns.scatterplot(
        data=df_gen,
        x='LOWER ALLELE',
        y='UPPER ALLELE',
        hue='Class',
        style='Class',
        palette={
            'Heterozygous expansion': '#2b5c8f',
            'Biallelic expansion': '#d62728'
        },
        markers={
            'Heterozygous expansion': 'o',
            'Biallelic expansion': 'X'
        },
        s=75,
        alpha=0.80,
        edgecolor='white',
        linewidth=0.5,
        ax=ax
    )

    # ------------------------------------------------------------------
    # Pathogenic thresholds
    # ------------------------------------------------------------------

    ax.axvline(
        43,
        color='#b22222',
        linestyle='--',
        linewidth=1.5,
        alpha=0.85
    )

    ax.axhline(
        43,
        color='#b22222',
        linestyle='--',
        linewidth=1.5,
        alpha=0.85
    )

    # ------------------------------------------------------------------
    # Highlight double-expansion region
    # ------------------------------------------------------------------

    ax.axvspan(
        43,
        75,
        color='#d62728',
        alpha=0.035
    )

    ax.axhspan(
        43,
        75,
        color='#d62728',
        alpha=0.035
    )

    # ------------------------------------------------------------------
    # Threshold labels
    # ------------------------------------------------------------------

    ax.text(
        44,
        73.5,
        '≥43 repeats',
        color='#b22222',
        fontsize=10,
        fontweight='semibold',
        ha='left',
        va='top'
    )

    ax.text(
        2,
        44,
        '≥43 repeats',
        color='#b22222',
        fontsize=10,
        fontweight='semibold',
        ha='left',
        va='bottom'
    )

    # ------------------------------------------------------------------
    # Annotate double-expansion individuals
    # ------------------------------------------------------------------

    # double_df = df_gen[
    #     double_expansion
    # ]

    # for _, row in double_df.iterrows():

    #     # participant = str(
    #     #     row['Participant ID']
    #     # )

    #     lower = int(
    #         row['LOWER ALLELE']
    #     )

    #     upper = int(
    #         row['UPPER ALLELE']
    #     )

    #     ax.annotate(
    #         f'({lower}/{upper})',
    #         xy=(lower, upper),
    #         xytext=(lower - 5, upper - 3),
    #         # arrowprops=dict(
    #         #     arrowstyle='->',
    #         #     lw=1,
    #         #     color='#d62728'
    #         # ),
    #         fontsize=8.5,
    #         fontweight='semibold',
    #         color='#b22222',
    #         ha='right'
    #     )

    # ------------------------------------------------------------------
    # Labels
    # ------------------------------------------------------------------

    ax.set_title(
        'LOWER vs UPPER CAG Repeat Alleles',
        fontsize=15,
        fontweight='bold',
        pad=12
    )

    ax.set_xlabel(
        'LOWER allele (CAG repeats)',
        fontsize=12,
        labelpad=8
    )

    ax.set_ylabel(
        'UPPER allele (CAG repeats)',
        fontsize=12,
        labelpad=8
    )

    # ------------------------------------------------------------------
    # Axis limits
    # ------------------------------------------------------------------

    ax.set_xlim(0, 75)
    ax.set_ylim(40, 75)

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------

    ax.legend(
        title=None,
        frameon=False,
        fontsize=10,
        loc='lower right'
    )

    # ------------------------------------------------------------------
    # Clean appearance
    # ------------------------------------------------------------------

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(
        axis='both',
        alpha=0.20
    )

    plt.tight_layout()

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------

    os.makedirs(
        os.path.dirname(output_path),
        exist_ok=True
    )

    plt.savefig(
        output_path,
        dpi=600,
        bbox_inches='tight'
    )

    plt.close()

    print(
        f"Allele scatter plot saved to: {output_path}"
    )


if __name__ == '__main__':
    plot_allele_scatter()