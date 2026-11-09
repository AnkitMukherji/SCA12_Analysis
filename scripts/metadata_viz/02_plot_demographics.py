import os
import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import load_and_clean_data


def plot_age_distribution(output_path=None):
    """
    Generates the age distribution plot showing:
    - Age at Assessment
    - Age of Onset (AO)
    - Mean values for both distributions
    """

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'plots',
            '01_age_distribution.png'
        )

    df = load_and_clean_data()

    # Set style
    sns.set_theme(style='whitegrid')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')

    # --- Age distributions ---

    sns.kdeplot(
        data=df['Age'].dropna(),
        fill=True,
        color='#2b5c8f',
        label='Age at Assessment',
        alpha=0.5,
        linewidth=2,
        ax=ax
    )

    sns.kdeplot(
        data=df['AO'].dropna(),
        fill=True,
        color='#d95f02',
        label='Age of Onset (AO)',
        alpha=0.5,
        linewidth=2,
        ax=ax
    )

    # --- Mean lines ---

    mean_age = df['Age'].mean()
    mean_ao = df['AO'].mean()

    ax.axvline(
        mean_age,
        color='#1b3a5b',
        linestyle='--',
        linewidth=1.5,
        label=f'Mean Age: {mean_age:.1f} y'
    )

    ax.axvline(
        mean_ao,
        color='#a63603',
        linestyle='--',
        linewidth=1.5,
        label=f'Mean Onset: {mean_ao:.1f} y'
    )

    # --- Formatting ---

    ax.set_title(
        'Participant Age & Age of Onset Distribution',
        fontsize=14,
        fontweight='bold',
        pad=15,
        color='#2c3e50'
    )

    ax.set_xlabel(
        'Years',
        fontsize=12,
        fontweight='semibold'
    )

    ax.set_ylabel(
        'Density',
        fontsize=12,
        fontweight='semibold'
    )

    ax.set_xlim(15, 90)

    ax.legend(
        frameon=True,
        facecolor='#ffffff',
        edgecolor='none',
        loc='upper left'
    )

    plt.tight_layout()

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    plt.savefig(
        output_path,
        dpi=600,
        bbox_inches='tight'
    )

    plt.close()

    print(f"Age distribution plot saved to: {output_path}")


if __name__ == '__main__':
    plot_age_distribution()

import os
import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import load_and_clean_data


def plot_gender_distribution(output_path=None):
    """
    Generates a donut chart showing the gender distribution
    of the cohort.
    """

    if output_path is None:
        output_path = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'plots',
            '02_gender_distribution.png'
        )

    df = load_and_clean_data()

    # Set style
    sns.set_theme(style='whitegrid')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11

    # Create figure
    fig, ax = plt.subplots(figsize=(7, 6))

    fig.patch.set_facecolor('#fafafa')
    ax.set_facecolor('#ffffff')

    # --- Gender counts ---

    gender_counts = df['Gender'].value_counts()

    genders = gender_counts.index.tolist()
    counts = gender_counts.values.tolist()

    total = sum(counts)

    percentages = [
        f'{(c / total) * 100:.1f}%'
        for c in counts
    ]

    labels = [
        f'{g}\n(n={c}, {p})'
        for g, c, p in zip(genders, counts, percentages)
    ]

    # --- Donut chart ---

    colors = ['#4a90e2', '#e26a6a']

    wedges, texts = ax.pie(
        counts,
        labels=labels,
        startangle=140,
        colors=colors,
        wedgeprops=dict(
            width=0.4,
            edgecolor='white',
            linewidth=2
        ),
        textprops=dict(
            fontsize=12,
            color='#2c3e50',
            fontweight='semibold'
        )
    )

    # --- Center text ---

    ax.text(
        0,
        0,
        f'Total\nn = {len(df)}',
        ha='center',
        va='center',
        fontsize=14,
        fontweight='bold',
        color='#2c3e50'
    )

    # --- Title ---

    ax.set_title(
        'Cohort Gender Distribution',
        fontsize=14,
        fontweight='bold',
        pad=15,
        color='#2c3e50'
    )

    plt.tight_layout()

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    plt.savefig(
        output_path,
        dpi=600,
        bbox_inches='tight'
    )

    plt.close()

    print(f"Gender distribution plot saved to: {output_path}")


if __name__ == '__main__':
    plot_gender_distribution()