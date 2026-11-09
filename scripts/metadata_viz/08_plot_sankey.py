import os
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from matplotlib.colors import to_rgba

# Ensure scripts directory is in python search path
scripts_dir = os.path.dirname(os.path.abspath(__file__))
if scripts_dir not in sys.path:
    sys.path.append(scripts_dir)

import importlib
plot_proportions = importlib.import_module("07_plot_proportions")
load_clean_data = plot_proportions.load_clean_data


def generate_sankey_data(df):
    """
    Prepares node and link information for the Sankey plot.
    """

    df_sankey = pd.DataFrame()

    # ---------------------------------------------------------
    # FIRST SYMPTOM
    # ---------------------------------------------------------
    df_sankey['First Symptom'] = df['FIRST_SYMPTOM'].replace({
        'Hand Tremor': 'Sx: Hand Tremor',
        'Gait Ataxia': 'Sx: Gait Ataxia',
        'Tremor And Gait Ataxia': 'Sx: Mixed (Tremor/Ataxia)',
        'Titubation': 'Sx: Titubation',
        'Tremulous Voice': 'Sx: Tremulous Voice'
    })

    # ---------------------------------------------------------
    # MRI
    # ---------------------------------------------------------
    df_sankey['MRI'] = 'MRI: ' + df['MRI_CATEGORIES'].astype(str)

    # ---------------------------------------------------------
    # PARASYMPATHETIC
    # ---------------------------------------------------------
    df_sankey['Parasympathetic'] = (
        df['PARASYMPATHETIC_ACTIVITY_CATEGORIES']
        .replace({
            'Normal': 'Para: Normal',
            'Mild Loss': 'Para: Mild/Mod Loss',
            'Moderate Loss': 'Para: Mild/Mod Loss',
            'Severe Loss': 'Para: Severe Loss',
            'Unknown': 'Para: Unknown'
        })
    )

    # ---------------------------------------------------------
    # SYMPATHETIC
    # ---------------------------------------------------------
    df_sankey['Sympathetic'] = (
        df['SYMPATHETIC_ACTIVITY_CATEGORIES']
        .replace({
            'Normal': 'Symp: Normal',
            'Mild Loss': 'Symp: Mild/Mod Loss',
            'Moderate Loss': 'Symp: Mild/Mod Loss',
            'Severe Loss': 'Symp: Severe Loss',
            'Unknown': 'Symp: Unknown'
        })
    )

    # ---------------------------------------------------------
    # BRADYKINESIA
    # ---------------------------------------------------------
    df_sankey['Bradykinesia'] = (
        df['BRADYKINESIA_CATEGORIES']
        .replace({
            'Present': 'Bradykinesia: Present',
            'Absent': 'Bradykinesia: Absent',
            'Unknown': 'Bradykinesia: Unknown'
        })
    )

    # ---------------------------------------------------------
    # PLANTAR
    # ---------------------------------------------------------
    df_sankey['Plantar'] = (
        df['PLANTAR_CATEGORIES']
        .replace({
            'Flexor (Normal)': 'Plantar: Flexor (Normal)',
            'Extensor (Abnormal)': 'Plantar: Extensor (Abnormal)',
            'Equivocal': 'Plantar: Equivocal',
            'Unknown': 'Plantar: Unknown'
        })
    )

    cols = [
        'First Symptom',
        'MRI',
        'Parasympathetic',
        'Sympathetic',
        'Bradykinesia',
        'Plantar'
    ]

    # ---------------------------------------------------------
    # CREATE UNIQUE NODE LIST
    # ---------------------------------------------------------
    all_nodes = []

    for col in cols:
        unique_vals = sorted(
            df_sankey[col].dropna().unique().tolist()
        )
        all_nodes.extend(unique_vals)

    node_to_idx = {
        node: idx for idx, node in enumerate(all_nodes)
    }

    # ---------------------------------------------------------
    # CREATE LINKS
    # ---------------------------------------------------------
    sources = []
    targets = []
    values = []

    for i in range(len(cols) - 1):

        col_src = cols[i]
        col_tgt = cols[i + 1]

        grouped = (
            df_sankey
            .groupby([col_src, col_tgt])
            .size()
            .reset_index(name='count')
        )

        for _, row in grouped.iterrows():

            sources.append(
                node_to_idx[row[col_src]]
            )

            targets.append(
                node_to_idx[row[col_tgt]]
            )

            values.append(
                row['count']
            )

    return all_nodes, sources, targets, values, df_sankey


# =============================================================
# COLOUR UTILITIES
# =============================================================

def hex_to_rgba(hex_color, alpha=0.38):
    """
    Convert #RRGGBB to rgba(...) for Plotly.
    """
    hex_color = hex_color.lstrip('#')

    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)

    return f'rgba({r}, {g}, {b}, {alpha})'


# =============================================================
# MAIN PLOT
# =============================================================

def plot_sankey_diagram(output_html=None, output_png=None):

    if output_html is None:
        output_html = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'plots',
            '09_sankey_diagram.html'
        )

    if output_png is None:
        output_png = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'plots',
            '09_sankey_diagram.png'
        )

    # =========================================================
    # LOAD DATA
    # =========================================================

    df = load_clean_data()

    nodes, sources, targets, values, df_sankey = generate_sankey_data(df)

    # =========================================================
    # COLOUR PALETTES
    # =========================================================

    # First symptom colours
    symptom_colors = {
        'Sx: Hand Tremor': '#2F6690',
        'Sx: Gait Ataxia': '#D95F02',
        'Sx: Mixed (Tremor/Ataxia)': '#7B3294',
        'Sx: Titubation': '#1B9E77',
        'Sx: Tremulous Voice': '#E7298A'
    }

    # MRI
    mri_colors = {
        'Normal': '#2ca02c',
        'Cerebellar Atrophy': '#d62728',
        'Cerebral Atrophy': '#ff7f0e',
        'Cerebral & Cerebellar Atrophy': '#9467bd',
        'Others': '#bcbd22',
        'Unknown': '#7f7f7f'
    }

    # Autonomic
    autonomic_colors = {
        'Normal': '#66A61E',
        'Mild/Mod Loss': '#E6AB02',
        'Severe Loss': '#D95F02',
        'Unknown': '#BDBDBD'
    }

    # Bradykinesia
    brady_colors = {
        'Present': '#C51B7D',
        'Absent': '#5E3C99',
        'Unknown': '#BDBDBD'
    }

    # Plantar
    plantar_colors = {
        'Flexor (Normal)': '#1B9E77',
        'Extensor (Abnormal)': '#D73027',
        'Equivocal': '#E6AB02',
        'Unknown': '#BDBDBD'
    }

    # =========================================================
    # NODE COLOURS
    # =========================================================

    node_colors = []

    for node in nodes:

        if node.startswith('Sx:'):

            color = symptom_colors.get(
                node,
                '#7F8C8D'
            )

        elif node.startswith('MRI:'):

            clean = node.replace('MRI: ', '')

            color = mri_colors.get(
                clean,
                '#BDBDBD'
            )

        elif node.startswith('Para:'):

            clean = node.replace('Para: ', '')

            color = autonomic_colors.get(
                clean,
                '#BDBDBD'
            )

        elif node.startswith('Symp:'):

            clean = node.replace('Symp: ', '')

            color = autonomic_colors.get(
                clean,
                '#BDBDBD'
            )

        elif node.startswith('Bradykinesia:'):

            clean = node.replace(
                'Bradykinesia: ',
                ''
            )

            color = brady_colors.get(
                clean,
                '#BDBDBD'
            )

        elif node.startswith('Plantar:'):

            clean = node.replace(
                'Plantar: ',
                ''
            )

            color = plantar_colors.get(
                clean,
                '#BDBDBD'
            )

        else:

            color = '#BDBDBD'

        node_colors.append(color)

    # =========================================================
    # CLEAN NODE LABELS
    # =========================================================

    clean_node_labels = []

    for node in nodes:

        clean_name = node

        prefixes = [
            'Sx: ',
            'MRI: ',
            'Para: ',
            'Symp: ',
            'Bradykinesia: ',
            'Plantar: '
        ]

        for prefix in prefixes:

            if clean_name.startswith(prefix):

                clean_name = clean_name[
                    len(prefix):
                ]

                break

        clean_node_labels.append(clean_name)

    # =========================================================
    # DEFINE GROUPS
    # =========================================================

    groups = [
        ('First Symptom', 'Sx:'),
        ('MRI Category', 'MRI:'),
        ('Parasympathetic Activity', 'Para:'),
        ('Sympathetic Activity', 'Symp:'),
        ('Bradykinesia', 'Bradykinesia:'),
        ('Plantar Response', 'Plantar:')
    ]

    # =========================================================
    # NODE POSITIONS
    # =========================================================

    # Horizontal position of each clinical domain
    group_x = {
        'First Symptom': 0.035,
        'MRI Category': 0.220,
        'Parasympathetic Activity': 0.405,
        'Sympathetic Activity': 0.590,
        'Bradykinesia': 0.775,
        'Plantar Response': 0.960
    }

    node_positions = {}
    y_min = 0.05
    y_max = 0.95
    padding = 0.04  # Vertical padding between nodes

    # Map group prefixes to columns in df_sankey
    prefix_to_col = {
        'Sx:': 'First Symptom',
        'MRI:': 'MRI',
        'Para:': 'Parasympathetic',
        'Symp:': 'Sympathetic',
        'Bradykinesia:': 'Bradykinesia',
        'Plantar:': 'Plantar'
    }

    for group_name, prefix in groups:
        group_nodes = [node for node in nodes if node.startswith(prefix)]
        group_nodes = sorted(group_nodes)

        col_name = prefix_to_col[prefix]
        counts = df_sankey[col_name].value_counts()
        counts = counts.reindex(group_nodes).fillna(0)

        k = len(group_nodes)
        s_total = counts.sum()

        H = y_max - y_min
        total_padding = (k - 1) * padding
        H_nodes = H - total_padding

        if H_nodes < 0.2:
            current_padding = max((H - 0.2) / (k - 1), 0.005)
            H_nodes = H - (k - 1) * current_padding
        else:
            current_padding = padding

        current_y = y_min
        for node in group_nodes:
            size = counts[node]
            h_i = (size / s_total) * H_nodes if s_total > 0 else (H_nodes / k)
            center_y = current_y + h_i / 2
            node_positions[node] = (group_x[group_name], center_y)
            current_y += h_i + current_padding

    node_x = [node_positions[node][0] for node in nodes]
    node_y = [node_positions[node][1] for node in nodes]


    # =========================================================
    # LINK COLOURS
    # =========================================================

    def hex_to_rgba(hex_color, alpha=0.28):

        hex_color = hex_color.lstrip('#')

        r = int(hex_color[0:2], 16)
        g = int(hex_color[2:4], 16)
        b = int(hex_color[4:6], 16)

        return (
            f'rgba({r}, {g}, {b}, {alpha})'
        )

    link_colors = []

    for source_idx in sources:

        link_colors.append(
            hex_to_rgba(
                node_colors[source_idx],
                alpha=0.30
            )
        )


    # =========================================================
    # BUILD SANKEY
    # =========================================================

    fig = go.Figure(
        data=[
            go.Sankey(
                
                arrangement='fixed',
                
                # ADD THIS: Make the default labels completely transparent
                textfont=dict(color='rgba(0,0,0,0)', size=1), 

                node=dict(

                    x=node_x,
                    y=node_y,

                    pad=38,
                    thickness=24,

                    line=dict(
                        color='rgba(0,0,0,0.25)',
                        width=0.7
                    ),

                    # FIX: Pass the unique nodes list to prevent D3 mapping bugs
                    label=[""] * len(nodes),

                    color=node_colors,

                    customdata=clean_node_labels,

                    hovertemplate=(
                        '<b>%{customdata}</b>'
                        '<br>'
                        'Patients: %{value}'
                        '<extra></extra>'
                    )
                ),

                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors,
                    hovertemplate=(
                        '<b>%{source.customdata}</b>'
                        ' → '
                        '<b>%{target.customdata}</b>'
                        '<br>'
                        'Patients: %{value}'
                        '<extra></extra>'
                    )
                )
            )
        ]
    )


    # =========================================================
    # ANNOTATIONS
    # =========================================================

    annotations = []


    # =========================================================
    # OVERALL TITLE
    # =========================================================

    annotations.append(
        dict(

            x=0.5,
            y=1.16,

            xref='paper',
            yref='paper',

            text=(
                '<b>Clinical Phenotypic Flow Across SCA12 Patients</b>'
            ),

            showarrow=False,

            xanchor='center',
            yanchor='middle',

            font=dict(
                family='Arial',
                size=22,
                color='#1F1F1F'
            )
        )
    )


    # =========================================================
    # GROUP / DOMAIN HEADINGS
    # =========================================================

    for group_name, prefix in groups:

        x = group_x[group_name]

        # Slightly below title
        y = 1.065

        annotations.append(
            dict(

                x=x,
                y=y,

                xref='paper',
                yref='paper',

                text=f'<b>{group_name}</b>',

                showarrow=False,

                xanchor='center',
                yanchor='middle',

                font=dict(
                    family='Arial',
                    size=15,
                    color='#333333'
                ),

                align='center'
            )
        )


    # =========================================================
    # CATEGORY LABELS
    # =========================================================
    #
    # Labels are positioned BELOW each bar.
    #
    # They use the exact same x/y coordinate as their node,
    # with a controlled offset so that every column has the
    # same visual treatment.
    # =========================================================

    for i, node in enumerate(nodes):

        x = node_x[i]
        y = node_y[i]

        label = clean_node_labels[i]

        # -----------------------------------------------------
        # Special formatting for long labels
        # -----------------------------------------------------

        if label == 'Mixed (Tremor/Ataxia)':

            label = (
                'Mixed<br>'
                '(Tremor/Ataxia)'
            )

        elif label == 'Extensor (Abnormal)':

            label = (
                'Extensor<br>'
                '(Abnormal)'
            )

        elif label == 'Flexor (Normal)':

            label = (
                'Flexor<br>'
                '(Normal)'
            )

        elif label == 'Mild/Mod Loss':

            label = (
                'Mild/Moderate<br>'
                'Loss'
            )

        # -----------------------------------------------------
        # Put label underneath node
        # -----------------------------------------------------

        text_color = '#000000'

        annotations.append(
            dict(

                x=x,
                y=1.0 - y,

                xref='paper',
                yref='paper',

                text=f"<b>{label}</b>",

                showarrow=False,

                xanchor='center',
                yanchor='middle',

                font=dict(
                    family='Arial',
                    size=14,
                    color=text_color
                ),

                align='center'
            )
        )


    # =========================================================
    # LAYOUT
    # =========================================================

    fig.update_layout(

        width=1500,
        height=900,

        paper_bgcolor='white',
        plot_bgcolor='white',

        font=dict(
            family='Arial',
            size=13,
            color='#222222'
        ),

        margin=dict(
            t=145,
            l=45,
            r=45,
            b=110
        ),

        annotations=annotations
    )

    # =========================================================
    # SAVE HTML
    # =========================================================

    os.makedirs(
        os.path.dirname(output_html),
        exist_ok=True
    )

    fig.write_html(
        output_html,
        include_plotlyjs=True
    )

    print(
        f"Interactive Sankey plot saved to: "
        f"{output_html}"
    )

    # =========================================================
    # SAVE PNG
    # =========================================================

    try:

        fig.write_image(
            output_png,
            width=1450,
            height=850,
            scale=3
        )

        print(
            f"Static PNG saved to: "
            f"{output_png}"
        )

    except Exception as e:

        print(
            "\nERROR: PNG export failed."
        )

        print(
            f"Reason: {e}"
        )

        print(
            "\nThe interactive HTML was successfully created."
        )


if __name__ == '__main__':
    plot_sankey_diagram()