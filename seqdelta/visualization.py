from __future__ import annotations

from collections import Counter
from math import ceil

import plotly.graph_objects as go
from plotly.offline import plot

from .models import AnalysisResult

TYPE_COLORS = {
    "substitution": "#f59e0b",
    "insertion": "#2563eb",
    "deletion": "#dc2626",
}

EFFECT_COLORS = {
    "silent": "#16a34a",
    "missense": "#f97316",
    "nonsense": "#7f1d1d",
    "frameshift": "#7c3aed",
    "inframe_insertion": "#2563eb",
    "inframe_deletion": "#dc2626",
    "unchanged": "#64748b",
}


def _to_div(figure: go.Figure, include_plotlyjs: str | bool = False) -> str:
    figure.update_layout(
        template="plotly_white",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#f8fafc",
        margin=dict(l=40, r=20, t=50, b=40),
        font=dict(family="IBM Plex Sans, Segoe UI, sans-serif", color="#0f172a"),
    )
    return plot(
        figure,
        include_plotlyjs=include_plotlyjs,
        output_type="div",
        config={"displayModeBar": False, "responsive": True},
    )


def mutation_map_div(result: AnalysisResult, include_plotlyjs: str | bool = False) -> str:
    mutations = result.nucleotide_mutations
    if not mutations:
        figure = go.Figure()
        figure.add_annotation(text="No mutations detected", x=0.5, y=0.5, showarrow=False)
        figure.update_xaxes(visible=False)
        figure.update_yaxes(visible=False)
        figure.update_layout(height=240, title="Linear Mutation Map")
        return _to_div(figure, include_plotlyjs)

    y_mapping = {"substitution": 3, "insertion": 2, "deletion": 1}
    figure = go.Figure()

    for mutation_type, lane in y_mapping.items():
        subset = [mutation for mutation in mutations if mutation.mutation_type == mutation_type]
        if not subset:
            continue
        figure.add_trace(
            go.Scatter(
                x=[mutation.position for mutation in subset],
                y=[lane] * len(subset),
                mode="markers",
                marker=dict(color=TYPE_COLORS[mutation_type], size=16, line=dict(width=2, color="#ffffff")),
                name=mutation_type.title(),
                hovertemplate=(
                    "Position %{x}<br>"
                    + "Type: "
                    + mutation_type
                    + "<br>Effect: %{customdata[0]}<br>Ref: %{customdata[1]}<br>Mut: %{customdata[2]}<extra></extra>"
                ),
                customdata=[
                    [mutation.effect or "n/a", mutation.ref_nt, mutation.mut_nt]
                    for mutation in subset
                ],
            )
        )

    figure.update_layout(
        title="Linear Mutation Map",
        height=320,
        xaxis_title="Reference position",
        yaxis=dict(
            title="Mutation class",
            tickmode="array",
            tickvals=[1, 2, 3],
            ticktext=["Deletion", "Insertion", "Substitution"],
            range=[0.5, 3.5],
        ),
        legend_title="Nucleotide layer",
    )
    return _to_div(figure, include_plotlyjs)


def mutation_density_div(result: AnalysisResult) -> str:
    sequence_length = max(result.reference.length, 1)
    bin_count = min(14, max(4, ceil(sequence_length / 20)))
    bin_size = max(1, ceil(sequence_length / bin_count))
    bins = list(range(1, sequence_length + bin_size, bin_size))
    counts = [0 for _ in range(len(bins) - 1)]

    for mutation in result.nucleotide_mutations:
        index = min((max(mutation.position, 1) - 1) // bin_size, len(counts) - 1)
        counts[index] += 1

    figure = go.Figure(
        data=[
            go.Bar(
                x=[f"{start}-{min(start + bin_size - 1, sequence_length)}" for start in bins[:-1]],
                y=counts,
                marker_color="#2563eb",
                hovertemplate="Window %{x}<br>Mutations %{y}<extra></extra>",
            )
        ]
    )
    figure.update_layout(
        title="Mutation Density Overview",
        height=320,
        xaxis_title="Reference windows",
        yaxis_title="Mutation count",
    )
    return _to_div(figure)


def effect_distribution_div(result: AnalysisResult) -> str:
    effects = Counter(mutation.effect or "unknown" for mutation in result.nucleotide_mutations)
    if not effects:
        effects["unchanged"] = 1

    labels = list(effects.keys())
    values = list(effects.values())
    colors = [EFFECT_COLORS.get(label, "#64748b") for label in labels]

    figure = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.58,
                marker=dict(colors=colors),
                sort=False,
                textinfo="label+percent",
                textposition="inside",
            )
        ]
    )
    figure.update_layout(title="Protein Effect Distribution", height=320)
    return _to_div(figure)


def build_visualizations(result: AnalysisResult, include_plotly_bundle: bool = True) -> dict[str, str]:
    include_plotlyjs = "inline" if include_plotly_bundle else False
    return {
        "mutation_map": mutation_map_div(result, include_plotlyjs),
        "mutation_density": mutation_density_div(result) if len(result.nucleotide_mutations) >= 3 else "",
        "effect_distribution": effect_distribution_div(result),
    }
