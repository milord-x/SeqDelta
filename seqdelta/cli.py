from __future__ import annotations

import json
from pathlib import Path

import typer
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from .mutation import analyze_files
from .report import build_key_mutation, describe_mutation, write_csv_report, write_html_report, write_json_report

app = typer.Typer(
    add_completion=False,
    rich_markup_mode="rich",
    help="SeqDelta mutation analysis toolkit for comparing reference and mutant FASTA sequences.",
)
console = Console()

TYPE_STYLES = {
    "substitution": "yellow",
    "insertion": "blue",
    "deletion": "red",
}

EFFECT_STYLES = {
    "silent": "green",
    "missense": "yellow",
    "nonsense": "bold red",
    "frameshift": "bold magenta",
    "inframe_insertion": "blue",
    "inframe_deletion": "red",
}


def _effect_text(effect: str | None) -> str:
    return effect.replace("_", " ") if effect else "n/a"


@app.callback()
def main() -> None:
    """Run SeqDelta commands."""


def _sequence_preview(label: str, aligned_sequence: str, width: int = 96) -> Text:
    preview = Text(f"{label:<10} ", style="bold white")
    for character in aligned_sequence[:width]:
        if character == "-":
            preview.append(character, style="bold red")
        else:
            preview.append(character, style="white")
    if len(aligned_sequence) > width:
        preview.append(" ...", style="dim")
    return preview


def _diff_preview(reference_aligned: str, mutant_aligned: str, width: int = 96) -> Text:
    preview = Text("delta      ", style="bold white")
    for ref_nt, mut_nt in zip(reference_aligned[:width], mutant_aligned[:width]):
        if ref_nt == mut_nt:
            preview.append("│", style="green")
        elif "-" in (ref_nt, mut_nt):
            preview.append(" ", style="dim")
        else:
            preview.append("•", style="yellow")
    if len(reference_aligned) > width:
        preview.append(" ...", style="dim")
    return preview


def _render_summary(result) -> None:
    summary = result.summary
    table = Table.grid(expand=True)
    table.add_column()
    table.add_column()
    table.add_column()
    table.add_column()
    table.add_row(
        f"[bold]Total mutations[/bold]\n{summary.total_mutations}",
        f"[bold]Substitutions[/bold]\n{summary.substitutions}",
        f"[bold]Insertions[/bold]\n{summary.insertions}",
        f"[bold]Deletions[/bold]\n{summary.deletions}",
    )
    table.add_row(
        f"[bold]Silent[/bold]\n{summary.silent}",
        f"[bold]Missense[/bold]\n{summary.missense}",
        f"[bold]Nonsense[/bold]\n{summary.nonsense}",
        f"[bold]Frameshift[/bold]\n{summary.frameshift}",
    )
    console.print(
        Panel(
            table,
            title=f"SeqDelta Summary  |  identity {summary.identity:.1%}",
            border_style="cyan",
        )
    )


def _render_warnings(result) -> None:
    if not result.warnings:
        return
    body = "\n".join(f"- {warning}" for warning in result.warnings)
    console.print(Panel(body, title="Analysis Warnings", border_style="magenta"))


def _render_input_context(reference_fasta: Path, mutant_fasta: Path, result) -> None:
    table = Table.grid(padding=(0, 2))
    table.add_column(style="bold white")
    table.add_column()
    table.add_row("Reference file", reference_fasta.name)
    table.add_row("Reference header", result.reference.description)
    table.add_row("Mutant file", mutant_fasta.name)
    table.add_row("Mutant header", result.mutant.description)
    table.add_row("Lengths", f"{result.reference.length} nt vs {result.mutant.length} nt")
    table.add_row("Detected mutations", str(result.summary.total_mutations))
    console.print(Panel(table, title="Input Context", border_style="blue"))


def _render_key_mutation(result) -> None:
    key_mutation = build_key_mutation(result)
    if not key_mutation:
        return
    mutation = key_mutation["mutation"]
    table = Table.grid(padding=(0, 2))
    table.add_column(style="bold white")
    table.add_column()
    table.add_row("Position", str(mutation.position))
    table.add_row("Nucleotide change", f"{mutation.ref_nt} -> {mutation.mut_nt}")
    table.add_row("Codon change", f"{mutation.ref_codon or '-'} -> {mutation.mut_codon or '-'}")
    table.add_row("Amino acid change", f"{mutation.ref_aa or '-'} -> {mutation.mut_aa or '-'}")
    table.add_row("Effect", _effect_text(mutation.effect))
    table.add_row("Interpretation", describe_mutation(mutation))
    console.print(Panel(table, title="Key Mutation", border_style="yellow"))


def _render_mutation_table(result) -> None:
    table = Table(title="Nucleotide-Level Mutation Table", box=box.SIMPLE_HEAVY, header_style="bold cyan")
    table.add_column("ID", style="bold")
    table.add_column("Pos", justify="right")
    table.add_column("Type")
    table.add_column("Ref")
    table.add_column("Mut")
    table.add_column("Ref codon")
    table.add_column("Mut codon")
    table.add_column("AA")
    table.add_column("Effect")

    for mutation in result.nucleotide_mutations:
        effect_style = EFFECT_STYLES.get(mutation.effect or "", "white")
        mutation_style = TYPE_STYLES.get(mutation.mutation_type, "white")
        aa_pair = f"{mutation.ref_aa or '-'}→{mutation.mut_aa or '-'}"
        table.add_row(
            mutation.mutation_id,
            str(mutation.position),
            f"[{mutation_style}]{mutation.mutation_type}[/{mutation_style}]",
            mutation.ref_nt,
            mutation.mut_nt,
            mutation.ref_codon or "-",
            mutation.mut_codon or "-",
            aa_pair,
            f"[{effect_style}]{_effect_text(mutation.effect)}[/{effect_style}]",
        )

    if not result.nucleotide_mutations:
        table.add_row("-", "-", "no change", "-", "-", "-", "-", "-", "unchanged")

    console.print(table)


def _render_protein_tables(result) -> None:
    codon_table = Table(title="Codon-Level Consequences", box=box.MINIMAL_DOUBLE_HEAD, header_style="bold green")
    codon_table.add_column("Codon")
    codon_table.add_column("Ref")
    codon_table.add_column("Mut")
    codon_table.add_column("AA")
    codon_table.add_column("Effect")

    changed_codons = [codon for codon in result.codon_changes if codon.changed][:20]
    for codon in changed_codons:
        style = EFFECT_STYLES.get(codon.effect, "white")
        codon_table.add_row(
            str(codon.codon_index),
            codon.ref_codon or "-",
            codon.mut_codon or "-",
            f"{codon.ref_aa or '-'}→{codon.mut_aa or '-'}",
            f"[{style}]{_effect_text(codon.effect)}[/{style}]",
        )

    if not changed_codons:
        codon_table.add_row("-", "-", "-", "-", "unchanged")

    aa_table = Table(title="Protein-Level Consequences", box=box.MINIMAL_DOUBLE_HEAD, header_style="bold magenta")
    aa_table.add_column("AA pos")
    aa_table.add_column("Ref AA")
    aa_table.add_column("Mut AA")
    aa_table.add_column("Effect")

    changed_amino_acids = [aa for aa in result.amino_acid_changes if aa.changed][:20]
    for amino_acid in changed_amino_acids:
        style = EFFECT_STYLES.get(amino_acid.effect, "white")
        aa_table.add_row(
            str(amino_acid.aa_index),
            amino_acid.ref_aa or "-",
            amino_acid.mut_aa or "-",
            f"[{style}]{_effect_text(amino_acid.effect)}[/{style}]",
        )

    if not changed_amino_acids:
        aa_table.add_row("-", "-", "-", "unchanged")

    console.print(codon_table)
    console.print(aa_table)


@app.command()
def compare(
    reference_fasta: Path = typer.Argument(..., exists=True, readable=True, help="Reference FASTA file."),
    mutant_fasta: Path = typer.Argument(..., exists=True, readable=True, help="Mutant FASTA file."),
    pretty: bool = typer.Option(True, "--pretty/--plain", help="Render polished Rich terminal output."),
    json_out: Path | None = typer.Option(None, "--json-out", help="Write the full structured analysis result as JSON."),
    csv_out: Path | None = typer.Option(None, "--csv-out", help="Write the nucleotide mutation table as CSV."),
    html_report: Path | None = typer.Option(None, "--html-report", help="Write a standalone HTML report with figures and tables."),
    protein_view: bool = typer.Option(False, "--protein-view/--no-protein-view", help="Include codon-level and protein-level tables in terminal output."),
) -> None:
    """Compare a reference FASTA and a mutant FASTA."""

    result = analyze_files(reference_fasta, mutant_fasta)

    if pretty:
        console.print(
            Panel.fit(
                f"[bold cyan]SeqDelta[/bold cyan]\n[white]{reference_fasta.name}[/white]  vs  [white]{mutant_fasta.name}[/white]",
                border_style="cyan",
            )
        )
        _render_input_context(reference_fasta, mutant_fasta, result)
        _render_warnings(result)
        _render_summary(result)
        _render_key_mutation(result)
        console.print(_sequence_preview("reference", result.alignment.reference_aligned))
        console.print(_diff_preview(result.alignment.reference_aligned, result.alignment.mutant_aligned))
        console.print(_sequence_preview("mutant", result.alignment.mutant_aligned))
        _render_mutation_table(result)
        if protein_view:
            _render_protein_tables(result)
    else:
        console.print(json.dumps(result.to_dict(), indent=2))

    if json_out:
        write_json_report(result, json_out)
        console.print(f"[green]JSON report saved:[/green] {json_out}")

    if csv_out:
        write_csv_report(result, csv_out)
        console.print(f"[green]CSV report saved:[/green] {csv_out}")

    if html_report:
        write_html_report(result, html_report)
        console.print(f"[green]HTML report saved:[/green] {html_report}")


if __name__ == "__main__":
    app()
