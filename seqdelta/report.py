from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from .models import AnalysisResult
from .visualization import build_visualizations

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TEMPLATES_DIR = PROJECT_ROOT / "templates"
STATIC_DIR = PROJECT_ROOT / "static"


def _resolve_resource_dir(preferred: Path, package_relative: str) -> Path:
    if preferred.exists():
        return preferred
    return Path(str(files("seqdelta").joinpath(package_relative)))


RESOLVED_TEMPLATES_DIR = _resolve_resource_dir(TEMPLATES_DIR, "assets/templates")
RESOLVED_STATIC_DIR = _resolve_resource_dir(STATIC_DIR, "assets/static")

environment = Environment(
    loader=FileSystemLoader(RESOLVED_TEMPLATES_DIR),
    autoescape=select_autoescape(["html", "xml"]),
)

AA_NAMES = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "Q": "Glutamine",
    "E": "Glutamic acid",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "*": "Stop codon",
    "X": "Unknown residue",
}

EFFECT_PRIORITY = {
    "frameshift": 0,
    "nonsense": 1,
    "missense": 2,
    "inframe_deletion": 3,
    "inframe_insertion": 4,
    "silent": 5,
    None: 6,
}


def _column_class(ref_nt: str, mut_nt: str) -> str:
    if ref_nt == mut_nt:
        return "match"
    if ref_nt == "-":
        return "insertion"
    if mut_nt == "-":
        return "deletion"
    return "substitution"


def _aa_label(amino_acid: str | None) -> str:
    if amino_acid is None:
        return "unknown residue"
    return f"{AA_NAMES.get(amino_acid, 'Unknown residue')} ({amino_acid})"


def _format_effect(effect: str | None) -> str:
    return (effect or "not classified").replace("_", " ").title()


def describe_mutation(mutation) -> str:
    if mutation.effect == "silent":
        return (
            f"Codon changed from {mutation.ref_codon or '-'} to {mutation.mut_codon or '-'}, "
            f"but both still translate to {_aa_label(mutation.ref_aa)}. Protein sequence is preserved at this site."
        )
    if mutation.effect == "missense":
        return (
            f"Amino acid changed from {_aa_label(mutation.ref_aa)} to {_aa_label(mutation.mut_aa)}, "
            "which may alter local protein structure or function."
        )
    if mutation.effect == "nonsense":
        return (
            f"The codon changed from {mutation.ref_codon or '-'} to {mutation.mut_codon or '-'} and now encodes "
            f"{_aa_label(mutation.mut_aa)}. This introduces a premature stop signal that can truncate the protein."
        )
    if mutation.effect == "frameshift":
        return (
            "This insertion/deletion shifts the reading frame. Downstream codons are redefined, which can trigger "
            "widespread amino acid changes and premature stop codons."
        )
    if mutation.effect == "inframe_insertion":
        return (
            "Inserted bases preserve the reading frame because the event size is divisible by three. "
            "The protein gains residues without shifting downstream codons."
        )
    if mutation.effect == "inframe_deletion":
        return (
            "Deleted bases preserve the reading frame because the event size is divisible by three. "
            "The protein loses residues without shifting downstream codons."
        )
    return mutation.note or "Sequence difference detected."


def _mutation_sort_key(mutation) -> tuple[int, int]:
    return EFFECT_PRIORITY.get(mutation.effect, 99), mutation.position


def _highlight_segments(sequence: str, start: int | None, end: int | None, highlight: str, window: int) -> dict[str, str]:
    if not sequence:
        return {"leading": "", "prefix": "", "focus": highlight or "-", "suffix": "", "trailing": ""}

    if start is None:
        anchor = max((end or 1) - 1, 0)
        left = max(0, anchor - window)
        right = min(len(sequence), anchor + window)
        return {
            "leading": "..." if left > 0 else "",
            "prefix": sequence[left:anchor],
            "focus": highlight or "-",
            "suffix": sequence[anchor:right],
            "trailing": "..." if right < len(sequence) else "",
        }

    zero_based_start = max(start - 1, 0)
    zero_based_end = max(end or start, start)
    left = max(0, zero_based_start - window)
    right = min(len(sequence), zero_based_end + window)
    return {
        "leading": "..." if left > 0 else "",
        "prefix": sequence[left:zero_based_start],
        "focus": highlight or "-",
        "suffix": sequence[zero_based_end:right],
        "trailing": "..." if right < len(sequence) else "",
    }


def build_mutation_focus_views(result: AnalysisResult, window: int = 12) -> list[dict]:
    views = []
    for mutation in sorted(result.nucleotide_mutations, key=lambda item: item.position):
        if mutation.mutation_type == "insertion":
            ref_segments = _highlight_segments(
                result.reference.sequence,
                None,
                mutation.ref_start,
                "-",
                window,
            )
            mut_segments = _highlight_segments(
                result.mutant.sequence,
                mutation.mut_start,
                mutation.mut_end,
                mutation.mut_nt,
                window,
            )
        elif mutation.mutation_type == "deletion":
            ref_segments = _highlight_segments(
                result.reference.sequence,
                mutation.ref_start,
                mutation.ref_end,
                mutation.ref_nt,
                window,
            )
            mut_segments = _highlight_segments(
                result.mutant.sequence,
                None,
                mutation.mut_start,
                "-",
                window,
            )
        else:
            ref_segments = _highlight_segments(
                result.reference.sequence,
                mutation.ref_start,
                mutation.ref_end,
                mutation.ref_nt,
                window,
            )
            mut_segments = _highlight_segments(
                result.mutant.sequence,
                mutation.mut_start,
                mutation.mut_end,
                mutation.mut_nt,
                window,
            )

        views.append(
            {
                "mutation": mutation,
                "reference": ref_segments,
                "mutant": mut_segments,
                "effect_label": _format_effect(mutation.effect),
                "effect_explanation": describe_mutation(mutation),
            }
        )
    return views


def build_input_summary(result: AnalysisResult) -> dict:
    return {
        "reference_file": result.reference.source_name or "reference.fasta",
        "mutant_file": result.mutant.source_name or "mutant.fasta",
        "reference_header": result.reference.description,
        "mutant_header": result.mutant.description,
        "reference_length": result.reference.length,
        "mutant_length": result.mutant.length,
        "mutation_count": len(result.nucleotide_mutations),
    }


def build_key_mutation(result: AnalysisResult) -> dict | None:
    if not result.nucleotide_mutations:
        return None
    mutation = sorted(result.nucleotide_mutations, key=_mutation_sort_key)[0]
    return {
        "mutation": mutation,
        "effect_label": _format_effect(mutation.effect),
        "effect_explanation": describe_mutation(mutation),
    }


def build_alignment_blocks(result: AnalysisResult, block_size: int = 60) -> list[dict]:
    reference_aligned = result.alignment.reference_aligned
    mutant_aligned = result.alignment.mutant_aligned
    blocks: list[dict] = []

    ref_cursor = 0
    mut_cursor = 0
    for start in range(0, len(reference_aligned), block_size):
        ref_start = ref_cursor + 1
        mut_start = mut_cursor + 1
        cells = []
        reference_slice = reference_aligned[start : start + block_size]
        mutant_slice = mutant_aligned[start : start + block_size]

        for ref_nt, mut_nt in zip(reference_slice, mutant_slice):
            cells.append(
                {
                    "ref_nt": ref_nt,
                    "mut_nt": mut_nt,
                    "marker": "|" if ref_nt == mut_nt and ref_nt != "-" else "•" if "-" not in (ref_nt, mut_nt) else " ",
                    "class_name": _column_class(ref_nt, mut_nt),
                }
            )
            if ref_nt != "-":
                ref_cursor += 1
            if mut_nt != "-":
                mut_cursor += 1

        blocks.append(
            {
                "ref_start": ref_start,
                "ref_end": ref_cursor,
                "mut_start": mut_start,
                "mut_end": mut_cursor,
                "cells": cells,
            }
        )

    return blocks


def _effect_interpretation(result: AnalysisResult) -> str:
    summary = result.summary
    if summary is None:
        return "No summary available."
    if summary.frameshift:
        return "Frameshift signal detected: downstream codons and amino acids are likely broadly altered."
    if summary.nonsense:
        return "Stop-gain signal detected: at least one mutation introduces a premature stop codon."
    if summary.missense:
        return "Protein-changing signal detected: one or more codons alter the translated amino acid."
    if summary.silent:
        return "Sequence differences are present, but at least one detected substitution is silent at the amino acid level."
    return "No sequence differences were detected between the two inputs."


def build_report_context(
    result: AnalysisResult,
    *,
    standalone_report: bool,
    report_title: str = "SeqDelta Mutation Report",
    include_plotly_bundle: bool = True,
    download_links: dict[str, str] | None = None,
) -> dict:
    changed_codons = [codon for codon in result.codon_changes if codon.changed]
    changed_amino_acids = [aa for aa in result.amino_acid_changes if aa.changed]
    visuals = build_visualizations(result, include_plotly_bundle=include_plotly_bundle)
    mutation_focus_views = build_mutation_focus_views(result)

    return {
        "report_title": report_title,
        "result": result,
        "summary": result.summary,
        "warnings": result.warnings,
        "input_summary": build_input_summary(result),
        "key_mutation": build_key_mutation(result),
        "alignment_blocks": build_alignment_blocks(result),
        "mutation_focus_views": mutation_focus_views,
        "visuals": visuals,
        "changed_codons": changed_codons,
        "changed_amino_acids": changed_amino_acids,
        "effect_interpretation": _effect_interpretation(result),
        "standalone_report": standalone_report,
        "download_links": download_links or {},
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "inline_css": (RESOLVED_STATIC_DIR / "style.css").read_text(encoding="utf-8"),
        "inline_js": (RESOLVED_STATIC_DIR / "app.js").read_text(encoding="utf-8"),
    }


def render_html_report(result: AnalysisResult, report_title: str = "SeqDelta Mutation Report") -> str:
    template = environment.get_template("results.html")
    context = build_report_context(
        result,
        standalone_report=True,
        report_title=report_title,
        include_plotly_bundle=True,
    )
    return template.render(**context)


def write_html_report(result: AnalysisResult, output_path: str | Path) -> Path:
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(render_html_report(result), encoding="utf-8")
    return target


def write_json_report(result: AnalysisResult, output_path: str | Path) -> Path:
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return target


def write_csv_report(result: AnalysisResult, output_path: str | Path) -> Path:
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mutation_id",
                "position",
                "ref_nt",
                "mut_nt",
                "mutation_type",
                "ref_codon",
                "mut_codon",
                "ref_aa",
                "mut_aa",
                "effect",
                "note",
            ],
        )
        writer.writeheader()
        for mutation in result.nucleotide_mutations:
            writer.writerow(
                {
                    "mutation_id": mutation.mutation_id,
                    "position": mutation.position,
                    "ref_nt": mutation.ref_nt,
                    "mut_nt": mutation.mut_nt,
                    "mutation_type": mutation.mutation_type,
                    "ref_codon": mutation.ref_codon or "",
                    "mut_codon": mutation.mut_codon or "",
                    "ref_aa": mutation.ref_aa or "",
                    "mut_aa": mutation.mut_aa or "",
                    "effect": mutation.effect or "",
                    "note": mutation.note or "",
                }
            )
    return target
