from __future__ import annotations

from collections import Counter
from pathlib import Path

from .alignment import align_sequences
from .models import AnalysisResult, AnalysisSummary, NucleotideMutation, SequenceRecordData
from .parser import parse_fasta_file
from .translation import (
    build_amino_acid_comparisons,
    build_codon_comparisons,
    codon_index_from_position,
    get_codon,
    translate_codon,
)


def _event_label(prefix: str, index: int) -> str:
    return f"{prefix}-{index:04d}"


def _classify_substitution_effect(reference_sequence: str, mutant_sequence: str, position: int) -> tuple[str | None, str | None, str | None, str | None, str]:
    codon_index = codon_index_from_position(position)
    ref_codon = get_codon(reference_sequence, codon_index)
    mut_codon = get_codon(mutant_sequence, codon_index)
    ref_aa = translate_codon(ref_codon)
    mut_aa = translate_codon(mut_codon)

    if ref_aa == mut_aa:
        effect = "silent"
    elif mut_aa == "*" and ref_aa != "*":
        effect = "nonsense"
    else:
        effect = "missense"

    return ref_codon, mut_codon, ref_aa, mut_aa, effect


def _classify_indel_effect(
    reference_sequence: str,
    mutant_sequence: str,
    mutation_type: str,
    position: int,
    ref_start: int | None,
    mut_start: int | None,
    ref_nt: str,
    mut_nt: str,
) -> tuple[str | None, str | None, str | None, str | None, str]:
    delta = len(mut_nt.replace("-", "")) - len(ref_nt.replace("-", ""))
    anchor_position = ref_start or position or 1
    codon_index = codon_index_from_position(anchor_position)
    ref_codon = get_codon(reference_sequence, codon_index)
    mut_codon = get_codon(mutant_sequence, codon_index_from_position(mut_start or anchor_position))
    ref_aa = translate_codon(ref_codon)
    mut_aa = translate_codon(mut_codon)

    if abs(delta) % 3 != 0:
        effect = "frameshift"
    else:
        effect = f"inframe_{mutation_type}"

    return ref_codon, mut_codon, ref_aa, mut_aa, effect


def _detect_nucleotide_mutations(alignment, reference_sequence: str, mutant_sequence: str) -> list[NucleotideMutation]:
    mutations: list[NucleotideMutation] = []
    reference_aligned = alignment.reference_aligned
    mutant_aligned = alignment.mutant_aligned

    ref_position = 0
    mut_position = 0
    substitution_index = 0
    insertion_index = 0
    deletion_index = 0
    column = 0

    while column < len(reference_aligned):
        ref_nt = reference_aligned[column]
        mut_nt = mutant_aligned[column]

        if ref_nt == mut_nt:
            if ref_nt != "-":
                ref_position += 1
                mut_position += 1
            elif mut_nt != "-":
                mut_position += 1
            column += 1
            continue

        if ref_nt != "-" and mut_nt != "-":
            ref_position += 1
            mut_position += 1
            substitution_index += 1
            mutation_id = _event_label("SUB", substitution_index)
            ref_codon, mut_codon, ref_aa, mut_aa, effect = _classify_substitution_effect(
                reference_sequence, mutant_sequence, ref_position
            )
            mutations.append(
                NucleotideMutation(
                    mutation_id=mutation_id,
                    position=ref_position,
                    ref_nt=ref_nt,
                    mut_nt=mut_nt,
                    mutation_type="substitution",
                    ref_start=ref_position,
                    ref_end=ref_position,
                    mut_start=mut_position,
                    mut_end=mut_position,
                    ref_codon=ref_codon,
                    mut_codon=mut_codon,
                    ref_aa=ref_aa,
                    mut_aa=mut_aa,
                    effect=effect,
                    note="Single-base substitution observed in pairwise alignment.",
                )
            )
            column += 1
            continue

        if ref_nt == "-" and mut_nt != "-":
            insertion_index += 1
            mutation_id = _event_label("INS", insertion_index)
            inserted: list[str] = []
            start_mut = mut_position + 1
            anchor = ref_position
            while column < len(reference_aligned) and reference_aligned[column] == "-" and mutant_aligned[column] != "-":
                inserted.append(mutant_aligned[column])
                mut_position += 1
                column += 1
            inserted_sequence = "".join(inserted)
            ref_codon, mut_codon, ref_aa, mut_aa, effect = _classify_indel_effect(
                reference_sequence,
                mutant_sequence,
                "insertion",
                anchor,
                None,
                start_mut,
                "-",
                inserted_sequence,
            )
            mutations.append(
                NucleotideMutation(
                    mutation_id=mutation_id,
                    position=anchor,
                    ref_nt="-",
                    mut_nt=inserted_sequence,
                    mutation_type="insertion",
                    ref_start=anchor if anchor > 0 else None,
                    ref_end=anchor if anchor > 0 else None,
                    mut_start=start_mut,
                    mut_end=mut_position,
                    ref_codon=ref_codon,
                    mut_codon=mut_codon,
                    ref_aa=ref_aa,
                    mut_aa=mut_aa,
                    effect=effect,
                    note="Inserted sequence is anchored to the preceding reference base.",
                )
            )
            continue

        deletion_index += 1
        mutation_id = _event_label("DEL", deletion_index)
        deleted: list[str] = []
        start_ref = ref_position + 1
        mut_anchor = mut_position if mut_position > 0 else 1
        while column < len(reference_aligned) and reference_aligned[column] != "-" and mutant_aligned[column] == "-":
            deleted.append(reference_aligned[column])
            ref_position += 1
            column += 1
        deleted_sequence = "".join(deleted)
        ref_codon, mut_codon, ref_aa, mut_aa, effect = _classify_indel_effect(
            reference_sequence,
            mutant_sequence,
            "deletion",
            start_ref,
            start_ref,
            mut_anchor,
            deleted_sequence,
            "-",
        )
        mutations.append(
            NucleotideMutation(
                mutation_id=mutation_id,
                position=start_ref,
                ref_nt=deleted_sequence,
                mut_nt="-",
                mutation_type="deletion",
                ref_start=start_ref,
                ref_end=ref_position,
                mut_start=mut_anchor if mut_anchor > 0 else None,
                mut_end=mut_anchor if mut_anchor > 0 else None,
                ref_codon=ref_codon,
                mut_codon=mut_codon,
                ref_aa=ref_aa,
                mut_aa=mut_aa,
                effect=effect,
                note="Deleted sequence spans one or more reference bases.",
            )
        )

    return mutations


def _build_summary(result: AnalysisResult) -> AnalysisSummary:
    mutation_types = Counter(mutation.mutation_type for mutation in result.nucleotide_mutations)
    effects = Counter(mutation.effect for mutation in result.nucleotide_mutations)

    return AnalysisSummary(
        total_mutations=len(result.nucleotide_mutations),
        substitutions=mutation_types.get("substitution", 0),
        insertions=mutation_types.get("insertion", 0),
        deletions=mutation_types.get("deletion", 0),
        silent=effects.get("silent", 0),
        missense=effects.get("missense", 0),
        nonsense=effects.get("nonsense", 0),
        frameshift=effects.get("frameshift", 0),
        inframe_insertion=effects.get("inframe_insertion", 0),
        inframe_deletion=effects.get("inframe_deletion", 0),
        identity=result.alignment.identity,
        alignment_columns=result.alignment.columns,
        reference_length=result.reference.length,
        mutant_length=result.mutant.length,
    )


def _first_frameshift_codon(mutations: list[NucleotideMutation]) -> int | None:
    candidates = [
        codon_index_from_position(mutation.position if mutation.position > 0 else 1)
        for mutation in mutations
        if mutation.effect == "frameshift"
    ]
    candidates = [candidate for candidate in candidates if candidate is not None]
    return min(candidates) if candidates else None


def analyze_records(reference: SequenceRecordData, mutant: SequenceRecordData) -> AnalysisResult:
    alignment = align_sequences(reference.sequence, mutant.sequence)
    nucleotide_mutations = _detect_nucleotide_mutations(alignment, reference.sequence, mutant.sequence)
    first_frameshift_codon = _first_frameshift_codon(nucleotide_mutations)

    result = AnalysisResult(
        reference=reference,
        mutant=mutant,
        alignment=alignment,
        nucleotide_mutations=nucleotide_mutations,
        codon_changes=build_codon_comparisons(reference.sequence, mutant.sequence, first_frameshift_codon),
        amino_acid_changes=build_amino_acid_comparisons(reference.sequence, mutant.sequence, first_frameshift_codon),
    )
    result.summary = _build_summary(result)
    return result


def analyze_sequences(
    reference_sequence: str,
    mutant_sequence: str,
    reference_id: str = "reference",
    mutant_id: str = "mutant",
) -> AnalysisResult:
    reference = SequenceRecordData(
        identifier=reference_id,
        description=reference_id,
        sequence=reference_sequence,
        length=len(reference_sequence),
        source_name=reference_id,
    )
    mutant = SequenceRecordData(
        identifier=mutant_id,
        description=mutant_id,
        sequence=mutant_sequence,
        length=len(mutant_sequence),
        source_name=mutant_id,
    )
    return analyze_records(reference, mutant)


def analyze_files(reference_path: str | Path, mutant_path: str | Path) -> AnalysisResult:
    reference = parse_fasta_file(reference_path)
    mutant = parse_fasta_file(mutant_path)
    return analyze_records(reference, mutant)
