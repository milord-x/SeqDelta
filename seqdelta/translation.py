from __future__ import annotations

from Bio.Seq import Seq

from .models import AminoAcidComparison, CodonComparison


def translate_codon(codon: str | None) -> str | None:
    if codon is None or len(codon) != 3 or "-" in codon:
        return None
    return str(Seq(codon).translate(table=1))


def codon_index_from_position(position: int | None) -> int | None:
    if position is None or position < 1:
        return None
    return ((position - 1) // 3) + 1


def get_codon(sequence: str, codon_index: int | None) -> str | None:
    if codon_index is None or codon_index < 1:
        return None
    start = (codon_index - 1) * 3
    stop = start + 3
    if start >= len(sequence):
        return None
    codon = sequence[start:stop]
    return codon if len(codon) == 3 else None


def get_peptide(sequence: str) -> str:
    coding_length = len(sequence) - (len(sequence) % 3)
    if coding_length == 0:
        return ""
    return str(Seq(sequence[:coding_length]).translate(table=1))


def build_codon_comparisons(
    reference_sequence: str,
    mutant_sequence: str,
    first_frameshift_codon: int | None = None,
) -> list[CodonComparison]:
    max_codons = max((len(reference_sequence) + 2) // 3, (len(mutant_sequence) + 2) // 3)
    codons: list[CodonComparison] = []

    for codon_index in range(1, max_codons + 1):
        ref_codon = get_codon(reference_sequence, codon_index)
        mut_codon = get_codon(mutant_sequence, codon_index)
        ref_aa = translate_codon(ref_codon)
        mut_aa = translate_codon(mut_codon)
        changed = ref_codon != mut_codon

        if first_frameshift_codon and codon_index >= first_frameshift_codon and changed:
            effect = "frameshift"
        elif not changed:
            effect = "unchanged"
        elif ref_aa == mut_aa:
            effect = "silent"
        elif mut_aa == "*" and ref_aa != "*":
            effect = "nonsense"
        else:
            effect = "missense"

        codons.append(
            CodonComparison(
                codon_index=codon_index,
                ref_start=((codon_index - 1) * 3) + 1 if ref_codon else None,
                mut_start=((codon_index - 1) * 3) + 1 if mut_codon else None,
                ref_codon=ref_codon,
                mut_codon=mut_codon,
                ref_aa=ref_aa,
                mut_aa=mut_aa,
                changed=changed,
                effect=effect,
            )
        )

    return codons


def build_amino_acid_comparisons(
    reference_sequence: str,
    mutant_sequence: str,
    first_frameshift_codon: int | None = None,
) -> list[AminoAcidComparison]:
    ref_peptide = get_peptide(reference_sequence)
    mut_peptide = get_peptide(mutant_sequence)
    max_length = max(len(ref_peptide), len(mut_peptide))
    amino_acids: list[AminoAcidComparison] = []

    for index in range(max_length):
        aa_index = index + 1
        ref_aa = ref_peptide[index] if index < len(ref_peptide) else None
        mut_aa = mut_peptide[index] if index < len(mut_peptide) else None
        changed = ref_aa != mut_aa

        if first_frameshift_codon and aa_index >= first_frameshift_codon and changed:
            effect = "frameshift"
        elif not changed:
            effect = "unchanged"
        elif ref_aa == mut_aa:
            effect = "silent"
        elif mut_aa == "*" and ref_aa != "*":
            effect = "nonsense"
        else:
            effect = "missense"

        amino_acids.append(
            AminoAcidComparison(
                aa_index=aa_index,
                ref_aa=ref_aa,
                mut_aa=mut_aa,
                changed=changed,
                effect=effect,
            )
        )

    return amino_acids
