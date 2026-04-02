from __future__ import annotations

from Bio import Align

from .models import AlignmentResult


def _reconstruct_aligned_strings(alignment) -> tuple[str, str]:
    target = str(alignment.target)
    query = str(alignment.query)
    coordinates = alignment.coordinates

    reference_chunks: list[str] = []
    mutant_chunks: list[str] = []

    for index in range(coordinates.shape[1] - 1):
        target_start = int(coordinates[0, index])
        target_end = int(coordinates[0, index + 1])
        query_start = int(coordinates[1, index])
        query_end = int(coordinates[1, index + 1])

        target_length = target_end - target_start
        query_length = query_end - query_start

        if target_length == query_length:
            reference_chunks.append(target[target_start:target_end])
            mutant_chunks.append(query[query_start:query_end])
            continue

        diagonal = min(target_length, query_length)
        if diagonal:
            reference_chunks.append(target[target_start : target_start + diagonal])
            mutant_chunks.append(query[query_start : query_start + diagonal])

        if target_length > diagonal:
            gap_length = target_length - diagonal
            reference_chunks.append(target[target_start + diagonal : target_end])
            mutant_chunks.append("-" * gap_length)

        if query_length > diagonal:
            gap_length = query_length - diagonal
            reference_chunks.append("-" * gap_length)
            mutant_chunks.append(query[query_start + diagonal : query_end])

    return "".join(reference_chunks), "".join(mutant_chunks)


def align_sequences(reference_sequence: str, mutant_sequence: str) -> AlignmentResult:
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(reference_sequence, mutant_sequence)[0]
    reference_aligned, mutant_aligned = _reconstruct_aligned_strings(alignment)

    match_line = []
    matches = 0
    for ref_nt, mut_nt in zip(reference_aligned, mutant_aligned):
        if ref_nt == mut_nt and ref_nt != "-":
            match_line.append("|")
            matches += 1
        elif "-" in (ref_nt, mut_nt):
            match_line.append(" ")
        else:
            match_line.append(".")

    identity = matches / max(len(reference_aligned), 1)

    return AlignmentResult(
        score=float(alignment.score),
        method="global_pairwise_alignment",
        reference_aligned=reference_aligned,
        mutant_aligned=mutant_aligned,
        match_line="".join(match_line),
        identity=identity,
        columns=len(reference_aligned),
    )
