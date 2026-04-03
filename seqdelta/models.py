from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(slots=True)
class SequenceRecordData:
    identifier: str
    description: str
    sequence: str
    length: int
    source_name: str | None = None


@dataclass(slots=True)
class AlignmentResult:
    score: float
    method: str
    reference_aligned: str
    mutant_aligned: str
    match_line: str
    identity: float
    columns: int


@dataclass(slots=True)
class NucleotideMutation:
    mutation_id: str
    position: int
    ref_nt: str
    mut_nt: str
    mutation_type: str
    ref_start: int | None = None
    ref_end: int | None = None
    mut_start: int | None = None
    mut_end: int | None = None
    ref_codon: str | None = None
    mut_codon: str | None = None
    ref_aa: str | None = None
    mut_aa: str | None = None
    effect: str | None = None
    note: str | None = None


@dataclass(slots=True)
class CodonComparison:
    codon_index: int
    ref_start: int | None
    mut_start: int | None
    ref_codon: str | None
    mut_codon: str | None
    ref_aa: str | None
    mut_aa: str | None
    changed: bool
    effect: str


@dataclass(slots=True)
class AminoAcidComparison:
    aa_index: int
    ref_aa: str | None
    mut_aa: str | None
    changed: bool
    effect: str


@dataclass(slots=True)
class AnalysisSummary:
    total_mutations: int
    substitutions: int
    insertions: int
    deletions: int
    silent: int
    missense: int
    nonsense: int
    frameshift: int
    inframe_insertion: int
    inframe_deletion: int
    identity: float
    alignment_columns: int
    reference_length: int
    mutant_length: int


@dataclass(slots=True)
class AnalysisResult:
    reference: SequenceRecordData
    mutant: SequenceRecordData
    alignment: AlignmentResult
    nucleotide_mutations: list[NucleotideMutation] = field(default_factory=list)
    codon_changes: list[CodonComparison] = field(default_factory=list)
    amino_acid_changes: list[AminoAcidComparison] = field(default_factory=list)
    summary: AnalysisSummary | None = None
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
