from __future__ import annotations

from io import StringIO
from pathlib import Path

from Bio import SeqIO

from .models import SequenceRecordData

VALID_DNA = set("ACGTRYSWKMBDHVN")


def normalize_sequence(sequence: str) -> str:
    cleaned = "".join(character for character in sequence.upper() if not character.isspace())
    invalid = sorted({character for character in cleaned if character not in VALID_DNA})
    if invalid:
        invalid_text = ", ".join(invalid)
        raise ValueError(f"Unsupported nucleotide symbols detected: {invalid_text}")
    return cleaned


def _record_to_model(record, default_identifier: str, source_name: str | None = None) -> SequenceRecordData:
    identifier = record.id or default_identifier
    description = record.description or identifier
    sequence = normalize_sequence(str(record.seq))
    if not sequence:
        raise ValueError("FASTA record is empty.")
    return SequenceRecordData(
        identifier=identifier,
        description=description,
        sequence=sequence,
        length=len(sequence),
        source_name=source_name,
    )


def parse_fasta_text(contents: str, default_identifier: str = "sequence", source_name: str | None = None) -> SequenceRecordData:
    handle = StringIO(contents)
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise ValueError("No FASTA records detected.")
    if len(records) > 1:
        raise ValueError("SeqDelta expects exactly one FASTA record per file.")
    return _record_to_model(records[0], default_identifier, source_name=source_name)


def parse_fasta_file(path: str | Path) -> SequenceRecordData:
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    with fasta_path.open("r", encoding="utf-8") as handle:
        contents = handle.read()
    return parse_fasta_text(contents, default_identifier=fasta_path.stem, source_name=fasta_path.name)
