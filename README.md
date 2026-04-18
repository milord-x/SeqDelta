# SeqDelta

SeqDelta is an educational, exploratory bioinformatics tool for sequence-level mutation analysis. It compares a reference DNA FASTA sequence against a mutant DNA FASTA sequence, performs explicit pairwise alignment, separates nucleotide-level differences from codon and protein consequences, and exports the result as terminal output, standalone HTML, JSON, and CSV.

SeqDelta is intended for:

- classroom demos
- student portfolios
- exploratory sequence comparison
- lightweight local analysis workflows

It is not a clinical or diagnostic tool.

## Features

- One-record FASTA input for reference and mutant sequences
- Global pairwise alignment
- Nucleotide-level mutation detection
- Codon-level consequence mapping
- Protein-level consequence comparison
- Effect classification for silent, missense, nonsense, frameshift, and in-frame indels
- CLI-first workflow with optional local web interface
- HTML, JSON, and CSV exports

## Installation

Install the CLI directly from GitHub:

```bash
pip install git+https://github.com/milord-x/SeqDelta.git
```

Then run:

```bash
seqdelta compare reference.fasta mutant.fasta --html-report report.html
```

## Download the Project

If you want the full repository locally:

```bash
git clone https://github.com/milord-x/SeqDelta.git
cd SeqDelta
```

## Full Local Setup

If you want both the CLI and the local web app in a clean environment:

```bash
git clone https://github.com/milord-x/SeqDelta.git
cd SeqDelta
python -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
```

Verify the install:

```bash
seqdelta --help
seqdelta compare --help
pytest
```

## Quick Demo

```bash
seqdelta compare examples/reference.fasta examples/missense.fasta --html-report demo.html
```

## CLI Usage

Basic comparison:

```bash
seqdelta compare reference.fasta mutant.fasta
```

Export all major formats:

```bash
seqdelta compare reference.fasta mutant.fasta --protein-view --json-out results.json --csv-out mutations.csv --html-report report.html
```

Current CLI options:

- `--json-out`
- `--csv-out`
- `--html-report`
- `--protein-view`
- `--pretty/--plain`

## Local Web App

Run the local FastAPI app:

```bash
python app.py
```

Then open:

```text
http://127.0.0.1:8000
```

The local web interface supports FASTA upload, browser-based result viewing, and downloadable HTML, JSON, and CSV outputs.

## Example Data

Ready-to-run FASTA files are included in [examples/](examples/):

- [reference.fasta](examples/reference.fasta)
- [substitution.fasta](examples/substitution.fasta)
- [missense.fasta](examples/missense.fasta)
- [silent.fasta](examples/silent.fasta)
- [frameshift.fasta](examples/frameshift.fasta)

Development fixtures are kept in [sample_data/](sample_data/).

## Screenshots

![Demo](screenshots/recording.gif)

The current recording shows the browser workflow end to end. Additional still screenshots can be added later for the upload page, results dashboard, mutation table, and alignment view.

## Project Structure

```text
SeqDelta/
├── app.py
├── cli.py
├── examples/
├── screenshots/
├── pyproject.toml
├── README.md
├── sample_data/
├── seqdelta/
├── static/
├── templates/
└── tests/
```

## Limitations

- Assumes the coding frame starts at nucleotide position 1
- No CDS annotation or reading-frame offset input
- No transcript or exon awareness
- Codon consequence mapping is sequence-based
- Expects one FASTA record per file
- No HGVS normalization or clinical interpretation
- Not suitable for diagnostic or medical use

## Positioning

SeqDelta should be presented as:

- an educational bioinformatics project
- an exploratory sequence analysis tool
- a sequence-level mutation comparison workflow

## Development

Fallback local command:

```bash
python cli.py compare examples/reference.fasta examples/missense.fasta --html-report demo.html
```

Run tests:

```bash
pytest
```
