# SeqDelta

SeqDelta is a bioinformatics tool for DNA mutation analysis and protein impact interpretation. It compares a reference FASTA against a mutant FASTA, performs explicit pairwise alignment, detects nucleotide-level differences, translates coding consequences, and renders both terminal and browser-grade reports.

The analysis is intentionally separated into three layers:

1. nucleotide-level mutations
2. codon-level changes
3. protein-level consequences

That separation carries through the CLI, the HTML report, the JSON/CSV exports, and the FastAPI web interface.

## Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/milord-x/SeqDelta.git
```

After installation, the primary CLI is:

```bash
seqdelta compare reference.fasta mutant.fasta --html-report report.html
```

`python cli.py` remains available as an internal development fallback, but it is not the primary public interface.

## What SeqDelta does

- parses arbitrary FASTA inputs, including real NCBI downloads
- runs explicit global pairwise alignment
- detects substitutions, insertions, and deletions
- classifies silent, missense, nonsense, frameshift, and in-frame indel effects
- distinguishes raw nucleotide events from codon and protein consequences
- exports results as terminal output, JSON, CSV, and standalone HTML
- serves a local web interface for interactive upload and presentation use

## CLI usage

Basic comparison:

```bash
seqdelta compare reference.fasta mutant.fasta
```

Generate HTML, JSON, and CSV outputs:

```bash
seqdelta compare reference.fasta mutant.fasta --protein-view --json-out results.json --csv-out mutations.csv --html-report report.html
```

Multiline form:

```bash
seqdelta compare \
  reference.fasta \
  mutant.fasta \
  --protein-view \
  --json-out results.json \
  --csv-out mutations.csv \
  --html-report report.html
```

CLI help:

```bash
seqdelta compare --help
```

## Web app usage

Launch the local FastAPI app:

```bash
python app.py
```

Then open:

```text
http://127.0.0.1:8000
```

The web interface supports:

- reference FASTA upload
- mutant FASTA upload
- full scientific results dashboard
- standalone HTML export
- JSON export
- CSV export

## Example sample data

Sample FASTA files are included in [`sample_data/`](/home/proxy/Projects/SeqDelta/sample_data):

- `reference.fasta`
- `mutant_silent.fasta`
- `mutant_missense.fasta`
- `mutant_nonsense.fasta`
- `mutant_insertion_frameshift.fasta`
- `mutant_deletion_frameshift.fasta`

Quick demo:

```bash
seqdelta compare sample_data/reference.fasta sample_data/mutant_missense.fasta --html-report demo.html
```

## Screenshots

Placeholder for:

- GitHub Pages landing page
- browser results dashboard
- mutation-focused alignment view
- full alignment view
- CLI terminal output

## Project structure

```text
SeqDelta/
├── app.py
├── cli.py
├── docs/
├── pyproject.toml
├── README.md
├── sample_data/
├── seqdelta/
│   ├── alignment.py
│   ├── cli.py
│   ├── models.py
│   ├── mutation.py
│   ├── parser.py
│   ├── report.py
│   ├── translation.py
│   ├── visualization.py
│   ├── web.py
│   └── assets/
├── static/
├── templates/
└── tests/
```

## GitHub Pages

The project website lives in [`docs/index.html`](/home/proxy/Projects/SeqDelta/docs/index.html). Enable GitHub Pages in the repository settings and set the source to the `/docs` folder.

## Limitations

- The current MVP assumes the coding frame begins at nucleotide 1.
- The standard genetic code is used for translation.
- One FASTA record per file is expected.
- Complex transcript-aware annotation and HGVS normalization are not yet implemented.
- Large structural variants and splice-aware consequences are out of scope.
- Frameshift interpretation is sequence-based, not transcript database-driven.

## Roadmap

- configurable coding frame and coding window
- transcript/exon-aware annotation
- HGVS-style variant naming
- protein domain overlays
- batch processing for multiple mutant samples
- packaged demo screenshots for the GitHub Pages site

## Development

Local fallback:

```bash
python cli.py compare sample_data/reference.fasta sample_data/mutant_nonsense.fasta --protein-view
```

Run tests:

```bash
pytest
```
