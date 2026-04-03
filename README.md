# SeqDelta

SeqDelta is an educational, exploratory bioinformatics tool for sequence-level mutation analysis. It compares a reference DNA FASTA sequence against a mutant DNA FASTA sequence, performs explicit pairwise alignment, separates nucleotide-level differences from codon and protein consequences, and exports the result as a polished terminal summary, standalone HTML report, JSON, and CSV.

SeqDelta is designed for:

- classroom demos
- student portfolios
- exploratory sequence comparison
- lightweight local analysis workflows

It is not positioned as a clinical, diagnostic, or transcript-annotation tool.

## What SeqDelta Does

- Parses one DNA FASTA record from a reference file and one DNA FASTA record from a mutant file
- Runs explicit global pairwise alignment
- Detects substitutions, insertions, and deletions at the nucleotide level
- Separates sequence events from codon-level and amino-acid-level consequences
- Classifies silent, missense, nonsense, frameshift, and in-frame indel outcomes
- Generates a modern local web report for presentation and screenshots
- Exports machine-readable JSON and CSV outputs for downstream inspection

## Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/milord-x/SeqDelta.git
```

After installation, the primary interface is the `seqdelta` command:

```bash
seqdelta compare reference.fasta mutant.fasta --html-report report.html
```

## Quick Demo

Run the included missense example:

```bash
seqdelta compare examples/reference.fasta examples/missense.fasta --html-report demo.html
```

This generates:

- a terminal summary
- a standalone HTML report at `demo.html`
- optional JSON and CSV if requested

## CLI Usage

Basic comparison:

```bash
seqdelta compare reference.fasta mutant.fasta
```

Rich terminal output plus all export formats:

```bash
seqdelta compare reference.fasta mutant.fasta --protein-view --json-out results.json --csv-out mutations.csv --html-report report.html
```

Useful help commands:

```bash
seqdelta --help
seqdelta compare --help
```

Current CLI options:

- `--json-out`
- `--csv-out`
- `--html-report`
- `--protein-view`
- `--pretty/--plain`

## Local Web App

Launch the local FastAPI app:

```bash
python app.py
```

Then open:

```text
http://127.0.0.1:8000
```

The local web interface supports:

- reference FASTA upload
- mutant FASTA upload
- browser-based results dashboard
- mutation summary cards
- mutation-centered alignment snippets
- full alignment view
- downloadable HTML, JSON, and CSV outputs

## Example Data

The repository includes ready-to-run FASTA files in [`examples/`](/home/proxy/Projects/SeqDelta/examples):

- [`reference.fasta`](/home/proxy/Projects/SeqDelta/examples/reference.fasta) - baseline reference sequence
- [`substitution.fasta`](/home/proxy/Projects/SeqDelta/examples/substitution.fasta) - single-nucleotide substitution example
- [`missense.fasta`](/home/proxy/Projects/SeqDelta/examples/missense.fasta) - codon change with amino-acid substitution
- [`silent.fasta`](/home/proxy/Projects/SeqDelta/examples/silent.fasta) - codon change without amino-acid change
- [`frameshift.fasta`](/home/proxy/Projects/SeqDelta/examples/frameshift.fasta) - insertion-based frameshift example

The original development fixtures remain in [`sample_data/`](/home/proxy/Projects/SeqDelta/sample_data).

## Screenshots

Placeholder assets are included in [`docs/screenshots/`](/home/proxy/Projects/SeqDelta/docs/screenshots) so the repository is ready for real captures later.

Upload page placeholder:

![Upload page placeholder](/home/proxy/Projects/SeqDelta/docs/screenshots/upload-page-placeholder.svg)

Results dashboard placeholder:

![Results dashboard placeholder](/home/proxy/Projects/SeqDelta/docs/screenshots/results-dashboard-placeholder.svg)

Mutation table placeholder:

![Mutation table placeholder](/home/proxy/Projects/SeqDelta/docs/screenshots/mutation-table-placeholder.svg)

Alignment view placeholder:

![Alignment view placeholder](/home/proxy/Projects/SeqDelta/docs/screenshots/alignment-view-placeholder.svg)

Recommended real captures to replace these placeholders:

- landing/upload page
- full results dashboard
- nucleotide mutation table
- mutation-focused alignment view

## Project Structure

```text
SeqDelta/
├── app.py
├── cli.py
├── docs/
│   ├── index.html
│   └── screenshots/
├── examples/
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

## Output Modes

SeqDelta currently supports:

- CLI-first terminal analysis
- local web dashboard for browser-based presentation
- standalone HTML report generation
- JSON export
- CSV export

## Limitations

SeqDelta is intentionally lightweight and sequence-centric. Current limitations are explicit:

- It assumes the coding frame starts at nucleotide position 1.
- It does not accept CDS annotations or reading-frame offsets.
- It has no transcript or exon awareness.
- Codon consequence mapping is sequence-based rather than annotation-based.
- It expects one FASTA record per file.
- It is designed for DNA sequence comparison, not full genome annotation workflows.
- It does not implement HGVS normalization or clinical variant interpretation.
- It is not suitable for diagnostic or clinical decision-making.

## Positioning

SeqDelta should be understood as:

- an educational bioinformatics project
- an exploratory sequence analysis tool
- a sequence-level mutation comparison workflow

It should not be described as a diagnostic system or medical interpretation platform.

## Roadmap

- configurable reading frame and coding window
- transcript-aware consequence mapping
- richer protein interpretation overlays
- stronger figure export and presentation assets
- more realistic example datasets

## Development

Internal fallback command:

```bash
python cli.py compare examples/reference.fasta examples/missense.fasta --html-report demo.html
```

Run tests:

```bash
pytest
```
