# Protein Function Prediction - Team 2

Meeting 2: prepare a reproducible protein-sequence / Gene Ontology dataset for
*Candida albicans*, *Candida parapsilosis*, and *Saccharomyces cerevisiae*.
No predictive model has been trained yet.

## Start here

The main report is **`meeting_02_data_preparation.ipynb`**, written in English with
saved outputs, explanations, code, diagrams, and seven figures. Open it in Jupyter,
VS Code, or GitHub to read the existing results; no dataset execution is needed.

The repository also contains the raw download snapshot and prepared dataset,
so teammates do not need to retrieve changing database releases themselves.

### Quick start: view outputs and use the prepared data

Use Ubuntu / WSL2. Python **3.14.4** is the environment tested for this snapshot.
You need Python venv support. Full pipeline execution also requires `/usr/bin/time`
(Ubuntu package `time`).

```bash
cd Protein-function-prediction
bash scripts/bootstrap.sh --skip-mmseqs
bash start_notebook.sh
```

Bootstrap creates `.venv`, installs the exact package versions from
`requirements-lock.txt`, registers the notebook kernel in that environment,
and verifies the dataset. Select **Python (BioInformatics)** in Jupyter.
Internet access is required to install packages.

### Run the complete preparation pipeline

```bash
bash scripts/bootstrap.sh
bash start_notebook.sh
```

Then run all notebook cells from top to bottom. Full clustering requires Linux
x86_64 with SSE4.1. Setup downloads the original MMseqs2 build and checks its
recorded archive hash and version. If the upstream `latest` archive changes,
setup stops and asks for the recorded archive instead of silently changing tools.

The notebook verifies and reuses the included raw snapshot. Because MMseqs2 run
records include local paths, moving the project can cause clustering/search to
run again. The 128-restart split is recomputed; allow several minutes for Run All.
Do not run multiple preparation notebooks against the same files concurrently.

### Verify the snapshot or execute the integration check

```bash
python3 scripts/verify_dataset.py
# After full bootstrap, using the isolated environment:
.venv/bin/python scripts/check_notebook.py
```

The first command uses only the Python standard library and checks hashes,
split/group integrity, full FASTA sequences, labels, and saved notebook outputs.
The second executes the entire notebook and checks that the dataset remains
identical to the committed snapshot; it updates saved outputs and run reports.

## What is included?

| Path | Purpose |
|---|---|
| `meeting_02_data_preparation.ipynb` | Canonical English meeting report |
| `src/data_loading.py` | Download, cache, and load GAF / UniProt records |
| `src/meeting2_preparation.py` | Evidence filtering, FASTA, grouping, splitting, validation |
| `scripts/` | Environment/tool setup and verification |
| `data/raw/` | Original GAF, reviewed TSV, GO evidence JSONL, and metadata |
| `data/processed/` | Frozen master tables, FASTA, labels, vocabulary, and splits |
| `data/interim/` | Small clustering/search results needed for the homology audit |
| `reports/` | Statistics, figures, provenance, and integrity results |
| `docs/source_manifest.json` | Original source versions, URLs, and checksums |
| `requirements.txt` | Direct dependencies with version ranges |
| `requirements-lock.txt` | Exact Python package versions from the tested environment |

`.venv`, external MMseqs2 binaries, large regenerated evidence audits, temporary
databases, and logs are deliberately ignored. Raw and prepared snapshot files
are tracked so the whole team has the same benchmark after cloning.

## Dataset snapshot

| Organism | Reviewed | Retained | Train | Validation | Test |
|---|---:|---:|---:|---:|---:|
| C. albicans | 1466 | 704 | 564 | 70 | 70 |
| C. parapsilosis | 23 | 9 | 7 | 1 | 1 |
| S. cerevisiae | 7919 | 5543 | 4434 | 554 | 555 |
| Total | 9408 | 6256 | 5005 | 625 | 626 |

There are 39,211 protein-GO pairs, 5,265 distinct direct GO terms, and 5,481 final
similarity groups. `eval` means validation. Dataset fingerprints are recorded
in `reports/dataset_manifest.json`. Cluster rows are sorted by protein ID for
deterministic exports; this normalization does not change any membership or
split assignment (see `reports/snapshot_normalization.json`).

Final labels must exist in UniProt and have allowed evidence in UniProt or a
matching GAF; exact GAF NOT conflicts are excluded. Sequence similarity grouping
uses MMseqs2 identity 0.5, coverage 0.8, cov-mode 0, plus duplicate and search-hit
merging. Whole groups stay in one split. Search is heuristic, not an exhaustive
proof of the absence of all similar cross-split pairs.

Limitations: strong species imbalance, only 9 retained C. parapsilosis proteins,
and GO terms in validation/test absent from training. Missing annotation is not
proof of absence of function. GO hierarchy propagation and final evaluation
policy belong to later project stages.

## Team workflow

Use this snapshot for the first models. Do not silently replace its raw files or
split assignments when adding organisms or new database releases. Create a new
dataset version and recheck similarity leakage. Keep predicted scores separate
from reference annotations.

The supplied repository LICENSE is retained. Public source data and third-party
tools retain their upstream provenance; see `docs/DATA_SOURCES.md`.
