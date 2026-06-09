# MitoConda - 12S Amplicon Nanopore Pipeline

![MitoConda Logo](resources/logo/MitoConda_Logo.jpeg)

[![Language](https://img.shields.io/badge/language-Python%20%26%20Shell-blue)](https://www.python.org/)
[![Workflow](https://img.shields.io/badge/workflow-Snakemake-orange)](https://snakemake.readthedocs.io/)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)

---

## Overview

MitoConda is a **modular, reproducible** Snakemake pipeline designed for the analysis of 12S amplicon Nanopore sequencing data. It integrates state‑of‑the‑art tools for quality control, taxonomic diversity, and phylogenetic reconstruction, all wrapped in a single, easy‑to‑use bash script.

---

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [1️⃣ Prepare input data](#prepare-input-data)
  - [2️⃣ Run the pipeline](#run-the-pipeline)
  - [3️⃣ Inspect results](#inspect-results)
- [Output Structure](#output-structure)
- [Contributing](#contributing)

---

## Features

The pipeline automates the following analysis steps:

- **Quality Control** – Filtering and trimming of raw Nanopore reads using `chopper`.
- **Diversity Analysis** – Taxonomic classification with `blast` + `vsearch`.
- **LCA Assignment** – Lowest Common Ancestor taxonomy inference.
- **Rarefaction Curves** – Visual assessment of sequencing depth.
- **Phylogeny** – MAFFT alignment, FastTree tree building, and post‑processing.

---

## Installation

```bash
# 1️⃣ Clone the repository
git clone https://github.com/matheus-cosentino/MitoConda.git
cd MitoConda

# 2️⃣ Create the base Conda environment (contains Snakemake)
conda env create -f envs/mitoconda.yaml
conda activate mitoconda
```
> **Note**: All other tool‑specific environments are generated automatically on the first pipeline run.

---

## Configuration

All parameters are stored in `config/config.yaml`. The most common options are:

| Parameter   | Description |
|------------|-------------|
| `output_dir` | Directory where results are written (default: `"results"`). |
| `data_dir`   | Folder with raw FASTQ files (default: `"data"`). |
| `phylo_target` | TaxID(s) to target for phylogenetic analysis (e.g., `40674`). |
| `phylo_genes`  | Genes for the phylogeny step (e.g., `12s`). |
| `chopper.*`   | Quality‑filter thresholds (see below). |

### Example – tuning Chopper
```yaml
chopper:
  min_quality: [10]   # lower quality threshold (default was 15)
  min_length:  [400]   # longer reads are kept (default 200)
  max_length:  [800]
```

---

## Usage

### 1️⃣ Prepare input data
Place your raw Nanopore FASTQ files (`sample1.fastq.gz`, `sample2.fastq.gz`, …) in the directory indicated by `data_dir` (or specify an alternative with `--input`).

### 2️⃣ Run the pipeline (recommended)
Use the bundled wrapper script to activate environments and toggle modules on‑the‑fly:

```bash
# Full analysis (QC + diversity + phylogeny)
bash MitoConda.sh \
    --input data \
    --output results \
    --diversity \
    --qc \
    --phylo
```

#### Common flags
- `--diversity` – run BLAST/VSEARCH diversity workflow.
- `--qc` – enable the Chopper QC step and generate a MultiQC report.
- `--phylo` – perform phylogenetic alignment & tree building.
- `--build-db-only` – only create the reference database for a given TaxID/gene.
- `--taxid <ID>` – override `phylo_target` from the config.
- `--gene <NAME>` – override `phylo_genes`.

### 3️⃣ Inspect results
- **MultiQC report** – `results/multiqc_all/<pident>_multiqc_report.html`
- **Abundance tables** – `results/*/Abundance/`
- **Phylogenetic trees** – `results/*/Phylo/*_Aligned.fasta` and corresponding Newick files.

You can also visualise the workflow graph:

```bash
snakemake --dag | dot -Tsvg > resources/logo/dag.svg
```

---

## Output Structure

```
results/
└── <sample_name>/
    ├── QC/                # Chopper outputs & NanoStat stats
    ├── Abundance/         # Taxonomic abundance tables
    ├── Blast/             # BLAST results per sample
    ├── Vsearch/           # Clustering outputs
    ├── Phylo/             # Renamed FASTA, aligned FASTA, cleaned alignments, trees
    └── multiqc_all/       # Combined MultiQC report
```

---

## Contributing

We welcome contributions! Please fork the repository, create a feature branch, and submit a pull request. For major changes, open an issue first to discuss your ideas.

---

*Happy analyzing!*
