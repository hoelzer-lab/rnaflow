
![](https://img.shields.io/badge/nextflow-19.10.0-brightgreen)
<!--![](https://img.shields.io/badge/uses-docker-blue.svg)-->

# RN(ext)A-Seq

A simple RNA-Seq differential gene expression pipeline using nextflow

```bash
nextflow run main.nf --help
```

Dependencies will automatically be installed via conda, just execute:

```bash
nextflow run main.nf --max_cores 6 --cores 2 --reads input.se.eco.csv --species eco
```

with `eco`, `mmu` or `hsa`, or define your own genome reference and annotation files in CSV files:

```bash
nextflow run main.nf --max_cores 6 --cores 2 --reads input.se.eco.csv --genome fastas.csv --annotation gtf.csv
```

Genomes and annotations from `--genome` and `--annotation` (and `--species`) are concatenated.
By default, all possible comparisons are performed. Use `--dge` to change this.

# Flow chart

![flow-chart](figures/chart.png)
