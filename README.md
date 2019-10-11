# rnaseq

A simple RNA-Seq differential gene expression pipeline using nextflow

```bash
nextflow main.nf --help
```

If all dependencies are installed (e.g. in a conda env) you should be able to execute:

```bash
nextflow main.nf --threads 2 --reads input.se.eco.csv --reference data/eco/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa --annotation data/eco/Escherichia_coli_k_12.ASM80076v1.44.gtf
```