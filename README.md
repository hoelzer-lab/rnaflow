
![](https://img.shields.io/badge/nextflow-20.07.1-brightgreen)
<!--![](https://img.shields.io/badge/uses-docker-blue.svg)-->

# RN(ext)A-Seq - An effective and simple RNA-Seq differential gene expression pipeline using Nextflow

## Quick installation

The pipeline is written in `Nextflow`, which can be used on any POSIX compatible system (Linux, OS X, etc). Windows system is supported through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

1. Install  [`nextflow`](https://nf-co.re/usage/installation)
    <details><summary>click here for a bash one-liner </summary>

    ```bash
    wget -qO- https://get.nextflow.io | bash
    ```

    </details>
1. Install [`conda`](https://conda.io/miniconda.html)
    <details><summary>click here for a bash two-liner for Miniconda3 Linux 64-bit</summary>

    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

    </details>
1. For transcriptome assembly: install [`docker`](https://docs.docker.com/engine/installation/)

All other dependencies and tools will be installed within the pipeline via `conda` or `docker`.

## Quick start

- Start a test run

```bash
nextflow run hoelzer/rnaseq -profile test,conda,local
```

- Call help

```bash
nextflow run hoelzer/rnaseq --help
```

## Usage

```bash
nextflow run hoelzer/rnaseq --reads test-data/input.se.hsa_small.csv --species hsa --max_cores 6 --cores 2
```

with `hsa`, `mmu`, `mau` or `eco` [build-in species](#build-in-species), or define your own genome reference and annotation files in CSV files:

```bash
nextflow run hoelzer/rnaseq --reads test-data/input.se.hsa_small.csv --genome fastas.csv --annotation gtf.csv --max_cores 6 --cores 2
```

Genomes and annotations from `--genome` and `--annotation` (and `--species`) are concatenated.
By default, all possible comparisons are performed. Use `--deg` to change this.

### Build-in species

| Species      | three-letter shortcut | Genome                              | Annotation                                    |
| ------------ | --------------------- | ----------------------------------- | --------------------------------------------- |
| Homo sapiens | `hsa`                 | Homo_sapiens.GRCh38.98              | Homo_sapiens.GRCh38.dna.primary_assembly      |
| Mus musculus | `mmu`                 | Mus_musculus.GRCm38.99              | Mus_musculus.GRCm38.dna.primary_assembly      |
| Homo sapiens | `mau`                 | Mesocricetus_auratus.MesAur1.0.100  | Mesocricetus_auratus.MesAur1.0.dna.toplevel   |
| Homo sapiens | `eco`                 | Escherichia_coli_k_12.ASM80076v1.45 | Escherichia_coli_k_12.ASM80076v1.dna.toplevel |

## Help message

```bash
Usage example:
nextflow run main.nf --cores 4 --reads input.csv --species eco
or
nextflow run main.nf --cores 4 --reads input.csv --species eco --assembly
or
nextflow run main.nf --cores 4 --reads input.csv --genome fastas.csv --annotation gtfs.csv
or
nextflow run main.nf --cores 4 --reads input.csv --genome fastas.csv --annotation gtfs.csv --species eco
Genomes and annotations from --genome, --annotation and --species are concatenated

Input:
--reads                  a CSV file following the pattern: Sample,R,Condition,Patient for single-end or Sample,R1,R2,Condition,Patient for paired-end
                                    (check terminal output if correctly assigned)
                                    In default all possible comparisons of conditions in one direction are made. Use --deg to change this.
--species                reference genome and annotation with automatic download.
                                    Currently supported are:
                                    - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                    - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                    - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]
                                    - mau [Ensembl: Mesocricetus_auratus.MesAur1.0.dna.toplevel | Mesocricetus_auratus.MesAur1.0.100]
--genome                 CSV file with genome reference FASTA files (one path in each line).
                                    If set, --annotation must also be set.
--annotation             CSV file with genome annotation GTF files (one path in each line)

Options:
--assembly               perform de novo and reference-based transcriptome assembly instead of DEG analysis [default $params.assembly]
--deg                    a CSV file following the pattern: conditionX,conditionY
                            Each line stands for one differential gene expression comparison.
--index                  the path to the hisat2 index prefix matching the genome provided via --species. 
                            If provided, no new index will be build. Must be named 'index.*.ht2'.  
                            Simply provide the path like 'data/db/index'. DEPRECATED
--mode                   either 'single' (single-end) or 'paired' (paired-end) sequencing [default $params.mode]
--strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default $params.strand]
--tpm                    threshold for TPM (transcripts per million) filter. A feature is discared, 
                            if in all conditions the mean TPM value of all libraries in this condition are below the threshold. [default $params.tpm]
--skip_sortmerna         skip rRNA removal via SortMeRNA [default $params.skip_sortmerna] 
--busco_db               the database used with BUSCO [default: $params.busco_db]
                            full list of available data sets at https://busco.ezlab.org/v2/frame_wget.html 
--dammit_uniref90        add UniRef90 to the dammit databases  [default: $params.dammit_uniref90]

Computing options:
--cores                  max cores per process for local use [default $params.cores]
--max_cores              max cores used on the machine for local use [default $params.max_cores]
--memory                 max memory in GB for local use [default $params.memory]
--output                 name of the result folder [default $params.output]

--permanentCacheDir      location for auto-download data like databases [default $params.permanentCacheDir]
--condaCacheDir          location for storing the conda environments [default $params.condaCacheDir]
--workdir                working directory for all intermediate results [default $params.workdir]

Nextflow options:
-with-report rep.html    cpu / ram usage (may cause errors)
-with-dag chart.html     generates a flowchart for the process tree
-with-timeline time.html timeline (may cause errors)

Profile:
-profile                 standard (local and conda),local, conda, slurm, ara (slurm, conda and customization) [default standard]
```

<!-- # Flow chart

![flow-chart](figures/chart.png) -->
