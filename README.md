
![](https://img.shields.io/badge/nextflow-20.07.1-brightgreen)
![](https://img.shields.io/badge/uses-conda-yellow.svg)
<!--![](https://img.shields.io/badge/uses-docker-blue.svg)-->

# RN(ext)A-Seq - An effective and simple RNA-Seq differential gene expression pipeline using Nextflow

![flow-chart](figures/workflow.jpg)
*Figure 1* Workflow. The user can decide after preprocessing to run a differential gene expression (DEG) analysis or a transcriptome assembly. Circles symbolize input data and download icons symbolize automated download of resources. Steps marked by asterisks are currently only available for some species.

## Quick installation

The pipeline is written in [`Nextflow`](https://nf-co.re/usage/installation), which can be used on any POSIX compatible system (Linux, OS X, etc). Windows system is supported through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

1. Install  `Nextflow`
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

OR

1. Install `conda`
    <details><summary>click here for a bash two-liner for Miniconda3 Linux 64-bit</summary>

    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

    </details>
1. Install `Nextflow` via `conda`
    <details><summary>click here to see how to do that</summary>

    ```bash
    conda create -n nextflow -c bioconda nextflow
    conda active nextflow
    ```

    </details>

For transcriptome assembly, please install also [`docker`](https://docs.docker.com/engine/installation/).

All other dependencies and tools will be installed within the pipeline via `conda` or `docker`.

## Quick start

### Start a test run

```bash
# conda active nextflow
nextflow run hoelzer/rnaseq -profile test,conda,local
```

... performs 
- a differential gene expression analysis with sub-sampled human read data,
- comparing two conditions, 
- with a local execution (uses max. 4 cores and 8GB) and 
- `conda` dependency management. 

### Call help

```bash
nextflow run hoelzer/rnaseq --help
```

## Usage

```bash
nextflow run hoelzer/rnaseq --reads input.csv --species hsa --max_cores 6 --cores 2
```

with `hsa`, `mmu`, `mau` or `eco` [build-in species](#build-in-species), or define your own genome reference and annotation files in CSV files:

```bash
nextflow run hoelzer/rnaseq --reads input.csv --genome fastas.csv --annotation gtfs.csv --max_cores 6 --cores 2
```

Genomes and annotations from `--genome` and `--annotation` (and `--species`) are concatenated.
By default, all possible comparisons are performed. Use `--deg` to change this.

### Input files

#### Read files (required)

Specify your read files in `FASTQ` format with `--reads input.csv`. The files `input.csv` has to look like this for single-end reads:

```
Sample,R,Condition,Patient
mock_rep1,/path/to/reads/mock1.fastq.gz,mock,,
mock_rep2,/path/to/reads/mock2.fastq.gz,mock,,
mock_rep3,/path/to/reads/mock3.fastq.gz,mock,,
treated_rep1,/path/to/reads/treat1.fastq.gz,treated,,
treated_rep2,/path/to/reads/treat2.fastq.gz,treated,,
treated_rep3,/path/to/reads/treat3.fastq.gz,treated,,
```

and for paired-end reads, like this:

```
Sample,R1,R2,Condition,Patient
mock_rep1,/path/to/reads/mock1_1.fastq,/path/to/reads/mock1_2.fastq,mock,A
mock_rep2,/path/to/reads/mock2_1.fastq,/path/to/reads/mock2_2.fastq,mock,B
mock_rep3,/path/to/reads/mock3_1.fastq,/path/to/reads/mock3_2.fastq,mock,C
treated_rep1,/path/to/reads/treat1_1.fastq,/path/to/reads/treat1_2.fastq,treated,A
treated_rep2,/path/to/reads/treat2_1.fastq,/path/to/reads/treat2_2.fastq,treated,B
treated_rep3,/path/to/reads/treat3_1.fastq,/path/to/reads/treat3_2.fastq,treated,C
```

Read files can be compressed (`.bz`). You need at least two replicates for each condition to run the pipeline. Patient labels are optional.

#### Genomes and annotation

If you don't use one of the [build-in species](#build-in-species), specify your genomes via `--genome fastas.csv`, with `fastas.csv` looking like this:

```
/path/to/reference_genome1.fasta
/path/to/reference_genome2.fasta
```

and `--annotation gtfs.csv` with `gtfs.csv` looking like this:

```
/path/to/reference_annotation_1.gtf
/path/to/reference_annotation_2.gtf
```

#### Build-in species

| Species      | three-letter shortcut | Genome                              | Annotation                                    |
| ------------ | --------------------- | ----------------------------------- | --------------------------------------------- |
| Homo sapiens | `hsa`                 | Homo_sapiens.GRCh38.98              | Homo_sapiens.GRCh38.dna.primary_assembly      |
| Mus musculus | `mmu`                 | Mus_musculus.GRCm38.99              | Mus_musculus.GRCm38.dna.primary_assembly      |
| Homo sapiens | `mau`                 | Mesocricetus_auratus.MesAur1.0.100  | Mesocricetus_auratus.MesAur1.0.dna.toplevel   |
| Homo sapiens | `eco`                 | Escherichia_coli_k_12.ASM80076v1.45 | Escherichia_coli_k_12.ASM80076v1.dna.toplevel |

#### Comparisons for DEG analysis

Per default all possible comparisons in one direction are performed. To change this, please define the needed comparison with `--deg comparisons.csv`, where each line contains a pairwise comparison:

```
conditionX,conditionY
conditionA,conditionB
conditionB,conditionA
```

## Workflow control

### Preprocessing

```bash
--mode                          # either 'single' (single-end) or 'paired' (paired-end) sequencing [default single]
--skip_sortmerna                # skip rRNA removal via SortMeRNA [default false]
--fastp_additional_params       # additional parameters for fastp [default '-5 -3 -W 4 -M 20 -l 15 -x -n 5 -z 6']
--histat2_additional_params     # additional parameters for HISAT2
```

###  DEG analysis

```bash
--strand                        # strandness for counting with featureCounts: 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default 0]
--tpm                           # threshold for TPM (transcripts per million) filter [default 1]
--deg                           # a CSV file following the pattern: conditionX,conditionY
```

### Transcriptome assembly

```bash
--assemly                       # switch to transcriptome assemly
--busco_db                      # BUSCO database ['euarchontoglires_odb9']
--dammit_uniref90               # add UniRef90to dammit databases [false]
```

## Profiles/configuration options

Per default the pipeline is a local execution with `conda` dependency management (corresponds to `-profile local,conda`). Adjust this setting by combining an executer option with a engine option, e.g. `-profile local,conda` or `-profile slurm,conda`.

### Executor options...
*... or how to schedule your workload.*

Currently implemented are `local` and `slurm` executions.

You can customize `local` with this parameters:

```bash
--cores                         # cores for one process [default 1]
--max_cores                     # max. cores used in total [default allAvailable]
--memory                        # max. memory in GB for local use [default 8 GB]
```

### Engine options... 
*... or in which environment to run the tools.*

Currently implemented is `conda`. For transcriptome assembly some tools need to be run with `docker`.

`docker` support for all steps is coming soon!

<!-- ## Help message
```
``` -->

<!-- # Flow chart

![flow-chart](figures/chart.png) -->
