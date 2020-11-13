
![](https://img.shields.io/github/v/release/hoelzer-lab/rnaseq)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)
![](https://github.com/hoelzer-lab/rnaseq/workflows/Syntax_check/badge.svg)

![](https://img.shields.io/badge/nextflow-20.07.1-brightgreen)
![](https://img.shields.io/badge/uses-conda-yellow.svg)
![](https://img.shields.io/badge/uses-docker-blue.svg)

<!--[![Generic badge](https://img.shields.io/badge/Publication-bioRxiv-red.svg)](https://www.)-->

# RN(ext)A-Seq - An effective and simple RNA-Seq differential gene expression pipeline using Nextflow

![flow-chart](figures/workflow.jpg)
*Figure 1.* Workflow. The user can decide after preprocessing to run a differential gene expression (DEG) analysis or a transcriptome assembly. Circles symbolize input data and download icons symbolize automated download of resources. Steps marked by asterisks are currently only available for some species.

## Quick installation

The pipeline is written in [`Nextflow`](https://nf-co.re/usage/installation), which can be used on any POSIX compatible system (Linux, OS X, etc). Windows system is supported through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux).

1. Install  `Nextflow`
    <details><summary>click here for a bash one-liner </summary>

    ```bash
    wget -qO- https://get.nextflow.io | bash
    # In the case you don’t have wget
    # curl -s https://get.nextflow.io | bash
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
nextflow run hoelzer-lab/rnaseq -profile test,conda,local
```

... performs

- a differential gene expression analysis with sub-sampled human read data,
- on a reduced human genome and annotation (chromosome 1, 10 and 11),
- comparing two conditions (MAQCA, MAQCB),
- with a local execution (uses max. 4 cores in total and 8GB) and
- `conda` dependency management.

### Call help

```bash
nextflow run hoelzer-lab/rnaseq --help
```

## Usage

```bash
nextflow run hoelzer-lab/rnaseq --reads input.csv --species hsa --include_species --max_cores 6 --cores 2
```

with `hsa`, `mmu`, `mau` or `eco` [build-in species](#build-in-species), or define your own genome reference and annotation files in CSV files:

```bash
nextflow run hoelzer-lab/rnaseq --reads input.csv --genome fastas.csv --annotation gtfs.csv --max_cores 6 --cores 2
```

Genomes and annotations from `--species`, if `--include_species` is set, `--genome` and `--annotation` are concatenated.
By default, all possible comparisons are performed. Use `--deg` to change this.

### Input files

#### Read files (required)

Specify your read files in `FASTQ` format with `--reads input.csv`. The file `input.csv` has to look like this for single-end reads:

```csv
Sample,R,Condition,Patient
mock_rep1,/path/to/reads/mock1.fastq.gz,mock,
mock_rep2,/path/to/reads/mock2.fastq.gz,mock,
mock_rep3,/path/to/reads/mock3.fastq.gz,mock,
treated_rep1,/path/to/reads/treat1.fastq.gz,treated,
treated_rep2,/path/to/reads/treat2.fastq.gz,treated,
treated_rep3,/path/to/reads/treat3.fastq.gz,treated,
```

and for paired-end reads, like this:

```csv
Sample,R1,R2,Condition,Patient
mock_rep1,/path/to/reads/mock1_1.fastq,/path/to/reads/mock1_2.fastq,mock,A
mock_rep2,/path/to/reads/mock2_1.fastq,/path/to/reads/mock2_2.fastq,mock,B
mock_rep3,/path/to/reads/mock3_1.fastq,/path/to/reads/mock3_2.fastq,mock,C
treated_rep1,/path/to/reads/treat1_1.fastq,/path/to/reads/treat1_2.fastq,treated,A
treated_rep2,/path/to/reads/treat2_1.fastq,/path/to/reads/treat2_2.fastq,treated,B
treated_rep3,/path/to/reads/treat3_1.fastq,/path/to/reads/treat3_2.fastq,treated,C
```

Read files can be compressed (`.gz`). You need at least two replicates for each condition to run the pipeline. Patient labels are optional and can be used to connect samples belonging to the same patient (or animal, origin, ...) for improved differential expression testing.

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

You can add a [build-in species](#build-in-species) to your defined genomes and annotation with `--species XXX --include_species`.

`--species` is also an identifier for the downstream pathway analysis. Available are WebGestalt set enrichment analysis (GSEA) for `hsa`, piano GSEA with different settings and consensus scoring for `hsa`, `mmu` and `mau`.

#### Build-in species

We provide a small set of build-in species for which the genome and annotation files are automatically downloaded from [Ensembl](https://www.ensembl.org/index.html) with `--species XXX --include_species`. Please let us know, we can easily add other species.

| Species      | three-letter shortcut | Genome                              | Annotation                                    |
| ------------ | --------------------- | ----------------------------------- | --------------------------------------------- |
| Homo sapiens | `hsa`                 | Homo_sapiens.GRCh38.98              | Homo_sapiens.GRCh38.dna.primary_assembly      |
| Mus musculus | `mmu`                 | Mus_musculus.GRCm38.99              | Mus_musculus.GRCm38.dna.primary_assembly      |
| Homo sapiens | `mau`                 | Mesocricetus_auratus.MesAur1.0.100  | Mesocricetus_auratus.MesAur1.0.dna.toplevel   |
| Homo sapiens | `eco`                 | Escherichia_coli_k_12.ASM80076v1.45 | Escherichia_coli_k_12.ASM80076v1.dna.toplevel |

#### Comparisons for DEG analysis

Per default, all possible pairwise comparisons _in one direction_ are performed. Thus, when _A_ is compared against _B_ the pipeline will not automatically compare _B_ vs. _A_ which will anyway only change the direction of the finally resulting fold changes. To change this, please define the needed comparison with `--deg comparisons.csv`, where each line contains a pairwise comparison:

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

### DEG analysis

```bash
--strand                        # strandness for counting with featureCounts: 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default 0]
--tpm                           # threshold for TPM (transcripts per million) filter [default 1]
--deg                           # a CSV file following the pattern: conditionX,conditionY
```

### Transcriptome assembly

```bash
--assemly                       # switch to transcriptome assemly
--busco_db                      # BUSCO database ['euarchontoglires_odb9']
--dammit_uniref90               # add UniRef90 to dammit databases, takes long [false]
```

## Profiles/configuration options

Per default, the pipeline is locally executed with `conda` dependency management (corresponds to `-profile local,conda`). Adjust this setting by combining an _executer_ option with an _engine_ option, e.g. `-profile local,conda` or `-profile slurm,conda`.

### Executor options...

*... or how to schedule your workload.*

Currently implemented are `local`, `slurm` and `lsf` executions.

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

## Output

The result folder is structured by each step and tool (`results/step/tool`) as follows:

```
results/
├── 01-Trimming
│   └── fastp                   trimmed reads
├── 02-rRNARemoval
│   └── SortMeRNA               rRNA-free (and trimmed) reads
├── 03-Mapping
│   └── HISAT2                  mapping results in BAM format with index files (BAI)
├── 04-Counting
│   └── featureCounts           counting table
├── 05-CountingFilter
│   └── TPM                     counting table with additional TPM value; formatted counting table filtered by TPM
├── 06-Annotation               filtered annotation; gene id, name and type mapping
├── 07-DifferentialExpression
│   └── DESeq2                  see below
├── 08-Assembly
│   └── de_novo
│      └── Trinity              Trinity assembly  (with --assembly)
├── 09-RNA-Seq_Annotation       BUSCO and dammit results (with --assembly)
├── Logs                        Nextflow execution timeline and workflow report
└── Summary                     MultiQC report
```

Please note, that `08-Assembly` and `09-RNA-Seq_Annotation` are part of the transcriptome assembly branch (`--assembly`). Here, steps `04` to `07` are not applicable.

### DESeq2 results

The `DESeq2` result is structured as follows:

```
07-DifferentialExpression/
└── DESeq2
   ├── data                         
   │   ├── counts                   normalized, transformed counts; size factors table
   │   └── input                    DESeq2 input summary
   ├── deseq2.Rout                  R log file
   ├── MAQCA_vs_MAQCB               results for pairwise comparison
   │   ├── downstream_analysis  
   │   │   ├── piano                piano results
   │   │   └── WebGestalt           WebGestalt results
   │   ├── input                    DESeq2 input summary
   │   ├── plots
   │   │   ├── heatmaps
   │   │   ├── MA
   │   │   ├── PCA
   │   │   ├── sample2sample
   │   │   └── volcano
   │   ├── reports                  DESeq2 result HTML table; summary report
   │   └── results                  raw and filtered DESeq2 result in CSV and XLSX format; DEG analysis summary
   └── plots                        heatmaps and PCA of all samples
```

We provide `DESeq2` normalized, regularized log (rlog), variance stabilizing (vsd) and log2(n+1) (ntd) transformed count tables (`DESeq2/data/counts`).

For each comparison (specified with `--deg` or per default all possible in one direction) a new folder `X_vs_Y` is created. This also describes the direction of the comparison, e.g. the log2FoldChange describes the change of a gene under condition B with respect to the gene under condition A.

Downstream analysis are currently provided for some species: GSEA consensus scoring with `piano` for *Homo sapiens*, *Mus musculus* and *Mesocricetus auratus*; and `WebGestalt` GSEA *Homo sapiens* and *Mus musculus*.

## Help message

<details><summary>click here to see the complete help message</summary>

```
Usage example:
nextflow run hoelzer-lab/rnaseq --cores 4 --reads input.csv --species eco
or
nextflow run hoelzer-lab/rnaseq --cores 4 --reads input.csv --species eco --assembly
or
nextflow run hoelzer-lab/rnaseq --cores 4 --reads input.csv --genome fasta_virus.csv --annotation gtf_virus.csv --species hsa --include_species
Genomes and annotations from --species, if --include_species is set, --genome and --annotation are concatenated.

Input:
--reads                  a CSV file following the pattern: Sample,R,Condition,Patient for single-end or Sample,R1,R2,Condition,Patient for paired-end
                                    (check terminal output if correctly assigned)
                                    In default all possible comparisons of conditions in one direction are made. Use --deg to change this.
--species                specifies the species identifier for downstream path analysis.
                         If `--include_species` is set, reference genome and annotation are added and automatically downloaded. [default hsa]
                                    Currently supported are:
                                    - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                    - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                    - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]
                                    - mau [Ensembl: Mesocricetus_auratus.MesAur1.0.dna.toplevel | Mesocricetus_auratus.MesAur1.0.100]
--genome                 CSV file with genome reference FASTA files (one path in each line)
                                    If set, --annotation must also be set.
--annotation             CSV file with genome annotation GTF files (one path in each line)
--include_species        Use genome and annotation of supproted species in addition to --genome and --annotation [default false]

Preprocessing options:
--mode                   either 'single' (single-end) or 'paired' (paired-end) sequencing [default single]
--skip_sortmerna         skip rRNA removal via SortMeRNA [default false] 
--index                  the path to the hisat2 index prefix matching the genome provided via --species. 
                         If provided, no new index will be build. Must be named 'index.*.ht2'.  
                         Simply provide the path like 'data/db/index'. DEPRECATED

DEG analysis options:
--strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default 0]
--tpm                    threshold for TPM (transcripts per million) filter. A feature is discared, 
                         if in all conditions the mean TPM value of all libraries in this condition are below the threshold. [default 1]
--deg                    a CSV file following the pattern: conditionX,conditionY
                         Each line stands for one differential gene expression comparison.    

Transcriptome assembly options:
--assembly               perform de novo and reference-based transcriptome assembly instead of DEG analysis [default false]
--busco_db               the database used with BUSCO [default: euarchontoglires_odb9]
                         full list of available data sets at https://busco.ezlab.org/v2/frame_wget.html 
--dammit_uniref90        add UniRef90 to the dammit databases  [default: false]

Computing options:
--cores                  max cores per process for local use [default 10]
--max_cores              max cores used on the machine for local use [default 30]
--memory                 max memory in GB for local use [default 8 GB]
--output                 name of the result folder [default results]

--permanentCacheDir      location for auto-download data like databases [default nextflow-autodownload-databases]
--condaCacheDir          location for storing the conda environments [default conda]
--workdir                working directory for all intermediate results [default /tmp/nextflow-work-ji57pog]
--softlink_results       softlink result files instead of copying

Nextflow options:
-with-report rep.html    cpu / ram usage (may cause errors)
-with-dag chart.html     generates a flowchart for the process tree
-with-timeline time.html timeline (may cause errors)

Execution/Engine profiles:
 The pipeline supports profiles to run via different Executers and Engines e.g.:
 -profile local,conda
  Executer (choose one):
  local
  slurm
  lsf
  Engines (choose one):
  conda
  docker [not supported yet]
  singularity [not supported yet]

For a test run (~ 1 h), add "test" to the profile, e.g. -profile test,local,conda.
Per default: local,conda is executed.

We also provide some pre-configured profiles for certain HPC environments:    
  ara (slurm, conda and parameter customization)
```

</details>

<!-- ## Help message
```
``` -->

<!-- # Flow chart

![flow-chart](figures/chart.png) -->
