#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*
* RNA-Seq-based detection of differentially expressed genes
*
* Author: martin.hoelzer@uni-jena.de
* Author: marie.lataretu@uni-jena.de
*/

// terminal prints
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.reads == '') {exit 1, "--reads is a required parameter"}
//if (params.reference == '') {exit 1, "--reference is a required parameter"}
//if (params.annotation == '') {exit 1, "--annotation is a required parameter"}

println "\nD I F F E R E N T I A L  G E N E  E X P R E S S I O N  A N A L Y S I S"
println "= = = = = = = = = = = =  = = = =  = = = = = = = = = =  = = = = = = = ="
println "Reference species:    $params.species"
println "SortMeRNA database:   $params.sortmerna_db"
println "Output path:          $params.output"
println "mode:                 $params.mode"
println "strandness:           $params.strand"

/************************** 
* INPUT CHANNELS 
**************************/

// illumina reads input via CSV
if (params.reads) { 
  illumina_input_ch = Channel
                .fromPath( params.reads, checkIfExists: true )
                .splitCsv()
                .map { row -> ["${row[0]}", [file("${row[1]}"), file("${row[2]}")]] }
                //.view() 
}

//if (params.index) {
//  index_ch = Channel.fromPath("${params.index}.*", checkIfExists: true)
//}


/************************** 
* MODULES
**************************/

/* Comment section: */
// This is the new DSL2 syntax that uses modules, exemplarily implemented for the trimming process
include './modules/referenceGet' params(reference: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/annotationGet' params(annotation: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/hisat2index' params(cores: params.cores, reference: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/fastp' params(cores: params.cores, output: params.output, dir: params.fastp_dir, mode: params.mode)
include './modules/hisat2' params(cores: params.cores, output: params.output, dir: params.hisat2_dir, mode: params.mode)
include './modules/featurecounts' params(cores: params.cores, output: params.output, dir: params.featurecounts_dir, mode: params.mode, strand: params.strand)
include './modules/prepare_annotation' params(output: params.output, dir: params.annotation_dir)
include './modules/deseq2' params(output: params.output, dir: params.deseq2_dir, species: params.species)

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_reference {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { referenceGet(); reference = referenceGet.out }
    // cloud storage file.exists()?
    if (params.cloudProcess) {
      reference_preload = file("${params.cloudDatabase}/genomes/${params.reference}/${params.reference}.fa.gz")
      if (reference_preload.exists()) { reference = reference_preload }
      else  { referenceGet(); reference = referenceGet.out } 
    }
  emit: reference
}

workflow download_annotation {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { annotationGet(); annotation = annotationGet.out }
    // cloud storage file.exists()?
    if (params.cloudProcess) {
      annotation_preload = file("${params.cloudDatabase}/annotations/${params.annotation}/${params.annotation}.gtf.gz")
      if (annotation_preload.exists()) { annotation = annotation_preload }
      else  { annotationGet(); annotation = annotationGet.out } 
    }
  emit: annotation
}

workflow hisat2_index_reference {
  get: reference
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { hisat2index(reference); index = hisat2index.out }
    // cloud storage file.exists()?
    if (params.cloudProcess) {
      index_preload = file("${params.cloudDatabase}/genomes/${params.reference}/${params.reference}*.ht2")
      if (index_preload.exists()) { index = index_preload }
      else  { hisat2index(reference); index = hisat2index.out } 
    }
  emit: index
}


/************************** 
* SUB WORKFLOWS
**************************/

/* Comment section: */


workflow analysis_reference_based {
  get:  illumina_input_ch
        index
        annotation

  main:
    //trim
    fastp(illumina_input_ch)

    //map
    hisat2(fastp.out, index)

    //count
    featurecounts(hisat2.out, annotation)

    //prepare annotation for R input
    prepare_annotation(annotation)

    //defs
    featurecounts.out
      .fork{tuple -> 
      name: tuple[0]
      fc_formated: tuple[1][1]
      }
      .set { fc_out }
    script = Channel.fromPath( "${params.scripts_dir}/deseq2.R" )

    fc_out.name
      .collect()
      .set { name }

    fc_out.fc_formated
      .collect()
      .set { fc_formated }

    deseq2(name, fc_formated, prepare_annotation.out, script)
} 

workflow analysis_de_novo {
/*WIP. maybe later...*/
} 


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
      // get the reference genome and index it for hisat2
      hisat2_index_reference(download_reference())
      index = hisat2_index_reference.out

      // get the annotation
      download_annotation()
      annotation = download_annotation.out

      // start reference-based analysis
      analysis_reference_based(illumina_input_ch, index, annotation)
}


def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    ${c_yellow}Usage example:${c_reset}
    nextflow run main.nf --cores 4 --reads input.csv

    ${c_yellow}Input:${c_reset}
    ${c_green}--reads${c_reset}         a CSV file following the pattern: mock_rep1,fastq1,fastq2 with fastq2 beeing optional 
                                        ${c_dim}(check terminal output if correctly assigned)${c_reset}
    ${c_green}--species${c_reset}       reference genome and annotation are selected based on this parameter [default $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]${c_reset}

    ${c_yellow}Options${c_reset}
    --sortmerna              the database used for SortMeRNA
    --index                  the path to the hisat2 index prefix matching the genome provided via --species. 
                             If provided, no new index will be build. Must be named 'index.*.ht2'.  
                             Simply provide the path like 'data/db/index'. DEPRECATED
    --mode                   either 'single' (single-end) or 'paired' (paired-end) sequencing [default $params.mode]
    --strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default $params.strand]
    --cores                  max cores for local use [default $params.cores]
    --memory                 max memory in GB for local use [default $params.memory]
    --output                 name of the result folder [default $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 standard (local, including conda) [default standard]
                             ${c_reset}
    """.stripIndent()
}
