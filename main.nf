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
if (params.dge == '') {exit 1, "--dge is a required parameter"}
//if (params.reference == '') {exit 1, "--reference is a required parameter"}
//if (params.annotation == '') {exit 1, "--annotation is a required parameter"}

println "\nD I F F E R E N T I A L  G E N E  E X P R E S S I O N  A N A L Y S I S"
println "= = = = = = = = = = = =  = = = =  = = = = = = = = = =  = = = = = = = ="
println "Reference species:    $params.species"
println "Output path:          $params.output"
println "mode:                 $params.mode"
println "strandness:           $params.strand"

/************************** 
* INPUT CHANNELS 
**************************/

/*
* read in sample sheet
*/
if (params.reads) { 
  if (params.mode == 'single') {
    Channel
      .fromPath( params.reads, checkIfExists: true)
      .splitCsv(header: true, sep: ',')
      .map{row ->
          def sample = row['Sample']
          def read = file(row['R'], checkIfExists: true)
          def condition = row['Condition']
          def patient = row['Patient']
          return [ sample, read, condition, patient ]
      }
      .tap { annotated_reads }
      .tap { illumina_input_ch }
      .map { sample, read, condition, patient ->
              return [ sample, [ read ] ] }
      .set { illumina_input_ch }

  } else {
    Channel
      .fromPath( params.reads, checkIfExists: true)
      .splitCsv(header: true, sep: ',')
      .map{row ->
          def sample = row['Sample']
          def read1 = file(row['R1'], checkIfExists: true)
          def read2 = file(row['R2'], checkIfExists: true)
          def condition = row['Condition']
          def patient = row['Patient']
          return [ sample, read1, read2, condition, patient ]
      }
      .tap { annotated_reads }
      .tap { illumina_input_ch }
      .map { sample, read1, read2, condition, patient ->
              return [ sample, [ read1, read2 ] ] }
      .set { illumina_input_ch }
  }
}

/*
* read in comparisons
*/
if (params.dge) {
  dge_comparisons_input_ch = Channel
      .fromPath( params.dge, checkIfExists: true)
      .splitCsv(header: true, sep: ',')
      .map{row ->
          def condition1 = row['Condition1']
          def condition2 = row['Condition2']
          return [ condition1, condition2 ]
      }
      // no further processing, in case other tools need this formatted in another way
}

//if (params.index) {
//  index_ch = Channel.fromPath("${params.index}.*", checkIfExists: true)
//}

/*
* CHECK INPUT
*/

annotated_reads
    .map{ row -> row[-2]}
    .collect()
    .subscribe onNext: { 
      for ( i in it ){
        assert 2 >= it.count(i)
      }
    }, onError: { exit 1, 'You need at least 2 samples per condition to perform a differential gene expression analysis.' }


sample_conditions = annotated_reads
    .map{row -> row[-2]}
    .toList()
    .toSet()
dge_comparisons_input_ch
    .collect()
    .flatten()
    .combine(sample_conditions)
    .subscribe onNext: {
        assert it[1].contains(it[0])
    }, onError: { exit 1, "The comparisons from ${params.dge} do not match the sample conitions in ${params.reads}." }

/************************** 
* MODULES
**************************/

// databases
include './modules/referenceGet' params(reference: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/annotationGet' params(annotation: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/sortmernaGet' params(cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/hisat2index' params(cores: params.cores, reference: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)

// analysis
include './modules/fastp' params(cores: params.cores, output: params.output, dir: params.fastp_dir, mode: params.mode)
include './modules/sortmerna' params(cores: params.cores, output: params.output, dir: params.sortmerna_dir, mode: params.mode)
include './modules/hisat2' params(cores: params.cores, output: params.output, dir: params.hisat2_dir, mode: params.mode)
include './modules/featurecounts' params(cores: params.cores, output: params.output, dir: params.featurecounts_dir, mode: params.mode, strand: params.strand)
include './modules/deseq2' params(output: params.output, dir: params.deseq2_dir, species: params.species)

// helpers
include './modules/prepare_annotation' params(output: params.output, dir: params.annotation_dir)
include './modules/prepare_annotation_gene_rows' params(output: params.output, dir: params.annotation_dir)

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

workflow download_sortmerna {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { sortmernaGet(); sortmerna = sortmernaGet.out }
    // cloud storage file.exists()?
    if (params.cloudProcess) {
      sortmerna_preload = file("${params.cloudDatabase}/databases/sortmerna/rRNA_databases")
      if (sortmerna_preload.exists()) { sortmerna = sortmerna_preload }
      else  { sortmernaGet(); sortmerna = sortmernaGet.out } 
    }
  emit: sortmerna
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
        hisat2_index
        annotation
        sortmerna_db
        dge_comparisons_input_ch

  main:
    //trim
    fastp(illumina_input_ch)

    //remove rRNA
    sortmerna(fastp.out, sortmerna_db)

    //map
    hisat2(sortmerna.out, hisat2_index)

    //count
    featurecounts(hisat2.out, annotation)

    //prepare annotation for R input
    prepare_annotation_gene_rows(annotation)
    prepare_annotation(annotation)

    //defs
    script = Channel.fromPath( "${params.scripts_dir}/deseq2.R", checkIfExists: true )
    script_refactor_reportingtools_table = Channel.fromPath( "${params.scripts_dir}/refactor_reportingtools_table.rb", checkIfExists: true )
    script_improve_deseq_table = Channel.fromPath( "${params.scripts_dir}/improve_deseq_table.rb", checkIfExists: true )

    deseq2_comparisons = dge_comparisons_input_ch
        .map { it.join(":") }
        .map { "\"${it}\"" }
        .collect()
        .map { it.join(",") }

    featurecounts.out
        .fork{tuple -> 
        sample: tuple[0]
        fc_counts_formated: tuple[1][1]
        }
        .set { fc_out }

    fc_out.fc_counts_formated
        .collect()
        .set { fc_counts_formated }

    fc_out.sample
        .join( annotated_reads )
        .fork{ it ->
        col_label: it[0]
        patient: it[-1]
        condition: it[-2]
        }
        .set { annotated_sample }

    annotated_sample.col_label
        .collect()
        .set { col_labels }

    annotated_sample.condition
        .collect()
        .set { conditions }

    annotated_sample.patient
        .collect{ 
          if (it) {
            return it
          } else {
            return
          }
         }
        .set { patients }

    deseq2(fc_counts_formated, col_labels, conditions, patients, deseq2_comparisons, prepare_annotation.out, prepare_annotation_gene_rows.out, script, script_refactor_reportingtools_table, script_improve_deseq_table)
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
      hisat2_index = hisat2_index_reference.out

      // get the annotation
      download_annotation()
      annotation = download_annotation.out

      // get sortmerna databases
      download_sortmerna()
      sortmerna_db = download_sortmerna.out

      // start reference-based analysis
      analysis_reference_based(illumina_input_ch, hisat2_index, annotation, sortmerna_db, dge_comparisons_input_ch)
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
    nextflow run main.nf --cores 4 --reads input.csv --dge comparisons.csv

    ${c_yellow}Input:${c_reset}
    ${c_green}--reads${c_reset}         a CSV file following the pattern: Sample,R,Condition,Patient for single-end or Sample,R1,R2,Condition,Patient for paired-end
                                        ${c_dim}(check terminal output if correctly assigned)${c_reset}
    ${c_green}--dge${c_reset}           a CSV file following the pattern: conditionX,conditionY
                                        Each line stands for one differential gene expression comparison.
    ${c_green}--species${c_reset}       reference genome and annotation are selected based on this parameter [default $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]${c_reset}

    ${c_yellow}Options${c_reset}
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
