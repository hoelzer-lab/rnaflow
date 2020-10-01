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
println "Starting time: $workflow.start"
println "Workdir location:"
println "  $workflow.workDir"
println "Launchdir location:"
println "  $workflow.launchDir"
println "Permanent cache directory:"
println "  $params.permanentCacheDir"
println "Configuration files:"
println "  $workflow.configFiles\u001B[0m"
println " "
if (workflow.profile == 'standard' || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPS to use: $params.max_cores\u001B[0m"
    println " "
}

if (params.assembly) {
    println "\u001B[32mPerform assembly (de novo and reference-based) instead of gene expression analysis."
    if (params.uniref90) {
        params.uniref90_dir = 'uniref90'
        println "Use UniRef90 instead of UniRefKB for annotation: yes\033[0m"
        println " "
    } else {
        println "Use UniRef90 instead of UniRefKB for annotation: no\033[0m"
        println " "
    }
}

Set species = ['hsa', 'eco', 'mmu', 'mau']

if ( params.profile ) { exit 1, "--profile is WRONG use -profile" }
if ( params.reads == '' ) { exit 1, "--reads is a required parameter" }
if ( params.species == '' && params.genome == '' ) { exit 1, "You need to set a genome for mapping and a annotation for counting: with --species " + species + " are provided and automatically downloaded; with --genome and --annotation set csv files for custom input." }
if ( (params.genome && params.annotation == '') || (params.genome == '' && params.annotation) ) { exit 1, "You need to provide genomes AND annotations (--genome and --annotation)." }
if ( params.species && ! (params.species in species) ) { exit 1, "Unsupported species. Use --species with " + species + " or --genome and --annotation." }
//if (params.reference == '') {exit 1, "--reference is a required parameter"}
//if (params.annotation == '') {exit 1, "--annotation is a required parameter"}

if ( params.dge ) { comparison = params.dge } else { comparison = 'all' }
log.info """\
    R N A - S E Q  A S S E M B L Y  &  D I F F E R E N T I A L  G E N E  E X P R E S S I O N  A N A L Y S I S
    = = = = = = =  = = = = = = = =  =  = = = = = = = = = = = =  = = = =  = = = = = = = = = =  = = = = = = = =
    Output path:          $params.output
    Mode:                 $params.mode
    Strandness:           $params.strand
    TPM threshold:        $params.tpm
    Comparisons:          $comparison 
    """
    .stripIndent()

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
* read in auto genome(s)
*/
if ( params.species ) {
    species_auto_ch = Channel.value( params.species )
} else {
    species_auto_ch = Channel.empty()
}

/*
* read in genome(s)
*/
if ( params.genome ) {
    reference_custom_ch = Channel
        .fromPath( params.genome, checkIfExists: true )
        .splitCsv()
        .map { it -> file("${it[0]}", checkIfExists: true) }
} else {
    reference_custom_ch = Channel.empty()
}

/*
* read in annotation(s)
*/
if ( params.annotation ) {
    annotation_custom_ch = Channel
        .fromPath( params.annotation, checkIfExists: true )
        .splitCsv()
        .map { it -> file("${it[0]}", checkIfExists: true) }
} else {
    annotation_custom_ch = Channel.empty()
}


annotated_reads
    .map{row -> row[-2]}
    .toList()
    .toSet()
    .into { sample_conditions; sample_conditions2 }
/*
* read in comparisons
*/
if (params.dge) {
    dge_comparisons_input_ch = Channel
        .fromPath( params.dge, checkIfExists: true )
        .splitCsv( header: true, sep: ',' )
        .map{ row ->
            def condition1 = row['Condition1']
            def condition2 = row['Condition2']
            return [ condition1, condition2 ]
        } // no further processing, in case other tools need this formatted in another way

        // check if conditions form dge and reads file match
        dge_comparisons_input_ch
            .collect()
            .flatten()
            .combine(sample_conditions)
            .subscribe onNext: {
                assert it[1].contains(it[0])
            }, onError: { exit 1, "The comparisons from ${params.dge} do not match the sample conditions in ${params.reads}." }
} else {
    // automatically use all possible comparisons
    dge_comparisons_input_ch = sample_conditions
        .flatten()
        .combine(sample_conditions2.flatten())
        .filter{ it[0] != it[1] }
        .map{ it -> it.sort() }
        .unique()
}
/*
* DESeq2 scripts
*/
deseq2_script = Channel.fromPath( workflow.projectDir + '/bin/deseq2.R', checkIfExists: true )
deseq2_script_refactor_reportingtools_table = Channel.fromPath( workflow.projectDir + '/bin/refactor_reportingtools_table.rb', checkIfExists: true )
deseq2_script_improve_deseq_table = Channel.fromPath( workflow.projectDir + '/bin/improve_deseq_table.rb', checkIfExists: true )

/*
* MultiQC config
*/
multiqc_config = Channel.fromPath( workflow.projectDir + '/assets/multiqc_config.yaml', checkIfExists: true )
regionReport_config = Channel.fromPath( workflow.projectDir + '/assets/regionReport_DESeq2Exploration_custom.Rmd', checkIfExists: true )

//if (params.index) {
//  index_ch = Channel.fromPath("${params.index}.*", checkIfExists: true)
//}

/*
* CHECK INPUT
*/

if (params.assembly) {
    annotated_reads
        .map{ row -> row[-2]}
        .collect()
        .subscribe onNext: {
            for ( i in it ){
            assert 0 <= it.count(i)
            }
        }, onError: { exit 1, 'You need at least one sample to perform an assembly.' }
} else {
    annotated_reads
        .map{ row -> row[-2]}
        .collect()
        .subscribe onNext: {
            for ( i in it ){
            assert 2 <= it.count(i)
            }
        }, onError: { exit 1, 'You need at least 2 samples per condition to perform a differential gene expression analysis.' }
}

if ( ! (params.tpm instanceof java.lang.Double || params.tpm instanceof java.lang.Float || params.tpm instanceof java.lang.Integer) ) {
    exit 1, "--tpm has to be numeric"
}

/************************** 
* MODULES
**************************/

// databases
include {referenceGet; concat_genome} from './modules/referenceGet'
include {annotationGet; concat_annotation} from './modules/annotationGet'
include {sortmernaGet} from './modules/sortmernaGet'
include {hisat2index} from './modules/hisat2'
include {buscoGetDB} from './modules/buscoGetDB'
include {dammitGetDB} from './modules/dammitGetDB'

// analysis
include {fastp} from './modules/fastp'
include {sortmerna} from './modules/sortmerna'
include {hisat2} from './modules/hisat2'
include {featurecounts} from './modules/featurecounts'
include {tpm_filter} from './modules/tpm_filter'
include {deseq2} from './modules/deseq2'
include {fastqc as fastqcPre; fastqc as fastqcPost} from './modules/fastqc'
include {multiqc; multiqc_sample_names} from './modules/multiqc'

// assembly & annotation
include {trinity} from './modules/trinity'
include {busco} from './modules/busco'
include {dammit} from './modules/dammit'
include {stringtie; stringtie_merge} from './modules/stringtie' 

// helpers
include {format_annotation; format_annotation_gene_rows} from './modules/prepare_annotation'

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_auto_reference {
    main:
        // local storage via storeDir
        if (!params.cloudProcess) { referenceGet( species_auto_ch ); reference_auto_ch = referenceGet.out }
        // cloud storage file.exists()?
        if (params.cloudProcess) {
            reference_preload = file("${params.permanentCacheDir}/genomes/${params.species}.fa")
            if (reference_preload.exists()) { reference_auto_ch = Channel.fromPath(reference_preload) }
            else { referenceGet( species_auto_ch ); reference_auto_ch = referenceGet.out } 
        }
    emit:
        reference_auto_ch
}

workflow download_auto_annotation {
    main:
        // local storage via storeDir
        if (!params.cloudProcess) { annotationGet( species_auto_ch ); annotation_auto_ch = annotationGet.out }
        // cloud storage file.exists()?
        if (params.cloudProcess) {
            annotation_preload = file("${params.permanentCacheDir}/annotations/${params.species}.gtf")
            if (annotation_preload.exists()) { annotation_auto_ch = Channel.fromPath(annotation_preload) }
            else { annotationGet( species_auto_ch ); annotation_auto_ch = annotationGet.out } 
        }
    emit:
        annotation_auto_ch
}

workflow download_sortmerna {
    main:
        // local storage via storeDir
        if (!params.cloudProcess) { sortmernaGet(); sortmerna = sortmernaGet.out }
        // cloud storage file.exists()?
        if (params.cloudProcess) {
            sortmerna_preload = file("${params.permanentCacheDir}/databases/sortmerna/data/rRNA_databases")
            if (sortmerna_preload.exists()) { sortmerna = sortmerna_preload }
            else { sortmernaGet(); sortmerna = sortmernaGet.out } 
        }
    emit:
        sortmerna
}

workflow download_busco {
    main:
        if (!params.cloudProcess) { buscoGetDB() ; database_busco = buscoGetDB.out }
        else if (params.cloudProcess) { 
            busco_db_preload = file("${params.permanentCacheDir}/databases/busco/${params.busco}/${params.busco}.tar.gz")
            if (busco_db_preload.exists()) { database_busco = busco_db_preload }
            else  { buscoGetDB(); database_busco = buscoGetDB.out }
        }
    emit: database_busco
}

workflow download_dammit {
    take: 
    busco_db_ch
    
    main:
    dammit_db_preload_path = "${params.permanentCacheDir}/databases/dammit/${params.busco}/dbs.tar.gz"
    if (params.uniref90) {
        dammit_db_preload_path = "${params.permanentCacheDir}/databases/dammit/uniref90/${params.busco}/dbs.tar.gz"
    }
    if (!params.cloudProcess) { dammitGetDB(busco_db_ch) ; database_dammit = dammitGetDB.out }
    if (params.cloudProcess) { 
        dammit_db_preload = file(dammit_db_preload_path)
        if (dammit_db_preload.exists()) { database_dammit = dammit_db_preload }
        else  { dammitGetDB(busco_db_ch); database_dammit = dammitGetDB.out }
    }
    emit: database_dammit
}


/************************** 
* SUB WORKFLOWS
**************************/

/***************************************
Preprocess RNA-Seq reads: qc, trimming, adapters, rRNA-removal, mapping
*/
workflow preprocess {
    take:
        illumina_input_ch
        reference
        sortmerna_db

    main:
        // initial QC of raw reads
        fastqcPre(illumina_input_ch)

        // trim with fastp
        fastp(illumina_input_ch)

        // QC after fastp
        fastqcPost(fastp.out.sample_trimmed)

        if ( params.skip_sortmerna ) {
            sortmerna_no_rna_fastq = fastp.out.sample_trimmed
            sortmerna_log = Channel.empty()
        } else {
            // remove rRNA with SortmeRNA
            sortmerna(fastp.out.sample_trimmed, sortmerna_db)
            sortmerna_no_rna_fastq = sortmerna.out.no_rna_fastq
            sortmerna_log = sortmerna.out.log
        }

        // HISAT2 index
        hisat2index(reference)
        // map with HISAT2
        hisat2(sortmerna_no_rna_fastq, hisat2index.out, params.histat2_additional_params)

    emit:
        sample_bam_ch = hisat2.out.sample_bam
        fastp_json_report = fastp.out.json_report
        sortmerna_log
        hisat2_log = hisat2.out.log  
        fastqcPre = fastqcPre.out.zip  
        fastqcPost = fastqcPost.out.zip
        cleaned_reads_ch = sortmerna_no_rna_fastq
} 


/******************************************
Differential gene expression analysis using a genome reference
*/
workflow expression_reference_based {
    take:
        sample_bam_ch
        fastp_json_report
        sortmerna_log
        hisat2_log
        fastqcPre
        fastqcPost
        annotation
        dge_comparisons_input_ch
        deseq2_script
        deseq2_script_refactor_reportingtools_table
        deseq2_script_improve_deseq_table
        multiqc_config

    main:
        // count with featurecounts
        featurecounts(sample_bam_ch, annotation)

        // prepare annotation for R input
        format_annotation_gene_rows(annotation)
        format_annotation(annotation)

        // filter by TPM value
        // prepare input channels
        tpm_prep_ch = featurecounts.out.counts
            .join( annotated_reads
                    .map{row -> [row[0], row[-2]]} 
            ).toSortedList { entry -> entry[0] } 
            .transpose()

        samples = tpm_prep_ch.first() 
        counts = tpm_prep_ch.take(2).last()
        conditions = tpm_prep_ch.last()

        tpm_filter(samples, counts, conditions)

        // prepare DEseq2
        tpm_filter.out.samples
            .flatMap()
            .join( annotated_reads
                .map{row -> [row[0], row[-2], row[-1]]}
            )
            .multiMap{ it ->
                col_label: it[0]
                condition: it[1]
                patient: it[2]
            }
            .set { annotated_sample }

        deseq2_comparisons = dge_comparisons_input_ch
            .map { it.join(":") }
            .map { "\"${it}\"" }
            .collect()
            .map { it.join(",") }

        // run DEseq2
        deseq2(regionReport_config, tpm_filter.out.filtered_counts, annotated_sample.condition.collect(), 
            annotated_sample.col_label.collect(), deseq2_comparisons, format_annotation.out, format_annotation_gene_rows.out, 
            annotated_sample.patient.collect(), species_auto_ch, deseq2_script, deseq2_script_refactor_reportingtools_table, 
            deseq2_script_improve_deseq_table)

        // run MultiQC
        multiqc_sample_names(annotated_reads.map{ row -> row[0..-3]}.collect())
        multiqc(multiqc_config, 
                multiqc_sample_names.out,
                fastp_json_report.collect(), 
                sortmerna_log.collect().ifEmpty([]), 
                hisat2_log.collect(), 
                featurecounts.out.log.collect(), 
                fastqcPre.collect(),
                fastqcPost.collect(),
                tpm_filter.out.stats,
                params.tpm
        )
} 

/*****************************************
De novo assembly of the preprocessed RNA-Seq reads. For now do co-assembly of all samples. 
ToDo: also do co-assembly of samples belonging to the same condition. 
*/
workflow assembly_denovo {
    take:
        cleaned_reads_ch
        busco_db
        dammit_db

    main:
        reads_ch = cleaned_reads_ch.map {sample, reads -> tuple reads}.collect()
        reads_input_csv = Channel.fromPath( params.reads, checkIfExists: true)
 
        // co-assembly
        trinity(reads_ch, reads_input_csv)

        // qc check
        busco(trinity.out.assembly, busco_db, 'trinity')    

        // transcript annotation 
        dammit(trinity.out.assembly, dammit_db, 'trinity')
} 

/*****************************************
Reference-based assembly of the preprocessed RNA-Seq reads. 
*/
workflow assembly_reference {
    take:
        genome_reference
        annotation_reference
        bams
        busco_db
        dammit_db

    main:
        // StringTie2 GTF-guided transcript prediction
        stringtie(genome_reference, annotation_reference, bams)

        // Merge each single GTF-guided StringTie2 GTF file 
        stringtie_merge(genome_reference, stringtie.out.gtf.collect(), 'stringtie')

        // qc check
        busco(stringtie_merge.out.transcripts, busco_db, 'stringtie')    

        // transcript annotation 
        dammit(stringtie_merge.out.transcripts, dammit_db, 'stringtie')
}


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
    // get the reference genome and index it for hisat2
    download_auto_reference()
    reference_auto = download_auto_reference.out

    // get the annotation
    download_auto_annotation()
    annotation_auto = download_auto_annotation.out

    // concatenate genomes and annotations
    concat_genome(reference_custom_ch.collect().mix(reference_auto).collect())
    reference = concat_genome.out
    concat_annotation(annotation_custom_ch.collect().mix(annotation_auto).collect())
    annotation = concat_annotation.out

    // get sortmerna databases
    download_sortmerna()
    sortmerna_db = download_sortmerna.out

    // preprocess RNA-Seq reads
    preprocess(illumina_input_ch, reference, sortmerna_db)

    // perform assembly & annotation
    if (params.assembly) {
        // dbs
        busco_db = download_busco()
        dammit_db = download_dammit(busco_db)
        // de novo
        assembly_denovo(preprocess.out.cleaned_reads_ch, busco_db, dammit_db)
        // reference-based
        assembly_reference(reference, annotation, preprocess.out.sample_bam_ch, busco_db, dammit_db)
    } else {
    // perform expression analysis
        // start reference-based differential gene expression analysis
        expression_reference_based(preprocess.out.sample_bam_ch,
                                preprocess.out.fastp_json_report,
                                preprocess.out.sortmerna_log,
                                preprocess.out.hisat2_log,
                                preprocess.out.fastqcPre,
                                preprocess.out.fastqcPost,
                                annotation,
                                dge_comparisons_input_ch, 
                                deseq2_script, 
                                deseq2_script_refactor_reportingtools_table, 
                                deseq2_script_improve_deseq_table, 
                                multiqc_config)
    }
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
    nextflow run main.nf --cores 4 --reads input.csv --species eco
    or
    nextflow run main.nf --cores 4 --reads input.csv --species eco --assembly
    or
    nextflow run main.nf --cores 4 --reads input.csv --genome fastas.csv --annotation gtfs.csv
    or
    nextflow run main.nf --cores 4 --reads input.csv --genome fastas.csv --annotation gtfs.csv --species eco
    ${c_dim}Genomes and annotatiosn from --genome, --annotation and --species are concatenated.${c_reset}

    ${c_yellow}Input:${c_reset}
    ${c_green}--reads${c_reset}         a CSV file following the pattern: Sample,R,Condition,Patient for single-end or Sample,R1,R2,Condition,Patient for paired-end
                                        ${c_dim}(check terminal output if correctly assigned)
                                        In default all possible comparisons of conditions are made. Use --dge to change this.${c_reset}
    ${c_green}--species${c_reset}       reference genome and annotation with automatic download.
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]${c_reset}
    ${c_green}--genome${c_reset}        CSV file with genome reference FASTA files (one path in each line).
                                        ${c_dim}If set, --annotation must also be set.${c_reset}
    ${c_green}--annotation${c_reset}    CSV file with genome annotation GTF files (one path in each line)

    ${c_yellow}Options${c_reset}
    --assembly               perform de novo and reference-based transcriptome assembly instead of DEG analysis [default $params.assembly]
    --uniref90               transcriptome annotation using UniRef90 instead of UniRefKB [default $params.uniref90]
    --dge                    a CSV file following the pattern: conditionX,conditionY
                             Each line stands for one differential gene expression comparison.
    --index                  the path to the hisat2 index prefix matching the genome provided via --species. 
                             If provided, no new index will be build. Must be named 'index.*.ht2'.  
                             Simply provide the path like 'data/db/index'. DEPRECATED
    --mode                   either 'single' (single-end) or 'paired' (paired-end) sequencing [default $params.mode]
    --strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default $params.strand]
    --tpm                    threshold for TPM (transcripts per million) filter. A feature is discared, 
                             if in all conditions the mean TPM value of all libraries in this condition are below the threshold. [default $params.tpm]
    --skip_sortmerna         Skip rRNA removal via SortMeRNA [default $params.skip_sortmerna] 
    --busco                 the database used with BUSCO [default: $params.busco]
                ${c_dim}full list of available data sets at https://busco.ezlab.org/v2/frame_wget.html ${c_reset}

    ${c_dim}Computing options:
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
                             ${c_reset}
    """.stripIndent()
}
