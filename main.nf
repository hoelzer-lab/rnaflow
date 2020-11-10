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
println "  $workflow.configFiles"
println "Cmd line:"
println "  $workflow.commandLine\u001B[0m"
if (workflow.repository != null){ println "\033[2mGit info: $workflow.repository - $workflow.revision [$workflow.commitId]\u001B[0m" }
println " "
if (workflow.profile.contains('standard') || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    println " "
}

if ( !workflow.revision ) { 
    println ""
    println "\033[0;33mWARNING: not a stable execution. Please use -r for full reproducibility.\033[0m\n"
}

def folder = new File(params.output)
if ( folder.exists() ) { 
    println ""
    println "\033[0;33mWARNING: Output folder already exists. Results might be overwritten! You can adjust the output folder via [--output]\033[0m\n"
}


if (params.assembly) {
    println "\u001B[32mPerform assembly (de novo and reference-based) instead of gene expression analysis."
    if (params.dammit_uniref90) {
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
if ( params.include_species && ! params.species ) { exit 1, "You need to select a species with --species " + species + " for automatic download." }
if ( params.include_species == '' && params.genome == '' ) { exit 1, "You need to set a genome for mapping and a annotation for counting: with --include_species --species " + species + " are provided and automatically downloaded; with --genome and --annotation set csv files for custom input." }
if ( (params.genome && params.annotation == '') || (params.genome == '' && params.annotation) ) { exit 1, "You need to provide genomes AND annotations (--genome and --annotation)." }
if ( (params.include_species && params.species) && ! (params.species in species) ) { exit 1, "Unsupported species for automatic download. Suported species as: " + species}

if ( params.deg ) { comparison = params.deg } else { comparison = 'all' }
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
            def read = (workflow.profile.contains('test')) ? file("$workflow.projectDir/" + row['R'], checkIfExists: true) : file(row['R'], checkIfExists: true)
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
            def read1 = (workflow.profile.contains('test')) ? file("$workflow.projectDir/" + row['R1'], checkIfExists: true) : file(row['R1'], checkIfExists: true)
            def read2 = (workflow.profile.contains('test')) ? file("$workflow.projectDir/" + row['R2'], checkIfExists: true) : file(row['R2'], checkIfExists: true)
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
species_auto_ch = Channel.value( params.species )

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

sample_conditions = annotated_reads
    .map{row -> row[-2]}
    .toList()
    .toSet()

/*
* read in comparisons
*/
if (params.deg) {
    deg_comparisons_input_ch = Channel
        .fromPath( params.deg, checkIfExists: true )
        .splitCsv( header: true, sep: ',' )
        .map{ row ->
            def condition1 = row['Condition1']
            def condition2 = row['Condition2']
            return [ condition1, condition2 ]
        } // no further processing, in case other tools need this formatted in another way

        // check if conditions form deg and reads file match
        deg_comparisons_input_ch
            .collect()
            .flatten()
            .combine(sample_conditions)
            .subscribe onNext: {
                assert it[1].contains(it[0])
            }, onError: { exit 1, "The comparisons from ${params.deg} do not match the sample conditions in ${params.reads}." }
} else {
    // automatically use all possible comparisons
    deg_comparisons_input_ch = sample_conditions
        .flatten()
        .combine(sample_conditions.flatten())
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

// test
include {get_reduced_genome_test; reduce_genome_test ; get_reduced_annotation_test ; reduce_annotation_test } from './modules/get_test_data.nf'

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
include {hisat2; index_bam} from './modules/hisat2'
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

workflow get_test_data {
    main:
        reference_test_preload = file("${params.permanentCacheDir}/genomes/${params.species}_small.fa")
        if ( reference_test_preload.exists()) { 
            reference_test = Channel.fromPath(reference_test_preload)
        } else {
            reference_complete_preload = file("${params.permanentCacheDir}/genomes/${params.species}.fa")
            if (reference_complete_preload.exists()) { 
                reference_auto_ch = Channel.fromPath( reference_complete_preload )
                reference_test = reduce_genome_test ( reference_auto_ch )
            } else {
                reference_test = get_reduced_genome_test( species_auto_ch )
            }
        }

        annotation_test_reload = file("${params.permanentCacheDir}/annotations/${params.species}_small.gtf")
        if ( annotation_test_reload.exists() ) {
            annotation_test = Channel.fromPath(annotation_test_reload)
        } else {
            annotation_complete_preload = file("${params.permanentCacheDir}/annotations/${params.species}.gtf")
            if (annotation_complete_preload.exists()) { 
                annotation_auto_ch = Channel.fromPath( annotation_complete_preload )
                annotation_test = reduce_annotation_test( annotation_auto_ch )
            } else { 
                annotation_test = get_reduced_annotation_test( species_auto_ch )
            } 
        }

    emit:
        reference_test
        annotation_test
}

workflow download_auto_reference {
    main:
        if (params.include_species){
            // local storage via storeDir
            if (!params.cloudProcess) { referenceGet( species_auto_ch ); reference_auto_ch = referenceGet.out }
            // cloud storage file.exists()?
            if (params.cloudProcess) {
                reference_preload = file("${params.permanentCacheDir}/genomes/${params.species}.fa")
                if (reference_preload.exists()) { reference_auto_ch = Channel.fromPath(reference_preload) }
                else { referenceGet( species_auto_ch ); reference_auto_ch = referenceGet.out } 
            }
        } else {
            reference_auto_ch = Channel.empty()
        }
    emit:
        reference_auto_ch
}

workflow download_auto_annotation {
    main:
        if (params.include_species){
            // local storage via storeDir
            if (!params.cloudProcess) { annotationGet( species_auto_ch ); annotation_auto_ch = annotationGet.out }
            // cloud storage file.exists()?
            if (params.cloudProcess) {
                annotation_preload = file("${params.permanentCacheDir}/annotations/${params.species}.gtf")
                if (annotation_preload.exists()) { annotation_auto_ch = Channel.fromPath(annotation_preload) }
                else { annotationGet( species_auto_ch ); annotation_auto_ch = annotationGet.out } 
            }
        } else {
            annotation_auto_ch = Channel.empty()
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
            busco_db_preload = file("${params.permanentCacheDir}/databases/busco/${params.busco_db}/${params.busco_db}.tar.gz")
            if (busco_db_preload.exists()) { database_busco = busco_db_preload }
            else  { buscoGetDB(); database_busco = buscoGetDB.out }
        }
    emit: database_busco
}

workflow download_dammit {
    take: 
    busco_db_ch
    
    main:
    dammit_db_preload_path = "${params.permanentCacheDir}/databases/dammit/${params.busco_db}/dbs.tar.gz"
    if (params.dammit_uniref90) {
        dammit_db_preload_path = "${params.permanentCacheDir}/databases/dammit/uniref90/${params.busco_db}/dbs.tar.gz"
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
        // index BAM files
        index_bam(hisat2.out.sample_bam)

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
        deg_comparisons_input_ch
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
        
        deseq2_comparisons = deg_comparisons_input_ch
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
    if ( workflow.profile.contains('test') ){
        get_test_data()
        reference = get_test_data.out.reference_test.collect()
        annotation = get_test_data.out.annotation_test.collect()
    } else {
        // get the reference genome
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
    }

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
                                deg_comparisons_input_ch, 
                                deseq2_script, 
                                deseq2_script_refactor_reportingtools_table, 
                                deseq2_script_improve_deseq_table, 
                                multiqc_config)
    }
}


workflow.onComplete { 
    if (workflow.success) {
        // copy execution and timeline HTML reports to output dir
        println (['bash', "${workflow.projectDir}/bin/reports.sh", "${params.output}", "${workflow.projectDir}/${params.runinfo}"].execute().text)
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
    nextflow run hoelzer-lab/rnaseq --cores 4 --reads input.csv --species eco
    or
    nextflow run hoelzer-lab/rnaseq --cores 4 --reads input.csv --species eco --assembly
    or
    nextflow run hoelzer-lab/rnaseq --cores 4 --reads input.csv --genome fasta_virus.csv --annotation gtf_virus.csv --species hsa --include_species
    ${c_dim}Genomes and annotations from --species, if --include_species is set, --genome and --annotation are concatenated.${c_reset}

    ${c_yellow}Input:${c_reset}
    ${c_green}--reads${c_reset}                  a CSV file following the pattern: Sample,R,Condition,Patient for single-end or Sample,R1,R2,Condition,Patient for paired-end
                                        ${c_dim}(check terminal output if correctly assigned)
                                        In default all possible comparisons of conditions in one direction are made. Use --deg to change this.${c_reset}
    ${c_green}--species${c_reset}                specifies the species identifier for downstream path analysis.
                             If `--include_species` is set, reference genome and annotation are added and automatically downloaded. [default $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]
                                        - mau [Ensembl: Mesocricetus_auratus.MesAur1.0.dna.toplevel | Mesocricetus_auratus.MesAur1.0.100]${c_reset}
    ${c_green}--genome${c_reset}                 CSV file with genome reference FASTA files (one path in each line)
                                        ${c_dim}If set, --annotation must also be set.${c_reset}
    ${c_green}--annotation${c_reset}             CSV file with genome annotation GTF files (one path in each line)
    ${c_green}--include_species${c_reset}        Use genome and annotation of supproted species in addition to --genome and --annotation [default $params.include_species]

    ${c_yellow}Preprocessing options:${c_reset}
    --mode                   either 'single' (single-end) or 'paired' (paired-end) sequencing [default $params.mode]
    --skip_sortmerna         skip rRNA removal via SortMeRNA [default $params.skip_sortmerna] 
    ${c_dim}--index                  the path to the hisat2 index prefix matching the genome provided via --species. 
                             If provided, no new index will be build. Must be named 'index.*.ht2'.  
                             Simply provide the path like 'data/db/index'. DEPRECATED${c_reset}

    ${c_yellow}DEG analysis options:${c_reset}
    --strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default $params.strand]
    --tpm                    threshold for TPM (transcripts per million) filter. A feature is discared, 
                             if in all conditions the mean TPM value of all libraries in this condition are below the threshold. [default $params.tpm]
    --deg                    a CSV file following the pattern: conditionX,conditionY
                             Each line stands for one differential gene expression comparison.    

    ${c_yellow}Transcriptome assembly options:${c_reset}
    --assembly               perform de novo and reference-based transcriptome assembly instead of DEG analysis [default $params.assembly]
    --busco_db               the database used with BUSCO [default: $params.busco_db]
                             ${c_dim}full list of available data sets at https://busco.ezlab.org/v2/frame_wget.html ${c_reset}
    --dammit_uniref90        add UniRef90 to the dammit databases  [default: $params.dammit_uniref90]

    ${c_dim}Computing options:
    --cores                  max cores per process for local use [default $params.cores]
    --max_cores              max cores used on the machine for local use [default $params.max_cores]
    --memory                 max memory in GB for local use [default $params.memory]
    --output                 name of the result folder [default $params.output]

    --permanentCacheDir      location for auto-download data like databases [default $params.permanentCacheDir]
    --condaCacheDir          location for storing the conda environments [default $params.condaCacheDir]
    --workdir                working directory for all intermediate results [default $params.workdir]

    Nextflow options:
    -with-tower              Activate monitoring via Nextflow Tower (needs TOWER_ACCESS_TOKEN set)
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}Execution/Engine profiles:${c_reset}
     The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.:
     -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
      ${c_green}Executer${c_reset} (choose one):
      local
      slurm
      lsf
      ${c_blue}Engines${c_reset} (choose one):
      conda
      docker [not supported yet]
      singularity [not supported yet]
    
    For a test run (~ 1 h), add "test" to the profile, e.g. -profile test,local,conda.
    Per default: local,conda is executed. 

    We also provide some pre-configured profiles for certain HPC environments:    
      ara (slurm, conda and parameter customization)
    ${c_reset}
    """.stripIndent()
}
