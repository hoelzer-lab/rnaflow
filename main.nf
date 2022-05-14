#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
* RNA-Seq-based detection of differentially expressed genes
*
* Authors: marie.lataretu@uni-jena.de, fischerd@rki.de, hoelzer.martin@gmail.com
*/

// Parameters sanity checking

Set valid_params = ['max_cores', 'cores', 'memory', 'profile', 'help', 'reads', 'genome', 'nanopore', 'minimap2_additional_params', 'minimap2_dir',  'annotation', 'deg', 'autodownload', 'pathway', 'species', 'include_species', 'strand', 'mode', 'tpm', 'fastp_additional_params', 'hisat2_additional_params', 'featurecounts_additional_params', 'feature_id_type', 'busco_db', 'dammit_uniref90', 'skip_sortmerna', 'skip_read_preprocessing', 'assembly', 'output', 'fastp_dir', 'sortmerna_dir', 'hisat2_dir', 'featurecounts_dir', 'tpm_filter_dir', 'annotation_dir', 'deseq2_dir', 'assembly_dir', 'rnaseq_annotation_dir', 'uniref90_dir', 'readqc_dir', 'multiqc_dir', 'nf_runinfo_dir', 'permanentCacheDir', 'condaCacheDir', 'singularityCacheDir', 'softlink_results', 'cloudProcess', 'permanent-cache-dir', 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process', 'setup', 'rna'] // don't ask me why there is 'permanent-cache-dir', 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process'

def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

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
if ( workflow.profile.contains('singularity') ) {
    println "Singularity cache directory:"
    println "  $params.singularityCacheDir"
}
if ( workflow.profile.contains('conda') ) { 
    println "Conda cache directory:"
    println "  $params.condaCacheDir"
}
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

if ( workflow.profile.contains('singularity') ) {
    println ""
    println "\033[0;33mWARNING: Singularity image building sometimes fails!"
    println "Multiple resumes (-resume) and --max_cores 1 --cores 1 for local execution might help.\033[0m\n"
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

if (params.nanopore) {
    println "\u001B[32mPerform processing of reads in Nanopore mode instead default short-read mode. After mapping, the same steps are used as for Illumina.\033[0m"
}


Set species = ['hsa', 'eco', 'mmu', 'mau']
Set autodownload = ['hsa', 'eco', 'mmu', 'mau']
Set pathway = ['hsa', 'mmu', 'mau']

if ( params.profile ) { exit 1, "--profile is WRONG use -profile" }

// required stuff
if ( ! params.reads && ! params.setup ) { exit 1, "--reads is a required parameter" }

// deprecated stuff
if ( params.mode ) { println "\033[0;33mWARNING: Parameter --mode is deprecated, read mode will automatically be detected from the sample file.\033[0m\n" }
if ( ((params.species || params.include_species) && (params.autodownload || params.pathway)) && ! workflow.profile.contains('test') ) {  exit 1, "Please use '--autodownload " + autodownload + " --pathway " + pathway + "' OR '--species " + species + " --include_species' (deprecated)." }
if ( ( params.species || params.include_species ) && ! workflow.profile.contains('test') ) { 
    println "\033[0;33mWARNING: --species " + species + " and --include_species are deprecated parameters. Please use --autodownload " + autodownload + " (corresponds to '--species " + species + " --include_species') and --pathway " + pathway + " (corresponds to '--species " + species + "') in the future.\033[0m\n" 
    if ( params.include_species && ! params.species ) { exit 1, "You need to select a species with --species " + species + " for automatic download." }
    if ( (! params.include_species) && params.genome == '' && ! workflow.profile.contains('test') ) { exit 1, "You need to select a species with --species " + species + " for automatic download." }
    if ( (params.genome && params.annotation == '') || (params.genome == '' && params.annotation) ) { exit 1, "You need to provide genomes AND annotations (--genome and --annotation)." }
    if ( (params.include_species && params.species) && ! params.species in species ) { exit 1, "Unsupported species for automatic download. Supported species are: " + species}
} else {
    if ( ! params.autodownload && ! params.genome && ! workflow.profile.contains('test') && !params.setup ) { exit 1, "You need to set a genome for mapping and an annotation for counting: with --autodownload " + autodownload + " are provided and automatically downloaded; with --genome and --annotation set csv files for custom input." }
    // logic stuff
    if ( params.genome && ! params.annotation ) { exit 1, "You need to provide genomes AND annotations (--genome and --annotation)." }
    if ( ! params.autodownload in autodownload ) { exit 1, "Unsupported species for automatic download. Supported species are: " + autodownload }
    if ( ! params.pathway in pathway ) { exit 1, "Unsupported species for downstream pathway analysis. Supported species are: " + pathway }
}

if ( params.deg ) { comparison = params.deg } else { comparison = 'all' }

/************************** 
* INPUT CHANNELS 
**************************/
if (params.setup) {
    Channel.fromPath( './configs/container.config' )
            .splitCsv(skip: 1, sep: '\t')
            .map{ row ->
                    if ( row[1] != null && row[2] != null) {
                        def tool = row[1]
                        def path = row[2].split('"')[1]
                        return [tool, path] 
                    }
            }
            .tap{ container_ch }
}


if (params.reads) { 
    Channel
        .fromPath( params.reads, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map{row ->
            def paired_end = row['R2'] ? true : false
            def read1 = row['R'] ? (workflow.profile.contains('test')) ? file("$workflow.projectDir/" + row['R'], checkIfExists: true) : file(row['R'], checkIfExists: true) : (workflow.profile.contains('test')) ? file("$workflow.projectDir/" + row['R1'], checkIfExists: true) : file(row['R1'], checkIfExists: true)
            def read2 = paired_end ? (workflow.profile.contains('test')) ? file("$workflow.projectDir/" + row['R2'], checkIfExists: true) : file(row['R2']) : "" 
            def meta = [:]
                meta.sample = row['Sample']
                meta.condition = row['Condition']
                meta.source = row['Source']
                meta.paired_end = paired_end
                meta.strandedness = row['Strandedness'] ? row['Strandedness'] : params.strand
            return meta.paired_end ? [ meta, [ read1, read2 ] ] : [ meta, [ read1 ] ]
        }
        .tap { annotated_reads }
        .tap { read_input_ch }
        .set { read_input_ch }
}else if ((!params.reads && params.setup) || params.setup) {
    println "\u001B[32mRunning in setup mode. Only necessary database and reference files will be downloaded.\033[0m\n"
    annotated_reads = Channel.empty()
    read_input_ch = Channel.empty()
}else{
    exit 1, "Parameter 'reads' undefined."
}

param_strand = ""
param_read_mode = ""

if (!params.setup) { 
    File csvFile = new File(params.reads)
    csvFile.eachLine { line ->
        def row = line.split(",")
        param_strand = row[5] ? row[5] : params.strand
        param_read_mode = row[2] ? "paired-end" : "single-end"
    }

    if ( param_strand == "0" ) { param_strand = "unstranded" }else if ( param_strand == "1" ) { param_strand = "stranded" }else if( param_strand == "2" ){ param_strand = "reversly stranded" }else{exit 1, "Could not detect strandedness of input file. Invalid strandedness parameter ${param_strand}."}
}

log.info """\
                R N A F L O W : R N A - S E Q  A S S E M B L Y  &  D I F F E R E N T I A L  G E N E  E X P R E S S I O N  A N A L Y S I S
                = = = = = = =   = = = = = = =  = = = = = = = =  =  = = = = = = = = = = = =  = = = =  = = = = = = = = = =  = = = = = = = =
                Output path:                    $params.output
                Strandedness                    $param_strand
                Read mode:                      $param_read_mode
                TPM threshold:                  $params.tpm
                Comparisons:                    $comparison 
                Nanopore mode:                  $params.nanopore
                """
                .stripIndent()  
/*
* read in autodownload genome(s)
*/
if ( params.species ) {
    species_auto_ch = Channel.value( params.species ) // deprecated reminder
} else if ( params.autodownload ) { 
    species_auto_ch = Channel.value( params.autodownload )
} else {
    species_auto_ch = Channel.value( '' )
}

/*
* read in species for downstream pathway analysis
*/
if ( params.species ) {
    species_pathway_ch = Channel.value( params.species ) // deprecated reminder
} else if ( params.pathway ) {
    species_pathway_ch = Channel.value( params.pathway )
} else {
    species_pathway_ch = Channel.value( '' )
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

sample_conditions = annotated_reads
    .map{ meta, reads -> meta.condition }
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
* FeatureCounts params
*/
if ( params.featurecounts_additional_params.contains('-t ') ) {
    for (param in params.featurecounts_additional_params.split('-')) {
        if (param.startsWith('t ')) {
            gtf_feature_type_ch = Channel.value( param.stripIndent().split(' ')[-1] )
        }
    }
} else {
    gtf_feature_type_ch = Channel.value('exon')
}
if ( params.featurecounts_additional_params.contains('-g ') ){
    for (param in params.featurecounts_additional_params.split('-')) {
        if (param.startsWith('g ')){
            value_of_g = param.stripIndent().split(' ')[-1]
            gtf_attr_type_ch = Channel.value( value_of_g )
            gtf_feature_type_of_attr_type_ch = Channel.value( value_of_g.split('_')[0] )
        }
    }
} else {
    gtf_attr_type_ch = Channel.value('gene_id')
    gtf_feature_type_of_attr_type_ch = Channel.value('gene')
}

/*
* DESeq2
*/
deseq2_script = Channel.fromPath( workflow.projectDir + '/bin/deseq2.R', checkIfExists: true )
deseq2_script_refactor_reportingtools_table = Channel.fromPath( workflow.projectDir + '/bin/refactor_reportingtools_table.rb', checkIfExists: true )
deseq2_script_improve_deseq_table = Channel.fromPath( workflow.projectDir + '/bin/improve_deseq_table.rb', checkIfExists: true )
deseq2_id_type_ch = Channel.value(params.feature_id_type)
species2prefix = Channel.fromPath( workflow.projectDir + '/assets/ens_species_mapping.tsv', checkIfExists: true)

/*
* MultiQC config
*/
multiqc_config = Channel.fromPath( workflow.projectDir + '/assets/multiqc_config.yaml', checkIfExists: true )
regionReport_config = Channel.fromPath( workflow.projectDir + '/assets/regionReport_DESeq2Exploration_custom.Rmd', checkIfExists: true )

/*
* CHECK INPUT
*/

if (params.assembly) {
    annotated_reads
        .map{ meta, reads -> meta.condition }
        .collect()
        .subscribe onNext: {
            for ( i in it ){
            assert 0 <= it.count(i)
            }
        }, onError: { exit 1, 'You need at least one sample to perform an assembly.' }
} else {
    annotated_reads
        .map{ meta, reads -> meta.condition }
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
include {minimap2index} from './modules/minimap2'
include {dammitGetDB} from './modules/dammitGetDB'
include {buscoGetDB} from './modules/buscoGetDB'

// analysis
include {fastp} from './modules/fastp'
include {sortmerna} from './modules/sortmerna'
include {hisat2; index_bam} from './modules/hisat2'
include {minimap2; index_bam_minimap2} from './modules/minimap2'
include {featurecounts} from './modules/featurecounts'
include {tpm_filter} from './modules/tpm_filter'
include {deseq2} from './modules/deseq2'
include {fastqc as fastqcPre; fastqc as fastqcPost} from './modules/fastqc'
include {nanoplot as nanoplot} from './modules/nanoplot'
include {multiqc; multiqc_sample_names} from './modules/multiqc'

// assembly & annotation
include {trinity} from './modules/trinity'
include {busco} from './modules/busco'
include {dammit} from './modules/dammit'
include {stringtie; stringtie_merge} from './modules/stringtie' 
include {rattle} from './modules/rattle'

// helpers
include {format_annotation; format_annotation_gene_rows} from './modules/prepare_annotation'
include {containerGet} from './modules/containerGet'
include {extract_tar_bz2} from './modules/utils'
include {format_read_name} from './modules/format_read_name'

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow get_test_data {
    main:
        reference_test_preload = file("${params.permanentCacheDir}/genomes/${params.species}_small.fa") // deprecated reminder
        if ( reference_test_preload.exists()) { 
            reference_test = Channel.fromPath(reference_test_preload)
        } else {
            reference_complete_preload = file("${params.permanentCacheDir}/genomes/${params.species}.fa") // deprecated reminder
            if (reference_complete_preload.exists()) { 
                reference_auto_ch = Channel.fromPath( reference_complete_preload )
                reference_test = reduce_genome_test ( reference_auto_ch )
            } else {
                reference_test = get_reduced_genome_test( species_auto_ch )
            }
        }

        annotation_test_reload = file("${params.permanentCacheDir}/annotations/${params.species}_small.gtf") // deprecated reminder
        if ( annotation_test_reload.exists() ) {
            annotation_test = Channel.fromPath(annotation_test_reload)
        } else {
            annotation_complete_preload = file("${params.permanentCacheDir}/annotations/${params.species}.gtf") // deprecated reminder
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
    take:
        setup_ch
    main:
        if (params.autodownload || params.include_species){ // deprecated reminder
            // local storage via storeDir
            if (!params.cloudProcess) { referenceGet( species_auto_ch ); reference_auto_ch = referenceGet.out }
            // cloud storage file.exists()?
            if (params.cloudProcess) {
                reference_preload = file("${params.permanentCacheDir}/genomes/${params.species}.fa") // deprecated reminder
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
    take:
        setup_ch
    main:
        if (params.autodownload || params.include_species){ // deprecated reminder
            // local storage via storeDir
            if (!params.cloudProcess) { annotationGet( species_auto_ch ); annotation_auto_ch = annotationGet.out }
            // cloud storage file.exists()?
            if (params.cloudProcess) {
                annotation_preload = file("${params.permanentCacheDir}/annotations/${params.species}.gtf") // deprecated reminder
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
    take:
        setup_ch
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
    take:
        setup_ch
    main:
        if (!params.cloudProcess) { buscoGetDB(); database_busco = buscoGetDB.out }
        if (params.cloudProcess) { 
            busco_db_preload = file("${params.permanentCacheDir}/databases/busco/${params.busco_db}.tar.gz")
            if (busco_db_preload.exists()) { database_busco = busco_db_preload }
            else  { buscoGetDB(); database_busco = buscoGetDB.out }
        }
    emit: 
        database_busco
}

workflow download_dammit {
    take:
        setup_ch
    main:
    dammit_db_preload_path = "${params.permanentCacheDir}/databases/dammit/${params.busco_db}/dbs.tar.gz"
    if (params.dammit_uniref90) {
        dammit_db_preload_path = "${params.permanentCacheDir}/databases/dammit/uniref90/${params.busco_db}/dbs.tar.gz"
    }
    if (!params.cloudProcess) { dammitGetDB() ; database_dammit = dammitGetDB.out }
    if (params.cloudProcess) { 
        dammit_db_preload = file(dammit_db_preload_path)
        if (dammit_db_preload.exists()) { database_dammit = dammit_db_preload }
        else  { dammitGetDB(); database_dammit = dammitGetDB.out }
    }
    emit: database_dammit
        
}


/************************** 
* SUB WORKFLOWS
**************************/
/***************************************
Set up all databse and reference files for pipeline run without network connection:
annotation, genome reference, buscoDB, dammitDB, sortmernaDB
*/
workflow setup {
    //take:

    main:
        containerGet(container_ch)
        download_auto_annotation(container_ch)
        download_auto_reference(container_ch)
        download_busco(container_ch)
        download_dammit(container_ch)
        download_sortmerna(container_ch)
} 

/***************************************
Preprocess Illumina RNA-Seq reads: qc, trimming, adapters, rRNA-removal, mapping
*/
workflow preprocess_illumina {
    take:
        read_input_ch
        reference
        sortmerna_db

    main:
        // initial QC of raw reads
        fastqcPre(read_input_ch)

        if ( params.skip_read_preprocessing ) {
            // skip fastp, set dependent ch to empty
            smr_in = read_input_ch
            fastp_json_report = Channel.empty()
            readqcPost = Channel.empty()
        } else {
            // trim with fastp
            fastp(read_input_ch)
            smr_in = fastp.out.sample_trimmed
            fastp_json_report = fastp.out.json_report
            // QC after fastp
            fastqcPost(fastp.out.sample_trimmed)
            readqcPost = fastqcPost.out.zip
        }

        if ( params.skip_sortmerna && !params.skip_read_preprocessing ) {
            // skip SMR but not fastp
            sortmerna_no_rna_fastq = fastp.out.sample_trimmed
            sortmerna_log = Channel.empty()
        } else if ( params.skip_sortmerna && params.skip_read_preprocessing ) {
            // skip SMR and fastp
            sortmerna_no_rna_fastq = format_read_name(read_input_ch)
            sortmerna_log = Channel.empty()
        } else {
            // remove rRNA with SortmeRNA
            sortmerna(smr_in, extract_tar_bz2(sortmerna_db))
            sortmerna_no_rna_fastq = sortmerna.out.no_rna_fastq
            sortmerna_log = sortmerna.out.log
        }

        // HISAT2 index
        hisat2index(reference)
        // map with HISAT2
        hisat2(sortmerna_no_rna_fastq, hisat2index.out, params.hisat2_additional_params)
        // index BAM files
        index_bam(hisat2.out.sample_bam)

    emit:
        sample_bam_ch = hisat2.out.sample_bam
        fastp_json_report 
        sortmerna_log
        mapping_log = hisat2.out.log  
        readqcPre = fastqcPre.out.zip  
        readqcPost
        cleaned_reads_ch = sortmerna_no_rna_fastq
} 

/***************************************
Preprocess Illumina RNA-Seq reads: qc, trimming, adapters, rRNA-removal, mapping
*/
workflow preprocess_nanopore {
    take:
        read_input_ch
        reference
        sortmerna_db

    main:
        // initial QC of raw reads
        nanoplot(read_input_ch)

        if ( params.skip_sortmerna ) {
            sortmerna_no_rna_fastq = read_input_ch
            sortmerna_log = Channel.empty()
        } else {
            // remove rRNA with SortmeRNA
            sortmerna(read_input_ch, extract_tar_bz2(sortmerna_db))
            sortmerna_no_rna_fastq = sortmerna.out.no_rna_fastq
            sortmerna_log = sortmerna.out.log
        }

        // Minimap2 index
        minimap2index(reference)
        // map with Minimap2
        minimap2(sortmerna_no_rna_fastq, minimap2index.out, params.minimap2_additional_params)
        // index BAM files
        index_bam_minimap2(minimap2.out.sample_bam)

    emit:
        sample_bam_ch = minimap2.out.sample_bam
        fastp_json_report = Channel.empty()
        sortmerna_log
        mapping_log = minimap2.out.log  
        readqcPre = nanoplot.out.zip
        readqcPost = Channel.empty()
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
        mapping_log
        readqcPre
        readqcPost
        annotation
        deg_comparisons_input_ch
        deseq2_script
        deseq2_script_refactor_reportingtools_table
        deseq2_script_improve_deseq_table
        multiqc_config
        species2prefix

    main:
        // count with featurecounts
        featurecounts(sample_bam_ch, annotation, params.featurecounts_additional_params)

        // prepare annotation for R input
        format_annotation_gene_rows(annotation, gtf_feature_type_ch)
        format_annotation(annotation, gtf_attr_type_ch, gtf_feature_type_ch)

        // filter by TPM value
        // prepare input channels
        tpm_prep_ch = featurecounts.out.counts
            .join( annotated_reads
                    .map{meta, reads -> [meta.sample, meta.condition]} 
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
                .map{ meta, reads -> [ meta.sample, meta.condition, meta.source ] }
            )
            .multiMap{ it ->
                col_label: it[0]
                condition: it[1]
                source: it[2]
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
           annotated_sample.source.collect(), species_pathway_ch, deseq2_script, deseq2_id_type_ch, deseq2_script_refactor_reportingtools_table, 
           deseq2_script_improve_deseq_table, species2prefix)

        // run MultiQC
        multiqc_sample_names( annotated_reads.map{ meta, reads -> meta }.unique{ it.paired_end }, annotated_reads.map{ meta, reads -> [ meta.sample, reads ].flatten() }.collect() )
        multiqc(multiqc_config, 
                multiqc_sample_names.out,
                fastp_json_report.collect().ifEmpty([]), 
                sortmerna_log.collect().ifEmpty([]), 
                mapping_log.collect(), 
                featurecounts.out.log.collect(), 
                readqcPre.collect().ifEmpty([]),
                readqcPost.collect().ifEmpty([]),
                tpm_filter.out.stats,
                params.tpm,
                [],
                []
        )
} 

/*****************************************
De novo assembly of the preprocessed RNA-Seq reads. For now do co-assembly of all samples. 
ToDo: also do co-assembly of samples belonging to the same condition. 
*/
workflow assembly_denovo {
    take:
        cleaned_reads_ch
        dammit_db
        busco_db

    main:
        reads_ch = cleaned_reads_ch.map {meta, reads -> tuple reads}.collect()
        reads_input_csv = Channel.fromPath( params.reads, checkIfExists: true)
 

        // co-assembly LR
        if ( params.nanopore ){
            rattle(reads_ch)

            // qc check
            tool_ch = Channel.value('rattle')
            busco(rattle.out.assembly, busco_db, tool_ch)    

            // transcript annotation 
            dammit(rattle.out.assembly, busco_db, dammit_db, tool_ch)
        }else{
            // co-assembly SR
            trinity(reads_ch, reads_input_csv)

            // qc check
            tool_ch = Channel.value('trinity')
            busco(trinity.out.assembly, busco_db, tool_ch)    

            // transcript annotation 
            dammit(trinity.out.assembly, busco_db, dammit_db, tool_ch)
        }
} 

/*****************************************
Reference-based assembly of the preprocessed RNA-Seq reads. 
*/
workflow assembly_reference {
    take:
        genome_reference
        annotation_reference
        bams
        dammit_db
        busco_db

    main:
        // StringTie2 GTF-guided transcript prediction
        stringtie(genome_reference, annotation_reference, bams)

        // Merge each single GTF-guided StringTie2 GTF file 
        tool_ch = Channel.value('stringtie')
        stringtie_merge(genome_reference, stringtie.out.gtf.collect(), tool_ch)

        // qc check
        busco(stringtie_merge.out.transcripts, busco_db, tool_ch)    

        // transcript annotation 
        dammit(stringtie_merge.out.transcripts, busco_db, dammit_db, tool_ch)
}


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
    if (params.setup) {
        setup()
    } else {
        if ( workflow.profile.contains('test') ){
            get_test_data()
            reference = get_test_data.out.reference_test.collect()
            annotation = get_test_data.out.annotation_test.collect()
        } else {
            // get the reference genome
            download_auto_reference([])
            reference_auto = download_auto_reference.out

            // get the annotation
            download_auto_annotation([])
            annotation_auto = download_auto_annotation.out

            // concatenate genomes and annotations
            concat_genome(reference_custom_ch.collect().mix(reference_auto).collect())
            reference = concat_genome.out
            concat_annotation(annotation_custom_ch.collect().mix(annotation_auto).collect())
            annotation = concat_annotation.out
        }

        // get sortmerna databases
        if ( ! params.skip_sortmerna ) { 
            download_sortmerna([])
            sortmerna_db = download_sortmerna.out
        } else {
            sortmerna_db = Channel.empty()
        }
        
        // preprocess RNA-Seq reads (Illumina or Nanopore)
        if (!params.nanopore) { 
            preprocess_illumina(read_input_ch, reference, sortmerna_db) 
        } else { 
            preprocess_nanopore(read_input_ch, reference, sortmerna_db) 
        }

        // perform assembly & annotation
        if (params.assembly) {
            // dbs
            busco_db = download_busco([])
            dammit_db = download_dammit([])
            // de novo
            if (!params.nanopore) {
                // de novo
                assembly_denovo(preprocess_illumina.out.cleaned_reads_ch, dammit_db, busco_db)
                // reference-based
                assembly_reference(reference, annotation, preprocess_illumina.out.sample_bam_ch, dammit_db, busco_db)
            } else {
                // de novo
                assembly_denovo(preprocess_nanopore.out.cleaned_reads_ch, dammit_db, busco_db)
                // reference-based
                assembly_reference(reference, annotation, preprocess_nanopore.out.sample_bam_ch, dammit_db, busco_db)
            }
        } else {
        // perform expression analysis
            // start reference-based differential gene expression analysis
            if (!params.nanopore) { 
            expression_reference_based(preprocess_illumina.out.sample_bam_ch,
                                    preprocess_illumina.out.fastp_json_report,
                                    preprocess_illumina.out.sortmerna_log,
                                    preprocess_illumina.out.mapping_log,
                                    preprocess_illumina.out.readqcPre,
                                    preprocess_illumina.out.readqcPost,
                                    annotation,
                                    deg_comparisons_input_ch, 
                                    deseq2_script, 
                                    deseq2_script_refactor_reportingtools_table, 
                                    deseq2_script_improve_deseq_table, 
                                    multiqc_config,
                                    species2prefix)
            } else {
            expression_reference_based(preprocess_nanopore.out.sample_bam_ch,
                                    preprocess_nanopore.out.fastp_json_report,
                                    preprocess_nanopore.out.sortmerna_log,
                                    preprocess_nanopore.out.mapping_log,
                                    preprocess_nanopore.out.readqcPre,
                                    preprocess_nanopore.out.readqcPost,
                                    annotation,
                                    deg_comparisons_input_ch, 
                                    deseq2_script, 
                                    deseq2_script_refactor_reportingtools_table, 
                                    deseq2_script_improve_deseq_table, 
                                    multiqc_config,
                                    species2prefix)
            }
        }
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

    ${c_yellow}Usage examples:${c_reset}
    nextflow run hoelzer-lab/rnaflow -profile test,local,conda
    nextflow run hoelzer-lab/rnaflow --cores 4 --reads input.csv --autodownload mmu --pathway mmu
    nextflow run hoelzer-lab/rnaflow --cores 4 --reads input.csv --autodownload eco --assembly
    nextflow run hoelzer-lab/rnaflow --cores 4 --reads input.csv --genome fasta_virus.csv --annotation gtf_virus.csv --autodownload hsa --pathway hsa
    ${c_dim}Genomes and annotations from --autodownload, --genome and --annotation are concatenated.${c_reset}

    ${c_yellow}Input:${c_reset}
    ${c_green}--reads${c_reset}                  A CSV file following the pattern: Sample,R,Condition,Source for single-end or Sample,R1,R2,Condition,Source for paired-end
                                        ${c_dim}(check terminal output if correctly assigned)
                                        Per default, all possible comparisons of conditions in one direction are made. Use --deg to change.${c_reset}
    ${c_green}--autodownload${c_reset}           Specifies the species identifier for automated download [default: $params.autodownload]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]
                                        - mau [Ensembl: Mesocricetus_auratus.MesAur1.0.dna.toplevel | Mesocricetus_auratus.MesAur1.0.100]${c_reset}
    ${c_dim}--species                Specifies the species identifier for downstream path analysis. (DEPRECATED)
                             If `--include_species` is set, reference genome and annotation are added and automatically downloaded. [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]
                                        - mau [Ensembl: Mesocricetus_auratus.MesAur1.0.dna.toplevel | Mesocricetus_auratus.MesAur1.0.100]${c_reset}
    ${c_green}--genome${c_reset}                 CSV file with genome reference FASTA files (one path in each line)
                                        ${c_dim}If set, --annotation must also be set.${c_reset}
    ${c_green}--annotation${c_reset}             CSV file with genome annotation GTF files (one path in each line)
    ${c_dim}--include_species        Either --species or --genome/--annotation need to be used. Both input seetings can be also combined to use genome and annotation of 
                             supported species in addition to --genome and --annotation (DEPRECATED) [default: $params.include_species]${c_reset}

    ${c_yellow}Preprocessing options:${c_reset}
    --mode                             Either 'single' (single-end) or 'paired' (paired-end) sequencing [default: $params.mode]
    --fastp_additional_params          Additional parameters for fastp [default: $params.fastp_additional_params]
    --skip_sortmerna                   Skip rRNA removal via SortMeRNA [default: $params.skip_sortmerna] 
    --skip_read_preprocessing          Skip preprocessing with fastp [default: $params.skip_read_preprocessing]
    --hisat2_additional_params         Additional parameters for HISAT2 [default: $params.hisat2_additional_params]
    --minimap2_additional_params       Additional parameters for minimap2 (Nanopore input) [default: $params.minimap2_additional_params]
    --featurecounts_additional_params  Additional parameters for FeatureCounts [default: $params.featurecounts_additional_params]
    --nanopore                         Activate Nanopore long-read mode (default is Illumina data) [default: $params.nanopore]  

    ${c_yellow}DEG analysis options:${c_reset}
    --strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default: $params.strand]
    --tpm                    Threshold for TPM (transcripts per million) filter. A feature is discared, if for all conditions the mean TPM value of all 
                             corresponding samples in this condition is below the threshold. [default: $params.tpm]
    --deg                    A CSV file following the pattern: conditionX,conditionY
                             Each line stands for one differential gene expression comparison.  
                             Must match the 'Condition' labels defined in the CSV file provided via --reads.  
    --pathway                Perform different downstream pathway analysis for the species. [default: $params.pathway]
                             ${c_dim}Currently supported are:
                                 - hsa | Homo sapiens
                                 - mmu | Mus musculus
                                 - mau | Mesocricetus auratus${c_reset}
    --feature_id_type        ID type for downstream analysis [default: $params.feature_id_type]

    ${c_yellow}Transcriptome assembly options:${c_reset}
    --assembly               Perform de novo and reference-based transcriptome assembly instead of DEG analysis [default: $params.assembly]
    --busco_db               The database used with BUSCO [default: $params.busco_db]
                             ${c_dim}Full list of available data sets at https://busco-data.ezlab.org/v5/data/lineages/ ${c_reset}
    --dammit_uniref90        Add UniRef90 to the dammit databases (time consuming!) [default: $params.dammit_uniref90]
    --rna                    Activate directRNA mode for ONT transcriptome assembly [default: $params.rna (cDNA)]


    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]
    --output                 Name of the result folder [default: $params.output]

    ${c_yellow}Caching:${c_reset}
    --permanentCacheDir      Location for auto-download data like databases [default: $params.permanentCacheDir]
    --condaCacheDir          Location for storing the conda environments [default: $params.condaCacheDir]
    --singularityCacheDir    Location for storing the singularity images [default: $params.singularityCacheDir]
    ${c_dim}--workdir                Working directory for all intermediate results [default: $params.workdir] (DEPRECATED: use `-w your/workdir` instead)${c_reset}
    --softlink_results       Softlink result files instead of copying.
    --setup                  Download all necessary DB, reference and image files without running the pipeline. [default: $params.setup]


    ${c_dim}Nextflow options:
    -with-tower              Activate monitoring via Nextflow Tower (needs TOWER_ACCESS_TOKEN set).
    -with-report rep.html    CPU / RAM usage (may cause errors).
    -with-dag chart.html     Generates a flowchart for the process tree.
    -with-timeline time.html Timeline (may cause errors).${c_reset}

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
      lsf
      latency

    
    ${c_blue}Engines${c_reset} (choose one):
      conda
      mamba
      docker
      singularity
    
    Per default: -profile local,conda is executed. 

    ${c_dim}For a test run (~ 15 min), add "test" to the profile, e.g. -profile test,local,conda.
    The command will create all conda environments and download and run test data.

    We also provide some pre-configured profiles for certain HPC environments:    
      ara (slurm, conda and parameter customization)
    ${c_reset}
    """.stripIndent()
}
