#!/usr/bin/env nextflow

/*
* RNA-Seq-based detection of differentially expressed genes
*
* Author: martin.hoelzer@uni-jena.de
* Author: marie.lataretu@uni-jena.de
*/

// nextflow run main.nf --threads 8 --reads input.se.csv --reference data/db/Rattus_norvegicus.Rnor_6.0.dna.toplevel.chr.fa --annotation data/db/Rattus_norvegicus.Rnor_6.0.91.chr.gtf --mode single --index data/db/index

def helpMSG() {
    log.info """

    Usage:
    nextflow run main.nf --threads 4 --reads input.csv --reference data/db/Rattus_norvegicus.Rnor_6.0.dna.toplevel.chr.fa --annotation data/db/Rattus_norvegicus.Rnor_6.0.91.chr.gtf

    Mandatory:
    --reads         a CSV file following the pattern: mock_rep1,fastq1,fastq2 with fastq2 beeing optional 
                    (check terminal output if correctly assigned)
    --reference     a reference genome with HISAT2 index for mapping
    --annotation    a annotation file in gtf format corresponding to the reference file

    Options
    --sortmerna              the database used for SortMeRNA
    --index                  the path to the hisat2 index prefix matching the genome provided via --reference. 
                             If provided, no new index will be build. Must be named 'index.*.ht2'. 
                             Simply provide the path like 'data/db/index'
    --mode                   either 'single' (single-end) or 'paired' (paired-end) sequencing [default $params.mode]
    --strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default $params.strand]
    --threads                max cores for local use [default $params.threads]
    --mem                    max memory in GB for local use [default $params.mem]
    --output                 name of the result folder [default $params.output]

    Profile:
    -profile                standard, gcloud (wip) [default: standard]

    """.stripIndent()
}

if (params.help) { exit 0, helpMSG() }
if (params.reads == '') {exit 1, "--reads is a required parameter"}
if (params.reference == '') {exit 1, "--reference is a required parameter"}
if (params.annotation == '') {exit 1, "--annotation is a required parameter"}

println "\nD I F F E R E N T I A L  G E N E  E X P R E S S I O N  A N A L Y S I S"
println "= = = = = = = = = = = =  = = = =  = = = = = = = = = =  = = = = = = = ="
println "Reference genome:     $params.reference"
println "Reference annotation: $params.annotation"
println "SortMeRNA database:   $params.sortmerna_db"
println "Output path:          $params.output\n"
println "mode:                 $params.mode\n"

// file channel from CSV
if (params.reads) { 
readFileList_ch = Channel
        .fromPath( params.reads , checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", [file("${row[1]}"), file("${row[2]}")]] }
        .into { reads_ch; reads_report_ch}
        reads_report_ch.subscribe { println "Got short reads: ${it}" }
}

reference_file = file(params.reference)

if (params.index) {
  index_ch = Channel.fromPath("${params.index}.*", checkIfExists: true)
}

/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process trimming {
  conda 'envs/fastp.yaml'
  publishDir "${params.output}/${params.trimming_dir}", mode: 'copy', pattern: "${name}*.trimmed.fastq"

  input:
  set val(name), file(reads) from reads_ch

  output:
  set val(name), file("${name}*.trimmed.fastq") into trimming_ch

  shell:
  if (params.mode == 'single') {
  """
  fastp -i !{reads[0]} -o !{name}.trimmed.fastq -n 5 --thread !{params.threads}  
  """
  }
  else {
  """
  fastp -i !{reads[0]} -I !{reads[1]} -o !{name}.R1.trimmed.fastq -O !{name}.R2.trimmed.fastq -n 5 --thread !{params.threads}  
  """         
  }
}

/************************************************************************
* INDEX & MAPPING
************************************************************************/
if (!params.index) {
  process index {
    conda 'envs/hisat2.yaml'

    input:
    file reference_file

    output:
    file("index.*") into index_ch

    shell:
    """
    hisat2-build -p !{params.threads} !{reference_file} index
    """
  }
}

process mapping {
  conda 'envs/hisat2.yaml'
  publishDir "${params.output}/${params.mapping_dir}", mode: 'copy', pattern: "${name}.sorted.bam"

  input:
  set val(name), file(reads) from trimming_ch
  file(database) from index_ch.collect()

  output:
  set val(name), file("${name}.sorted.bam") into mapping_ch

  shell:
  if (params.mode == 'single') {
  """
  hisat2 -x index -U !{reads[0]} -p !{params.threads} | samtools view -bS | samtools sort -o !{name}.sorted.bam -T tmp --threads !{params.threads}
  """
  }
  else {
  """
  hisat2 -x index -1 !{reads[0]} -2 !{reads[1]} -p !{params.threads} | samtools view -bS | samtools sort -o !{name}.sorted.bam -T tmp --threads !{params.threads}
  """
  } 
}

/************************************************************************
* COUNTING 
* TODO: check that we do not miss a gene due to the 'sed 1d' remove of the first two lines
************************************************************************/
process counting {
  conda 'envs/subread.yaml'
  publishDir "${params.output}/${params.counting_dir}", mode: 'copy', pattern: "${name}.counts*"

  input:
  set val(name), file(bam) from mapping_ch
  file(annotation) from file(params.annotation)

  output:
  set val(name), file("${name}.counts*") into counting_ch // [mock_rep1, [/home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.counts, /home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.counts.formated, /home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.counts.summary]]

  shell:
  if (params.mode == 'single') {
  '''
  featureCounts -T !{params.threads} -s !{params.strand} -a !{annotation} -o !{name}.counts !{bam}
  awk '{print $1"\t"$7}' !{name}.counts | sed 1d | sed 1d > !{name}.counts.formated
  '''
  }
  else {
  '''
  featureCounts -pBP -T !{params.threads} -s !{params.strand} -a !{annotation} -o !{name}.counts !{bam}
  awk '{print $1"\t"$7}' !{name}.counts | sed 1d | sed 1d > !{name}.counts.formated
  '''
  }
}

/************************************************************************
* DIFFERENTIAL EXPRESSION 
************************************************************************/

// maybe translate this ruby script into python
process prepare_annotation {
input: 
file(annotation) from file(params.annotation)

output:
file("${annotation}.id2ensembl") into prepare_annotation_ch

shell:
'''
#!/usr/bin/env python3
with open("!{annotation}", 'r') as gtf, open("!{annotation}.id2ensembl", 'a') as out:
  for line in gtf:
    if not line.startswith('#'):
      split_line = line.split('\\t')
      if split_line[2] == 'gene':
        desc = split_line[8]
        gene_id = line.split('gene_id')[1].split(';')[0].replace('"', '').strip()
        if 'gene_name' in line:
          gene_name = desc.split('gene_name')[1].split(';')[0].replace('"','').strip()
        else:
          gene_name = gene_id
        gene_biotype = desc.split('gene_biotype')[1].split(';')[0].replace('"','').strip()
        out.write('\\t'.join([gene_id, gene_name, gene_biotype]) + '\\n')
'''
}

process diff {
  publishDir "${params.output}/${params.diff_dir}", mode: 'copy', pattern: "*"

  input:
  set val(name), file(counts) from counting_ch
  file(annotation) from file(params.annotation)
  file(ensembl2id) from prepare_annotation_ch

  output:
  set val(name), file(counts) into final_ch 

  shell:
  '''
  echo !{name}.counts.formated
  '''
}

/*
project_dir <- "/home/hoelzer/git/nanozoo/wf_gene_expression/results/" 
samples <- c("/home/hoelzer/git/nanozoo/wf_gene_expression/results/03-Counting/mock_rep1.counts.formated","/home/hoelzer/git/nanozoo/wf_gene_expression/results/03-Counting/mock_rep2.counts.formated","/home/hoelzer/git/nanozoo/wf_gene_expression/results/03-Counting/mock_rep3.counts.formated","/home/hoelzer/git/nanozoo/wf_gene_expression/results/03-Counting/treated_rep1.counts.formated","/home/hoelzer/git/nanozoo/wf_gene_expression/results/03-Counting/treated_rep2.counts.formated","/home/hoelzer/git/nanozoo/wf_gene_expression/results/03-Counting/treated_rep3.counts.formated") 
conditions <- c("mock","mock","mock","treated","treated","treated") 
col.labels <- c("mock_rep1","mock_rep2","mock_rep3","treated_rep1","treated_rep2","treated_rep3") 
levels <- c("mock","treated") 
comparisons <- c("mock:treated")
ensembl2genes <- "data/db/Rattus_norvegicus.Rnor_6.0.91.chr.id2name"
species <- "rno"
patients <- c()
*/