process stringtie {
    label 'stringtie'  

    if ( params.softlink_results ) {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/stringtie", pattern: "${sample_name}_stringtie.gtf"
      //publishDir "${params.output}/${name}/${params.rnaseq_annotation_dir}/stringtie", pattern: "${sample_name}_stringtie.fna"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/stringtie", pattern: "${sample_name}_gene_abundance.txt"
    } else {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/stringtie", mode: 'copy', pattern: "${sample_name}_stringtie.gtf"
      //publishDir "${params.output}/${name}/${params.rnaseq_annotation_dir}/stringtie", mode: 'copy', pattern: "${sample_name}_stringtie.fna"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/stringtie", mode: 'copy', pattern: "${sample_name}_gene_abundance.txt"
    }

  input:
    path genome
    path gtf
    tuple val(sample_name), file(bam)

  output:
    path "${sample_name}_stringtie.gtf", emit: gtf
    path "${sample_name}_gene_abundance.txt", emit: abundance
    //tuple val(name), file("${sample_name}_stringtie.fna"), emit: transcripts

  script:
    """
      stringtie -G ${gtf} -p ${task.cpus} -o ${sample_name}_stringtie.gtf -A ${sample_name}_gene_abundance.txt ${bam} 
      #gffread -w ${sample_name}_stringtie.fna -g ${genome} ${sample_name}_stringtie.gtf
    """
  }

process stringtie_merge {
    label 'stringtie'  

    errorStrategy { task.exitStatus = 1 ? 'ignore' :  'terminate' }

    if ( params.softlink_results ) { 
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/${mode}", pattern: "${mode}.gtf"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/${mode}", pattern: "${mode}.fna"
    } else {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/${mode}", mode: 'copy', pattern: "${mode}.gtf"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/${mode}", mode: 'copy', pattern: "${mode}.fna"
    }

  input:
    path(genome)
    file(stringtie_gtfs)
    val(mode)

  output:
    path "${mode}.gtf", emit: gtf
    path "${mode}.fna", emit: transcripts

  script:
    """
      stringtie --merge -o ${mode}.gtf *.gtf
      gffread -w ${mode}.fna -g ${genome} ${mode}.gtf
    """
  }

/*
-G reference_guide.gff (gtf/gff3)

THIS COULD BE INTERESTING TO USE INFORMATION FROM MULTIPLE RNA-SEQ FILES
Transcript merge usage mode: 
  stringtie --merge [Options] { gtf_list | strg1.gtf ...}
With this option StringTie will assemble transcripts from multiple
input files generating a unified non-redundant set of isoforms. In this mode
the following options are available:
  -G <guide_gff>   reference annotation to include in the merging (GTF/GFF3)
  -o <out_gtf>     output file name for the merged transcripts GTF
                    (default: stdout)
  -m <min_len>     minimum input transcript length to include in the merge
                    (default: 50)
  -c <min_cov>     minimum input transcript coverage to include in the merge
                    (default: 0)
  -F <min_fpkm>    minimum input transcript FPKM to include in the merge
                    (default: 1.0)
  -T <min_tpm>     minimum input transcript TPM to include in the merge
                    (default: 1.0)
  -f <min_iso>     minimum isoform fraction (default: 0.01)
  -g <gap_len>     gap between transcripts to merge together (default: 250)
  -i               keep merged transcripts with retained introns; by default
                   these are not kept unless there is strong evidence for them
  -l <label>       name prefix for output transcripts (default: MSTRG)
*/