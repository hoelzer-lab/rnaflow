process trinity {
    label 'trinity'  
    publishDir "${params.output}/${params.assembly_dir}/", mode: 'copy', pattern: "trinity.fasta"

  input:
    path reads 
    path csv

  output:
    path "trinity.fasta", emit: assembly

  script:
    if (params.mode == 'paired')
    """
      # Update the original CSV file to match quality controlled reads and Trinity input
      grep -v Condition ${csv} | awk 'BEGIN{FS=","}{print \$4"\\t"\$1"\\t"\$1".R1.other.fastq.gz\\t"\$1".R2.other.fastq.gz"}' > \$(basename \$PWD)_input.csv
      MEM=\$(echo ${task.memory} | awk '{print \$1}')
      Trinity --seqType fq --samples_file \$(basename \$PWD)_input.csv --max_memory \${MEM}G --CPU ${task.cpus}
      mv trinity_out_dir/Trinity.fasta trinity.fasta
    """
    else if (params.mode == 'single')
    """
      # Update the original CSV file to match quality controlled reads and Trinity input
      grep -v Condition ${csv} | awk 'BEGIN{FS=","}{print \$3"\\t"\$1"\\t"\$1".other.fastq.gz"}' > \$(basename \$PWD)_input.csv
      MEM=\$(echo ${task.memory} | awk '{print \$1}')
      Trinity --seqType fq --samples_file \$(basename \$PWD)_input.csv --max_memory \${MEM}G --CPU ${task.cpus}
      mv trinity_out_dir/Trinity.fasta trinity.fasta
    """
    else 
      error "Invalid read mode: ${params.mode}"
  }





/*
#      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                      # if single-end instead of paired-end, then leave the 4th column above empty.
*/