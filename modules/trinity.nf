process trinity {
    label 'trinity'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.assembly_dir}/Trinity", pattern: "trinity.fasta" }
    else { publishDir "${params.output}/${params.assembly_dir}/Trinity", mode: 'copy', pattern: "trinity.fasta" }

  input:
    val (meta)
    path (reads)
    path (csv)

  output:
    path "trinity.fasta", emit: assembly

  script:
    if (meta.paired)
    """
      # check if sortmerna was used and adjust file names accordingly
      TYPE='other'
      if (${params.skip_sortmerna} == 'true'); then
        TYPE='trimmed'
        if (${params.skip_read_preprocessing} == 'true'); then
          TYPE='raw';
        fi
      fi

      # Update the original CSV file to match quality controlled reads and Trinity input
      for SAMPLE in \$(grep -v Sample ${csv} | awk 'BEGIN{FS=","};{print \$1}');
        do CONDITION=\$(grep \$SAMPLE ${csv} | awk 'BEGIN{FS=","};{print \$4}'); 
        printf \$CONDITION"\\t"\$SAMPLE"\\t"\$SAMPLE".\$TYPE.fastq.gz\\n"; 
        done > \$(basename \$PWD)_input.csv
      
      MEM=\$(echo ${task.memory} | awk '{print \$1}')
      Trinity --seqType fq --samples_file \$(basename \$PWD)_input.csv --max_memory \${MEM}G --bflyCalculateCPU --CPU ${task.cpus} --output trinity_out_dir/
      mv trinity_out_dir*.fasta trinity.fasta
    """
    else
    """
      # check if sortmerna was used and adjust file names accordingly
      TYPE='other'
      if (${params.skip_sortmerna} == 'true'); then
        TYPE='trimmed'
        if (${params.skip_read_preprocessing} == 'true'); then
          TYPE='raw';
        fi
      fi

      # Update the original CSV file to match quality controlled reads and Trinity input
      for SAMPLE in \$(grep -v Sample ${csv} | awk 'BEGIN{FS=","};{print \$1}'); 
        do CONDITION=\$(grep \$SAMPLE ${csv} | awk 'BEGIN{FS=","};{print \$4}'); 
        printf \$CONDITION"\\t"\$SAMPLE"\\t"\$SAMPLE".\$TYPE.fastq.gz\\n"; 
        done > \$(basename \$PWD)_input.csv

      MEM=\$(echo ${task.memory} | awk '{print \$1}')
      Trinity --seqType fq --samples_file \$(basename \$PWD)_input.csv --max_memory \${MEM}G --bflyCalculateCPU --CPU ${task.cpus}  --output trinity_out_dir/
      mv trinity_out_dir*.fasta trinity.fasta
    """
  }

//In addition, the --bflyHeapSpaceMax is available. If you are confident that no instances of Butterfly will use all 10GB of memory, setting this to a smaller value may allow more Butterfly processes to run.

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
