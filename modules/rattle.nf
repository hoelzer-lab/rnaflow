process rattle {
    label 'rattle'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.assembly_dir}/rattle", pattern: "transcriptome.fa" }
    else { publishDir "${params.output}/${params.assembly_dir}/rattle", mode: 'copy', pattern: "transcriptome.fa" }

  input:
    path (reads) 

  output:
    path "transcriptome.fa", emit: assembly

  script:
    """
    # unzip and filter for sequences shorter than k-mer size
    zcat -c ${reads} | awk 'BEGIN{OFS="\\n"} /^@/ {getline s; getline a; getline q; if (length(s) >= 100) {print \$0,s,a,q}}' > filtered.fastq

    # clustering step
    mkdir -p output/clusters/

    if ( ${params.rna} == 'true' ); then
      rattle cluster -i filtered.fastq -t ${task.cpus} -o ./output/ --fastq --rna 
    else
      rattle cluster -i filtered.fastq -t ${task.cpus} -o ./output/ --fastq 
    fi

    # extract clusters
    rattle extract_clusters -i filtered.fastq -c ./output/clusters.out -o ./output/clusters --fastq 

    # error correction
    rattle correct -i filtered.fastq -c ./output/clusters.out -o ./output/ -t ${task.cpus} --fastq

    # polish
    if ( ${params.rna} == 'true' ); then
      rattle polish -i ./output/consensi.fq -o ./output/  -t ${task.cpus} --rna
    else
      rattle polish -i ./output/consensi.fq -o ./output/  -t ${task.cpus}
    fi

    mv output/transcriptome.fq transcriptome.fq

    # convert to fasta for busco etc.
    cat transcriptome.fq | awk '{if(NR%4==1) {printf(">%s\\n",substr(\$0,2));} else if(NR%4==2) print;}' > transcriptome.fa
    rm *fastq
    """
  }
