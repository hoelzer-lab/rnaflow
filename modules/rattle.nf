process rattle {
    label 'rattle'
    time '48h'

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
    rattle cluster -i filtered.fastq -t ${task.cpus} -o ./output/ --fastq --rna #(direct rna seq, disables checking both strands)

    # extract clusters
    rattle extract_clusters -i filtered.fastq -c ./output/clusters.out -o ./output/clusters --fastq 

    # error correction
    rattle correct -i filtered.fastq -c ./output/clusters.out -o ./output/ -t ${task.cpus} --fastq

    # polish
    rattle polish -i ./output/consensi.fq -o ./output/  -t ${task.cpus} --rna

    mv output/transcriptome.fq transcriptome.fq

    # convert to fasta for busco etc.
    cat transcriptome.fq | awk '{if(NR%4==1) {printf(">%s\\n",substr(\$0,2));} else if(NR%4==2) print;}' > transcriptome.fa
    rm *fastq
    """
  }