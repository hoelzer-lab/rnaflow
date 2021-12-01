process rattle {
    label 'rattle'
    time '48h'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.assembly_dir}/rattle", pattern: "transcriptome.fq" }
    else { publishDir "${params.output}/${params.assembly_dir}/rattle", mode: 'copy', pattern: "transcriptome.fq" }

  input:
    path (reads) 

  output:
    path "transcriptome.fq", emit: assembly

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
    rm *fastq
    """
  }