process busco {
    label 'busco'

    if ( params.softlink_results ) {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", pattern: "busco_*_summary.txt"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", pattern: "busco_*_figure.pdf"
    } else {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", mode: 'copy', pattern: "busco_*_summary.txt"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", mode: 'copy', pattern: "busco_*_figure.pdf"
    }
    
    { publishDir "${params.permanentCacheDir}/databases/busco", mode: 'copy', pattern: "${params.busco_db}_odb10.tar.gz" }


    input:
      path fasta
      val preload
      val tool
      

    output:
      path("busco_${tool}_summary.txt"), emit: summary
      path("busco_${tool}_figure.pdf")
      path "full_table_${tool}_results.tsv"
      path "${params.busco_db}_odb10.tar.gz", emit: odb10
    
    script:
      if ( preload )
        """
        # run busco
        busco -i ${fasta} -o results  -m tran -c ${task.cpus} -l ${params.permanentCacheDir}/databases/busco/${params.busco_db}_odb10.tar.gz 
        cp results/run_${params.busco_db}_odb10/short_summary.txt busco_${tool}_summary.txt

        # generate Plot and rehack Rscript
        generate_plot.py -wd results/
        sed -i 's/busco_figure.png/busco_figure.pdf/g' results/busco_figure.R
        Rscript results/busco_figure.R
        cp results/busco_figure.pdf busco_${tool}_figure.pdf
        cp results/run_${params.busco_db}_odb10/full_table.tsv full_table_${tool}_results.tsv

        """
      else
        """
        # run busco
        busco -i ${fasta} -o results  -m tran -c ${task.cpus} -l ${params.busco_db} #--download_path ${params.busco_db}.tar.gz busco_downloads/lineages/..
        cp results/run_${params.busco_db}_odb10/short_summary.txt busco_${tool}_summary.txt

        # generate Plot and rehack Rscript
        generate_plot.py -wd results/
        sed -i 's/busco_figure.png/busco_figure.pdf/g' results/busco_figure.R
        Rscript results/busco_figure.R
        cp results/busco_figure.pdf busco_${tool}_figure.pdf
        cp results/run_${params.busco_db}_odb10/full_table.tsv full_table_${tool}_results.tsv

        mv busco_downloads/lineages/${params.busco_db}_odb10 .
        tar -czf ${params.busco_db}_odb10.tar.gz ${params.busco_db}_odb10 
        """
}

/* Comments:
Buscos plotting feature does not work in Docker by default.
The "hack" modifies via sed the Rscript that gets generated after generate_plot.py
After changing to pdf (circumventing resolutions) we rerun the script via Rscript. Tadaa... we have a plot
*/
