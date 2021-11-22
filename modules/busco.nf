process busco {
    label 'busco'

    if ( params.softlink_results ) {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", pattern: "busco_*_summary.txt"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", pattern: "busco_*_figure.pdf"
    } else {
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", mode: 'copy', pattern: "busco_*_summary.txt"
      publishDir "${params.output}/${params.rnaseq_annotation_dir}/BUSCO", mode: 'copy', pattern: "busco_*_figure.pdf"
    }

    busco_db_preload = file("${params.permanentCacheDir}/databases/busco/${params.busco_db}/${params.busco_db}*")

    if ( busco_db_preload.exists() ) { storeDir "${params.permanentCacheDir}/databases/busco/${params.busco_db}" }


    input:
      path fasta
      //path database
      val tool
      
    output:
      tuple path("busco_${tool}_summary.txt"), path("busco_${tool}_figure.pdf")
      path "full_table_${tool}_results.tsv"
    
    script:

      if ( busco_db_preload.exists() ) {
        """
        tar -zxf ${database}

        # run busco
        busco -i ${fasta} -o results  -m tran -c ${task.cpus} -l ${busco_db_preload} 
        cp run_results/short_summary_results.txt busco_${tool}_summary.txt

        # generate Plot and rehack Rscript
        generate_plot.py -wd run_results/
        sed -i 's/busco_figure.png/busco_figure.pdf/g' run_results/busco_figure.R
        Rscript run_results/busco_figure.R
        cp run_results/busco_figure.pdf busco_${tool}_figure.pdf
        cp run_results/full_table_results.tsv full_table_${tool}_results.tsv
        """ 
      } else  {
        """
        tar -zxf ${database}

        # run busco
        busco -i ${fasta} -o results  -m tran -c ${task.cpus} -l ${params.busco_db} #--download-path
        cp run_results/short_summary_results.txt busco_${tool}_summary.txt

        # generate Plot and rehack Rscript
        generate_plot.py -wd run_results/
        sed -i 's/busco_figure.png/busco_figure.pdf/g' run_results/busco_figure.R
        Rscript run_results/busco_figure.R
        cp run_results/busco_figure.pdf busco_${tool}_figure.pdf
        cp run_results/full_table_results.tsv full_table_${tool}_results.tsv
        """
      }
}

/* Comments:
Buscos plotting feature does not work in Docker by default.
The "hack" modifies via sed the Rscript that gets generated after generate_plot.py
After changing to pdf (circumventing resolutions) we rerun the script via Rscript. Tadaa... we have a plot
*/
