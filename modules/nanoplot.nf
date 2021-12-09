/************************************************************************
* Nanoplot
************************************************************************/
process nanoplot {
    label 'nanoplot'
    tag "$meta.sample"
    publishDir "${params.output}/${params.readqc_dir}/${meta.sample}/", mode: 'copy', pattern: "${meta.sample}_read_quality_report.html"
    publishDir "${params.output}/${params.readqc_dir}/${meta.sample}/", mode: 'copy', pattern: "${meta.sample}_read_quality.txt"
    publishDir "${params.output}/${params.readqc_dir}/${meta.sample}/figures", mode: 'copy', pattern: "*.png"
    publishDir "${params.output}/${params.readqc_dir}/${meta.sample}/vector_figures", mode: 'copy', pattern: "*.pdf"
    // addet ignore here - to avoid pipeline breaking
    errorStrategy 'ignore'
    input:
      tuple val(name), path(reads)
    output:
      path("${meta.sample}_nanoplot.zip", emit: zip) optional true
      tuple val(name), path("*.html"), path("*.pdf") optional true
      tuple val(name), path("${meta.sample}_read_quality.txt"), path("*.png") optional true
    script:
      """
      NanoPlot -t ${task.cpus} --fastq ${reads} --title '${meta.sample}' --color darkslategrey --N50 --plots hex --loglength -f png --store
      NanoPlot -t ${task.cpus} --pickle NanoPlot-data.pickle --title '${meta.sample}' --color darkslategrey --N50 --plots hex --loglength -f pdf
      mv NanoPlot-report.html ${meta.sample}_read_quality_report.html
      mv NanoStats.txt ${meta.sample}_read_quality.txt
      zip ${meta.sample}_nanoplot.zip *.html *.pdf *.png ${meta.sample}_read_quality.txt
      """
}

/* Comments:
We run nanoplot 2 times to get png and pdf files.
The second time its done via the pickle file of the previous run, to save computing time
*/


