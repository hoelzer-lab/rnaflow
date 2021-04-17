/************************************************************************
* Nanoplot
************************************************************************/
process nanoplot {
    label 'nanoplot'
    publishDir "${params.output}/${params.readqcdir}/${name}/", mode: 'copy', pattern: "${name}_read_quality_report.html"
    publishDir "${params.output}/${params.readqcdir}/${name}/", mode: 'copy', pattern: "${name}_read_quality.txt"
    publishDir "${params.output}/${params.readqcdir}/${name}/figures", mode: 'copy', pattern: "*.png"
    publishDir "${params.output}/${params.readqcdir}/${name}/vector_figures", mode: 'copy', pattern: "*.pdf"
    // addet ignore here - to avoid pipeline breaking
    errorStrategy 'ignore'
    input:
      tuple val(name), path(reads)
    output:
      path("${name}_nanoplot.zip", emit: zip) optional true
      tuple val(name), path("*.html"), path("*.pdf") optional true
      tuple val(name), path("${name}_read_quality.txt"), path("*.png") optional true
    script:
      """
      NanoPlot -t ${task.cpus} --fastq ${reads} --title '${name}' --color darkslategrey --N50 --plots hex --loglength -f png --store
      NanoPlot -t ${task.cpus} --pickle NanoPlot-data.pickle --title '${name}' --color darkslategrey --N50 --plots hex --loglength -f pdf
      mv NanoPlot-report.html ${name}_read_quality_report.html
      mv NanoStats.txt ${name}_read_quality.txt
      zip ${name}_nanoplot.zip *.html *.pdf *.png ${name}_read_quality.txt
      """
}

/* Comments:
We run nanoplot 2 times to get png and pdf files.
The second time its done via the pickle file of the previous run, to save computing time
*/


