#! /usr/bin/env nextflow

process REPORT{
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(starsolo_outdir)
    tuple val(meta), path(qualimap_outdir)
    tuple val(meta), path(saturation_outJSON)

    output:
    tuple val(meta), path("*report.html"), emit: report

    script:
    """
    Rscript -e 'rmarkdown::render("$baseDir/bin/scRNA_report.Rmd", params = list(sampleName = "${meta.id}", starsolo_out = "$starsolo_outdir/Summary.csv", qualimap_out = "$qualimap_outdir/rnaseq_qc_results.txt", starsolo_bc = "$starsolo_outdir/UMIperCellSorted.txt", starsolo_matrixDir="$starsolo_outdir/filtered", nCPUs = "$task.cpus", saturation_json = "${saturation_outJSON}"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${meta.id}_report.html")'

    """
}