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
    tuple val(meta), path("*metrics.json"), emit: metrics
    tuple val(meta), path("*raw.h5seurat"), emit: h5seurat

    script:
    """
    Rscript -e 'rmarkdown::render("$baseDir/bin/scRNA_report.Rmd", params = list(sampleName = "${meta.id}", starsolo_out = "$starsolo_outdir/Summary.csv", qualimap_out = "$qualimap_outdir/rnaseq_qc_results.txt", qualimap_gene_coverage = "$qualimap_outdir/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt", starsolo_bc = "$starsolo_outdir/UMIperCellSorted.txt", starsolo_matrixDir="$starsolo_outdir/filtered", nCPUs = "$task.cpus", saturation_json = "${saturation_outJSON}"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${meta.id}_report.html")'

    """
}