#! /usr/bin/env nextflow

process REPORT{
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(starsolo_summary)
    tuple val(meta), path(starsolo_UMI_file)
    tuple val(meta), path(starsolo_filteredDir)
    tuple val(meta), path(qualimap_outdir)
    tuple val(meta), path(saturation_outJSON)
    tuple val(meta), path(version_json)

    output:
    tuple val(meta), path("*report.html") , emit: report
    tuple val(meta), path("*metrics.json"), emit: metrics
    tuple val(meta), path("*raw.h5seurat"), emit: h5seurat
    tuple val(meta), path("*_DEG.tsv")    , optional: true, emit: DEGlist

    script:
    // Different input files names when including multi-gene reasds
    def summaryFile = params.soloMultiMappers == "Unique" ? "Summary.csv" : "Summary.multiple.csv"
    def matrixDir = params.soloMultiMappers == "Unique" ? "filtered" : "filtered_mult"
    """
    Rscript -e 'rmarkdown::render("$baseDir/bin/scRNA_report.Rmd", params = list(sampleName = "${meta.id}", starsolo_out = "${starsolo_summary}", qualimap_out = "$qualimap_outdir/rnaseq_qc_results.txt", qualimap_gene_coverage = "$qualimap_outdir/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt", starsolo_bc = "$starsolo_UMI_file", starsolo_matrixDir="${starsolo_filteredDir}", nCPUs = "$task.cpus", saturation_json = "${saturation_outJSON}", version_json = "${version_json}"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${meta.id}_report.html")'

    """
}