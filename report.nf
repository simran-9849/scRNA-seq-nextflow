#! /usr/bin/env nextflow

process REPORT{
    tag "$meta.id"
    label 'process_medium'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(starsolo_summary)
    tuple val(meta), path(starsolo_UMI_file)
    tuple val(meta), path(starsolo_filteredDir)
    tuple val(meta), path(featureStats)
    tuple val(meta), path(geneCoverage)
    tuple val(meta), path(saturation_outJSON)
    tuple val(meta), path(version_json)

    output:
    tuple val(meta), path("*report.html") , emit: report
    tuple val(meta), path("*metrics.json"), emit: metrics_json
    tuple val(meta), path("*metrics.tsv"),  emit: metrics_tsv
    tuple val(meta), path("*_DEG.tsv")    , optional: true, emit: DEGlist

    script:
    // Different input files names when including multi-gene reasds
    def summaryFile = params.soloMultiMappers == "Unique" ? "Summary.csv" : "Summary.multiple.csv"
    def matrixDir = params.soloMultiMappers == "Unique" ? "filtered" : "filtered_mult"
    def nMem = "${task.memory.toBytes()}"
    """
    #! /usr/bin/env Rscript

    rmarkdown::render(
        "$baseDir/bin/scRNA_report.Rmd",
        params = list(
            sampleName = "${meta.id}",
            starsolo_out = "${starsolo_summary}",
            featureStats = "${meta.id}.featureCoverage_stats.json",
            geneCoverage = "${meta.id}.scaled.tab",
            starsolo_bc = "$starsolo_UMI_file",
            starsolo_matrixDir="${starsolo_filteredDir}",
            nCPUs = "$task.cpus",
            nMem = "${nMem}",
            saturation_json = "${saturation_outJSON}", version_json = "${version_json}"
        ),
        intermediates_dir = getwd(),
        knit_root_dir = getwd(),
        output_dir = getwd(),
        output_file = "${meta.id}_report.html"
    )
    """
}