#! /usr/bin/env nextflow
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './modules/nf-core_rnaseq/functions'

params.options = [:]
options        = initOptions(params.options)

process REPORT{
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'report', meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(starsolo_outdir)
    tuple val(meta), path(qualimap_outdir)

    output:
    tuple val(meta), path("*report.html"), emit: report

    script:
    """
    Rscript -e 'rmarkdown::render("$baseDir/bin/scRNA_report.Rmd", params = list(starsolo_out = "$starsolo_outdir/Summary.csv", qualimap_out = "$qualimap_outdir/rnaseq_qc_results.txt", starsolo_bc = "$starsolo_outdir/UMIperCellSorted.txt"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${meta.id}_report.html")'

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
    ${getSoftwareName(task.process)}: \$(Rscript -e "sessionInfo()")
    END_VERSIONS
    """
}