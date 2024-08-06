process REPORT_VDJ {
    tag "${meta.id}"
    label 'process_medium'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/${meta.id}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta),
          val(starsolo_summary),
          val(starsolo_UMI_file),
          val(starsolo_filteredDir),
          val(featureStats),
          val(geneCoverage),
          val(saturation_outJSON),
          val(vdj_report),
          val(vdj_airr),
          val(vdj_kneeOut),
          val(vdj_finalOut),
          val(trust4_cells),
          val(trust4_metrics),
          val(trust4_cloneType)
    path(version_json)    

    output:
    tuple val(meta), path("*report.html") ,                                         emit: report
    //tuple val(meta), path("*metrics.json"),                                       emit: metrics
    //tuple val(meta), path("*raw.h5seurat"),                                       emit: h5seurat
    //tuple val(meta), path("${starsolo_filteredDir}"),               optional: true, emit: gex_filteredDir
    tuple val(meta), path("*_DEG.tsv"),                             optional: true, emit: DEGlist
    tuple val(meta), path("${meta.id}_GEX.Summary.unique.csv"),     optional: true, emit: GEX_Summary
    tuple val(meta), path("${meta.id}_*.metrics.tsv"),                              emit: metrics_tsv
    tuple val(meta), path("${meta.id}_*.metrics.json"),                             emit: metrics_json
    tuple val(meta), path("${meta.id}*_clonotypes.tsv"),                            emit: vdj_clonotype
    tuple val(meta), path("${meta.id}*_results.tsv"),                               emit: vdj_results
    tuple val(meta), path("${meta.id}*_results.productiveOnly_withLineage.tsv"),    emit: vdj_lineage
    tuple val(meta), path("${version_json}"),                       optional: true, emit: versions

    script:
    // Different input files names when including multi-gene reasds
    //def summaryFile = params.soloMultiMappers == "Unique" ? "Summary.csv" : "Summary.multiple.csv"
    //def matrixDir = params.soloMultiMappers == "Unique" ? "filtered" : "filtered_mult"
    def GEX_summaryFile = starsolo_summary
    def featureStatsFile = featureStats
    def geneCoverageFile = geneCoverage
    def GEX_UMI_file = starsolo_UMI_file
    def GEX_matrixDir = starsolo_filteredDir
    def GEX_saturation = saturation_outJSON

    // VDJ_B inputs
    def VDJ_B_report   = vdj_report[0]
    def VDJ_B_airr     = vdj_airr[0]
    def VDJ_B_kneeOut  = vdj_kneeOut[0]
    def VDJ_B_finalOut = vdj_finalOut[0]
    def VDJ_B_cells    = trust4_cells[0]
    def VDJ_B_trust4_metrics = trust4_metrics[0]
    def VDJ_B_cloneType = trust4_cloneType[0]

    // VDJ_T inputs
    def VDJ_T_report   = vdj_report[1]
    def VDJ_T_airr     = vdj_airr[1]
    def VDJ_T_kneeOut  = vdj_kneeOut[1]
    def VDJ_T_finalOut = vdj_finalOut[1]
    def VDJ_T_cells    = trust4_cells[1]
    def VDJ_T_trust4_metrics = trust4_metrics[1]
    def VDJ_T_cloneType = trust4_cloneType[1]

    // Using UMI or Reads
    def withUMI = params.soloType == "CB_samTagOut" ? "FALSE" : "TRUE"
    """
    #! /usr/bin/env Rscript

    rmarkdown::render(
        "$baseDir/bin/scRNA_vdj_gex_report.Rmd",
        params = list(
            sampleName = "${meta.id}",
            starsolo_out = "${GEX_summaryFile}",
            featureStats = "${featureStatsFile}",
            geneCoverage = "${geneCoverageFile}",
            starsolo_bc = "${GEX_UMI_file}",
            starsolo_matrixDir="${GEX_matrixDir}",
            nCPUs = "$task.cpus",
            saturation_json = "${GEX_saturation}",
            version_json = "${version_json}",
            VDJ_B_report = "${VDJ_B_report}",
            VDJ_B_airr = "${VDJ_B_airr}",
            VDJ_B_kneeOut = "${VDJ_B_kneeOut}",
            VDJ_B_finalOut = "${VDJ_B_finalOut}",
            VDJ_B_cells = "${VDJ_B_cells}",
            VDJ_B_metrics = "${VDJ_B_trust4_metrics}",
            VDJ_B_cloneType = "${VDJ_B_cloneType}",
            VDJ_T_report = "${VDJ_T_report}",
            VDJ_T_airr = "${VDJ_T_airr}",
            VDJ_T_kneeOut = "${VDJ_T_kneeOut}",
            VDJ_T_finalOut = "${VDJ_T_finalOut}",
            VDJ_T_cells = "${VDJ_T_cells}",
            VDJ_T_metrics = "${VDJ_T_trust4_metrics}",
            VDJ_T_cloneType = "${VDJ_T_cloneType}",
            withUMI = ${withUMI}
        ),
        intermediates_dir = getwd(),
        knit_root_dir = getwd(),
        output_dir = getwd(),
        output_file = "${meta.id}_VDJ_report.html"
    )
    """
}
