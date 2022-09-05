process REPORT_VDJ {
    tag "${meta.id}:VDJ_report"
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), val(feature_types), path(starsolo_summary)
    tuple val(meta), val(feature_types), path(starsolo_UMI_file)
    tuple val(meta), val(feature_types), path(starsolo_filteredDir)
    tuple val(meta), val(feature_types), path(qualimap_outdir)
    tuple val(meta), val(feature_types), path(saturation_outJSON)
    tuple val(meta), val(feature_types), path(trust4_metrics)
    tuple val(meta), val(feature_types), path(trust4_kneeData)
    tuple val(meta), val(feature_types), path(trust4_cloneType)
    path(version_json)    

    output:
    tuple val(meta), path("*report.html") , emit: report
    //tuple val(meta), path("*metrics.json"), emit: metrics
    //tuple val(meta), path("*raw.h5seurat"), emit: h5seurat
    //tuple val(meta), path("*_DEG.tsv")    , optional: true, emit: DEGlist

    script:
    // Different input files names when including multi-gene reasds
    def summaryFile = params.soloMultiMappers == "Unique" ? "Summary.csv" : "Summary.multiple.csv"
    def matrixDir = params.soloMultiMappers == "Unique" ? "filtered" : "filtered_mult"
    def GEX_summaryFile = ''
    def GEX_qualimapDir = ''
    def GEX_UMI_file = ''
    def GEX_matrixDir = ''
    def GEX_saturation = ''
    def VDJ_B_trust4_metrics = ''
    def VDJ_T_trust4_metrics = ''
    def VDJ_B_kneeData = ''
    def VDJ_T_kneeData = ''
    def VDJ_B_cloneType = ''
    def VDJ_T_cloneType = ''

    def feature_types_list = feature_types.collect{ it.toString() }
    def starsolo_summary_list = starsolo_summary.collect{ it.toString() }
    def starsolo_UMI_file_list = starsolo_UMI_file.collect{ it.toString() }
    def starsolo_filteredDir_list = starsolo_filteredDir.collect{ it.toString() }
    def qualimap_outdir_list = qualimap_outdir.collect{ it.toString() }
    def saturation_outJSON_list = saturation_outJSON.collect{ it.toString() }
    def trust4_metrics_list = trust4_metrics.collect{ it.toString() }
    def trust4_kneeData_list = trust4_kneeData.collect{ it.toString() }
    def trust4_cloneType_list = trust4_cloneType.collect{ it.toString() }
    feature_types_list.eachWithIndex{ v, ix ->
      if(v == "GEX"){
          GEX_summaryFile = starsolo_summary_list[ix]
          GEX_qualimapDir = qualimap_outdir_list[ix]
          GEX_UMI_file = starsolo_UMI_file_list[ix]
          GEX_matrixDir = starsolo_filteredDir_list[ix]
          GEX_saturation = saturation_outJSON_list[ix]
      }else if(v == "VDJ-B"){
          VDJ_B_trust4_metrics = trust4_metrics_list[ix]
          VDJ_B_kneeData = trust4_kneeData_list[ix]
          VDJ_B_cloneType = trust4_cloneType_list[ix]
      }else if(v == "VDJ-T"){
          VDJ_T_trust4_metrics = trust4_metrics_list[ix]
          VDJ_T_kneeData = trust4_kneeData_list[ix]
          VDJ_T_cloneType = trust4_cloneType_list[ix]
      }
    }
    """
    Rscript -e 'rmarkdown::render("$baseDir/bin/scRNA_vdj_gex_report.Rmd", params = list(sampleName = "${meta.id}", starsolo_out = "${GEX_summaryFile}", qualimap_out = "${GEX_qualimapDir}/rnaseq_qc_results.txt", qualimap_gene_coverage = "${GEX_qualimapDir}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt", starsolo_bc = "${GEX_UMI_file}", starsolo_matrixDir="${GEX_matrixDir}", nCPUs = "$task.cpus", saturation_json = "${GEX_saturation}", version_json = "${version_json}", VDJ_T_metrics = "${VDJ_T_trust4_metrics}", VDJ_B_metrics = "${VDJ_B_trust4_metrics}", VDJ_T_kneeData = "${VDJ_T_kneeData}", VDJ_B_kneeData = "${VDJ_B_kneeData}", VDJ_T_cloneType = "${VDJ_T_cloneType}", VDJ_B_cloneType = "${VDJ_B_cloneType}"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${meta.id}_VDJ_report.html")'
    """
}
