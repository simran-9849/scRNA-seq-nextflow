process REPORT_VDJ {
    tag "${meta.id}:VDJ_report"
    label 'process_medium'
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), val(starsolo_summary_feature_types), path(starsolo_summary)
    tuple val(meta), val(starsolo_UMI_file_feature_types), path(starsolo_UMI_file)
    tuple val(meta), val(starsolo_filteredDir_feature_types), path(starsolo_filteredDir)
    tuple val(meta), val(qualimap_outdir_feature_types), path(qualimap_outdir)
    tuple val(meta), val(saturation_outJSON_feature_types), path(saturation_outJSON)
    tuple val(meta), val(trust4_metrics_feature_types), path(trust4_metrics)
    tuple val(meta), val(trust4_kneeData_feature_types), path(trust4_kneeData)
    tuple val(meta), val(trust4_cloneType_feature_types), path(trust4_cloneType)
    path(version_json)    

    output:
    tuple val(meta), path("*report.html") , emit: report
    //tuple val(meta), path("*metrics.json"), emit: metrics
    //tuple val(meta), path("*raw.h5seurat"), emit: h5seurat
    //tuple val(meta), path("*_DEG.tsv")    , optional: true, emit: DEGlist

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        if(feature_types.size() > 1){
            map = [feature_types, data_list].transpose().collectEntries()
        }else if(feature_types.size() == 1){
            map = [feature_types, [data_list]].transpose().collectEntries()
        }//else{
            // return empty, them program will stop and throw the warnning
            //def map = [:]
        //}
        return map
    }

    def starsolo_summary_map = associate_feature_type(starsolo_summary_feature_types, starsolo_summary)
    def starsolo_UMI_file_map = associate_feature_type(starsolo_UMI_file_feature_types, starsolo_UMI_file)
    def starsolo_filteredDir_map = associate_feature_type(starsolo_filteredDir_feature_types, starsolo_filteredDir)
    def qualimap_outdir_map = associate_feature_type(qualimap_outdir_feature_types, qualimap_outdir)
    def saturation_outJSON_map = associate_feature_type(saturation_outJSON_feature_types, saturation_outJSON)
    def trust4_metrics_map = associate_feature_type(trust4_metrics_feature_types, trust4_metrics)
    def trust4_kneeData_map = associate_feature_type(trust4_kneeData_feature_types, trust4_kneeData)
    def trust4_cloneType_map = associate_feature_type(trust4_cloneType_feature_types, trust4_cloneType)

    // Different input files names when including multi-gene reasds
    //def summaryFile = params.soloMultiMappers == "Unique" ? "Summary.csv" : "Summary.multiple.csv"
    //def matrixDir = params.soloMultiMappers == "Unique" ? "filtered" : "filtered_mult"
    def GEX_summaryFile = starsolo_summary_map["GEX"] ?  starsolo_summary_map["GEX"] : ""
    //def GEX_qualimapDir = qualimap_outdir_map["GEX"] ? qualimap_outdir_map["GEX"] : ""
    def qualimapInput =  qualimap_outdir_map["GEX"] ? qualimap_outdir_map["GEX"] + "/rnaseq_qc_results.txt" : ""
    def qualimapGeneCoverage = qualimap_outdir_map["GEX"] ? qualimap_outdir_map["GEX"] + "/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" : ""
    def GEX_UMI_file = starsolo_UMI_file_map["GEX"] ? starsolo_UMI_file_map["GEX"] : ""
    def GEX_matrixDir = starsolo_filteredDir_map["GEX"] ? starsolo_filteredDir_map["GEX"] : ""
    def GEX_saturation = saturation_outJSON_map["GEX"] ? saturation_outJSON_map["GEX"] : ""
    def VDJ_B_trust4_metrics = trust4_metrics_map["VDJ-B"] ? trust4_metrics_map["VDJ-B"] : ""
    def VDJ_T_trust4_metrics = trust4_metrics_map["VDJ-T"] ? trust4_metrics_map["VDJ-T"] : ""
    def VDJ_B_kneeData = trust4_kneeData_map["VDJ-B"] ? trust4_kneeData_map["VDJ-B"] : ""
    def VDJ_T_kneeData = trust4_kneeData_map["VDJ-T"] ? trust4_kneeData_map["VDJ-T"] : ""
    def VDJ_B_cloneType = trust4_cloneType_map["VDJ-B"] ? trust4_cloneType_map["VDJ-B"] : ""
    def VDJ_T_cloneType = trust4_cloneType_map["VDJ-T"] ? trust4_cloneType_map["VDJ-T"] : "" 


    """
    #! /usr/bin/env Rscript

    rmarkdown::render(
        "$baseDir/bin/scRNA_vdj_gex_report.Rmd",
        params = list(
            sampleName = "${meta.id}",
            starsolo_out = "${GEX_summaryFile}",
            qualimap_out = "${qualimapInput}",
            qualimap_gene_coverage = "${qualimapGeneCoverage}",
            starsolo_bc = "${GEX_UMI_file}",
            starsolo_matrixDir="${GEX_matrixDir}",
            nCPUs = "$task.cpus",
            saturation_json = "${GEX_saturation}",
            version_json = "${version_json}",
            VDJ_T_metrics = "${VDJ_T_trust4_metrics}",
            VDJ_B_metrics = "${VDJ_B_trust4_metrics}",
            VDJ_T_kneeData = "${VDJ_T_kneeData}",
            VDJ_B_kneeData = "${VDJ_B_kneeData}",
            VDJ_T_cloneType = "${VDJ_T_cloneType}",
            VDJ_B_cloneType = "${VDJ_B_cloneType}"
        ),
        intermediates_dir = getwd(), knit_root_dir = getwd(),
        output_dir = getwd(),
        output_file = "${meta.id}_VDJ_report.html"
    )
    """
}
