process REPORT_VDJ {
    tag "${meta.id}"
    label 'process_medium'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), val(starsolo_summary_featureTypes),     path(starsolo_summary)
    tuple val(meta), val(starsolo_UMI_file_featureTypes),    path(starsolo_UMI_file)
    tuple val(meta), val(starsolo_filteredDir_featureTypes), path(starsolo_filteredDir)
    tuple val(meta), val(qualimap_outdir_featureTypes),      path(qualimap_outdir)
    tuple val(meta), val(saturation_outJSON_featureTypes),   path(saturation_outJSON)
    tuple val(meta), val(vdj_report_featureTypes),           path(vdj_report)
    tuple val(meta), val(vdj_airr_featureTypes),             path(vdj_airr)
    tuple val(meta), val(vdj_kneeOut_featureTypes),          path(vdj_kneeOut)
    tuple val(meta), val(vdj_finalOut_featureTypes),         path(vdj_finalOut)
    tuple val(meta), val(trust4_cells_featureTypes),         path(trust4_cells)
    tuple val(meta), val(trust4_metrics_featureTypes),       path(trust4_metrics)
    tuple val(meta), val(trust4_cloneType_featureTypes),     path(trust4_cloneType)
    path(version_json)    

    output:
    tuple val(meta), path("*report.html") , emit: report
    //tuple val(meta), path("*metrics.json"), emit: metrics
    //tuple val(meta), path("*raw.h5seurat"), emit: h5seurat
    tuple val(meta), path("*_DEG.tsv")    , optional: true, emit: DEGlist
    tuple val(meta), path("${meta.id}_GEX.Summary.unique.csv"), optional: true, emit: GEX_Summary
    tuple val(meta), path("${meta.id}_*.metrics.tsv"), emit: metrics_tsv 

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        dList = []
        if(feature_types.size() == 1 && data_list.getClass() == nextflow.processor.TaskPath){
            dList = [data_list]
        }else{
            dList = data_list
        }
        map = [feature_types, dList].transpose().collectEntries()
        //if(feature_types.size() > 1){
        //map = [feature_types, data_list].transpose().collectEntries()
        //}else if(feature_types.size() == 1){
        //    map = [feature_types, [data_list]].transpose().collectEntries()
        //}//else{
            // return empty, them program will stop and throw the warnning
            //def map = [:]
        //}
        return map
    }

    def starsolo_summary_map = associate_feature_type(starsolo_summary_featureTypes, starsolo_summary)
    def starsolo_UMI_file_map = associate_feature_type(starsolo_UMI_file_featureTypes, starsolo_UMI_file)
    def starsolo_filteredDir_map = associate_feature_type(starsolo_filteredDir_featureTypes, starsolo_filteredDir)
    def qualimap_outdir_map = associate_feature_type(qualimap_outdir_featureTypes, qualimap_outdir)
    def saturation_outJSON_map = associate_feature_type(saturation_outJSON_featureTypes, saturation_outJSON)
    def vdj_report_map   = associate_feature_type(vdj_report_featureTypes, vdj_report)
    def vdj_airr_map     = associate_feature_type(vdj_airr_featureTypes, vdj_airr)
    def vdj_kneeOut_map  = associate_feature_type(vdj_kneeOut_featureTypes, vdj_kneeOut)
    def vdj_finalOut_map = associate_feature_type(vdj_finalOut_featureTypes, vdj_finalOut)
    def trust4_cells_map = associate_feature_type(trust4_cells_featureTypes, trust4_cells)
    def trust4_metrics_map = associate_feature_type(trust4_metrics_featureTypes, trust4_metrics)
    def trust4_cloneType_map = associate_feature_type(trust4_cloneType_featureTypes, trust4_cloneType)

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

    // VDJ_B inputs
    def VDJ_B_report   = vdj_report_map["VDJ-B"] ? vdj_report_map["VDJ-B"] : ""
    def VDJ_B_airr     = vdj_airr_map["VDJ-B"] ? vdj_airr_map["VDJ-B"] : ""
    def VDJ_B_kneeOut  = vdj_kneeOut_map["VDJ-B"] ? vdj_kneeOut_map["VDJ-B"] : ""
    def VDJ_B_finalOut = vdj_finalOut_map["VDJ-B"] ? vdj_finalOut_map["VDJ-B"] : ""
    def VDJ_B_cells    = trust4_cells_map["VDJ-B"] ? trust4_cells_map["VDJ-B"] : ""
    def VDJ_B_trust4_metrics = trust4_metrics_map["VDJ-B"] ? trust4_metrics_map["VDJ-B"] : ""
    def VDJ_B_cloneType = trust4_cloneType_map["VDJ-B"] ? trust4_cloneType_map["VDJ-B"] : ""

    // VDJ_T inputs
    def VDJ_T_report   = vdj_report_map["VDJ-T"] ? vdj_report_map["VDJ-T"] : ""
    def VDJ_T_airr     = vdj_airr_map["VDJ-T"] ? vdj_airr_map["VDJ-T"] : ""
    def VDJ_T_kneeOut  = vdj_kneeOut_map["VDJ-T"] ? vdj_kneeOut_map["VDJ-T"] : ""
    def VDJ_T_finalOut = vdj_finalOut_map["VDJ-T"] ? vdj_finalOut_map["VDJ-T"] : ""
    def VDJ_T_cells    = trust4_cells_map["VDJ-T"] ? trust4_cells_map["VDJ-T"] : ""
    def VDJ_T_trust4_metrics = trust4_metrics_map["VDJ-T"] ? trust4_metrics_map["VDJ-T"] : ""
    def VDJ_T_cloneType = trust4_cloneType_map["VDJ-T"] ? trust4_cloneType_map["VDJ-T"] : ""

    // Using UMI or Reads
    def withUMI = params.soloType == "CB_samTagOut" ? "FALSE" : "TRUE"
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
