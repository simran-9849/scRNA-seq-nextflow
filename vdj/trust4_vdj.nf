process TRUST4_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_high'
    publishDir "${params.outdir}/trust4/${meta.id}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/_GEX/){
            return null
        }else if(filename=~/Solo.out/){
            return filename.split("/")[-1]
        }else{
            return filename
        }
    }

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsoloSummary)
    path(trust4_vdj_refGenome_fasta)
    path(trust4_vdj_imgt_fasta)
    path(whitelist)

    output:
    tuple val(meta), path("TRUST_${meta.id}_${meta.feature_types}_toassemble_bc.fa"), emit: toassemble_bc
    tuple val(meta), path("TRUST_${meta.id}_${meta.feature_types}_barcode_report.tsv"), emit: trust4_report
    tuple val(meta), path("TRUST_${meta.id}_${meta.feature_types}_barcode_airr.tsv"), emit: trust4_airr

    script:
    if(meta.feature_types =~ /VDJ-[TB]/){
    """
    run-trust4 -f ${trust4_vdj_refGenome_fasta} --ref ${trust4_vdj_imgt_fasta} \\
    -b ${starsoloBAM} \\
    --barcode CB \\
    --barcodeWhitelist $whitelist \\
    --UMI UB \\
    --outputReadAssignment \\
    -t ${task.cpus}
    """
    }else{
    """
    touch TRUST_${meta.id}_${meta.feature_types}_toassemble_bc.fa
    touch TRUST_${meta.id}_${meta.feature_types}_barcode_report.tsv
    touch TRUST_${meta.id}_${meta.feature_types}_barcode_airr.tsv
    """
    }
}

process VDJ_METRICS {
    tag "${meta.id}:VDJ_metrics"
    label 'process_low'
    publishDir "${params.outdir}/trust4/${meta.id}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/_GEX/){
            return null
        }else if(filename=~/Solo.out/){
            return filename.split("/")[-1]
        }else{
            return filename
        }
    }

    input:
    tuple val(meta), val(report_featureTypes),          path(report)
    tuple val(meta), val(airr_featureTypes),            path(airr)
    tuple val(meta), val(toassemble_bc_featureTypes),   path(toassemble_bc)
    tuple val(meta), val(starsoloSummary_featureTypes), path(starsoloSummary)
    tuple val(meta), val(filteredDir_featureTypes),     path(filteredDir)

    output:
    tuple val(meta), path("${meta.id}_*.vdj_metrics.json"), emit: metricsJSON
    tuple val(meta), path("${meta.id}_*.knee_input.tsv"), emit: kneeData
    tuple val(meta), path("${meta.id}_*.cloneType_out.tsv"), emit: cloneType
    tuple val(meta), path("TRUST_${meta.id}_*_barcode_report.filtered.tsv"), emit: trust4_report
    tuple val(meta), path("TRUST_${meta.id}_*_barcode_airr.filtered.tsv"), emit: trust4_airr

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        if(feature_types.size() > 1){
            map = [feature_types, data_list].transpose().collectEntries()
        }else if(feature_types.size() == 1){
            map = [feature_types, [data_list]].transpose().collectEntries()
        }//else{
            // return empty, the program will stop and throw the warnning
            //def map = [:]
        //}
        return map
    }
    def report_map          = associate_feature_type(report_featureTypes, report)
    def airr_map            = associate_feature_type(airr_featureTypes, airr)
    def toassemble_bc_map   = associate_feature_type(toassemble_bc_featureTypes, toassemble_bc)
    def starsoloSummary_map = associate_feature_type(starsoloSummary_featureTypes, starsoloSummary)
    def filteredDir_map     = associate_feature_type(filteredDir_featureTypes, filteredDir)

    if(filteredDir_featureTypes.contains("GEX")){
        if(report_featureTypes.contains("VDJ-B") && report_featureTypes.contains("VDJ-T")){
        """
        zcat ${filteredDir_map["GEX"]}/barcodes.tsv.gz > gex_cells.tsv
        ## VDJ-B
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-B_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-B"]} \\
                                 VDJ-B
        mv VDJ-B.vdj_metrics.json ${meta.id}_VDJ-B.vdj_metrics.json
        mv VDJ-B.knee_input.tsv ${meta.id}_VDJ-B.knee_input.tsv
        mv VDJ-B.cloneType_out.tsv ${meta.id}_VDJ-B.cloneType_out.tsv
        ## VDJ-T
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-T_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-T"]} \\
                                 VDJ-T
        mv VDJ-T.vdj_metrics.json ${meta.id}_VDJ-T.vdj_metrics.json
        mv VDJ-T.knee_input.tsv ${meta.id}_VDJ-T.knee_input.tsv
        mv VDJ-T.cloneType_out.tsv ${meta.id}_VDJ-T.cloneType_out.tsv
        """
        }else if(report_featureTypes.contains("VDJ-B") && !report_featureTypes.contains("VDJ-T")){
        """
        ## Filter trust4 result by gex cells
        zcat ${filteredDir_map["GEX"]}/barcodes.tsv.gz > gex_cells.tsv
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-B_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-B"]} \\
                                 VDJ-B
        mv VDJ-B.vdj_metrics.json ${meta.id}_VDJ-B.vdj_metrics.json
        mv VDJ-B.knee_input.tsv ${meta.id}_VDJ-B.knee_input.tsv
        mv VDJ-B.cloneType_out.tsv ${meta.id}_VDJ-B.cloneType_out.tsv
        """
        }else if(!report_featureTypes.contains("VDJ-B") && report_featureTypes.contains("VDJ-T")){
        """
        ## Filter trust4 result by gex cells
        zcat ${filteredDir_map["GEX"]}/barcodes.tsv.gz > gex_cells.tsv
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-T_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-T"]} \\
                                 VDJ-T
        mv VDJ-T.vdj_metrics.json ${meta.id}_VDJ-T.vdj_metrics.json
        mv VDJ-T.knee_input.tsv ${meta.id}_VDJ-T.knee_input.tsv
        mv VDJ-T.cloneType_out.tsv ${meta.id}_VDJ-T.cloneType_out.tsv
        """
        }
    }else{
        if(report_featureTypes.contains("VDJ-B") && report_featureTypes.contains("VDJ-T")){
        """
        ## No GEX, no filter, use all the barcode in the report as "gex_cells"
        ## VDJ-B
        awk 'NR>1 && \$1!="-"{print \$1}' TRUST_${meta.id}_VDJ-B_barcode_report.tsv > gex_cells.tsv
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-B_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-B"]} \\
                                 VDJ-B
        mv VDJ-B.vdj_metrics.json ${meta.id}_VDJ-B.vdj_metrics.json
        mv VDJ-B.knee_input.tsv ${meta.id}_VDJ-B.knee_input.tsv
        mv VDJ-B.cloneType_out.tsv ${meta.id}_VDJ-B.cloneType_out.tsv
        ## VDJ-T
        awk 'NR>1 && \$1!="-"{print \$1}' TRUST_${meta.id}_VDJ-T_barcode_report.tsv > gex_cells.tsv
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-T_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-T"]} \\
                                 VDJ-T
        mv VDJ-T.vdj_metrics.json ${meta.id}_VDJ-T.vdj_metrics.json
        mv VDJ-T.knee_input.tsv ${meta.id}_VDJ-T.knee_input.tsv
        mv VDJ-T.cloneType_out.tsv ${meta.id}_VDJ-T.cloneType_out.tsv
        """
        }else if(report_featureTypes.contains("VDJ-B") && !report_featureTypes.contains("VDJ-T")){
        """
        ## No GEX, no filter, use all the barcode in the report as "gex_cells"
        awk 'NR>1 && \$1!="-"{print \$1}' TRUST_${meta.id}_VDJ-B_barcode_report.tsv > gex_cells.tsv
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-B_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-B_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-B"]} \\
                                 VDJ-B
        mv VDJ-B.vdj_metrics.json ${meta.id}_VDJ-B.vdj_metrics.json
        mv VDJ-B.knee_input.tsv ${meta.id}_VDJ-B.knee_input.tsv
        mv VDJ-B.cloneType_out.tsv ${meta.id}_VDJ-B.cloneType_out.tsv
        """
        }else if(!report_featureTypes.contains("VDJ-B") && report_featureTypes.contains("VDJ-T")){
        """
        ## No GEX, no filter, use all the barcode in the report as "gex_cells"
        awk 'NR>1 && \$1!="-"{print \$1}' TRUST_${meta.id}_VDJ-T_barcode_report.tsv > gex_cells.tsv
        trust4_metrics_filter.sh gex_cells.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_report.tsv \\
                                 TRUST_${meta.id}_VDJ-T_barcode_airr.tsv \\
                                 TRUST_${meta.id}_VDJ-T_toassemble_bc.fa \\
                                 ${starsoloSummary_map["VDJ-T"]} \\
                                 VDJ-T
        mv VDJ-T.vdj_metrics.json ${meta.id}_VDJ-T.vdj_metrics.json
        mv VDJ-T.knee_input.tsv ${meta.id}_VDJ-T.knee_input.tsv
        mv VDJ-T.cloneType_out.tsv ${meta.id}_VDJ-T.cloneType_out.tsv
        """
        }
    }

}