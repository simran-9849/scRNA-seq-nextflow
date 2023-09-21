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

    def gex_cells           = filteredDir_featureTypes.contains("GEX") ? "${filteredDir_map['GEX']}/barcodes.tsv.gz" : "None"
    def featureArgument     = report_map.findAll{it.key != "GEX"}.keySet()
    def summaryArgument     = featureArgument.collect{ starsoloSummary_map[it] }.join(",")
    featureArgument         = featureArgument.join(",")
    """
    
    ## VDJ-B
    trust4_metrics_filter.sh ${meta.id} \\
                             ${featureArgument} \\
                             ${gex_cells} \\
                             ${summaryArgument}

    """
}