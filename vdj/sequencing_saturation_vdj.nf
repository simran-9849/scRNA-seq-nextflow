process CHECK_SATURATION_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_high'
    publishDir "${params.outdir}/saturation/",
        mode: "${params.publish_dir_mode}",
        enabled: params.publishSaturation as boolean

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsolo_filteredDir)
    path(whitelist)

    output:
    tuple val(meta), path("${meta.id}_${meta.feature_types}.saturation_out.json"), emit: outJSON

    script:
    def multiMapper = params.soloMultiMappers == "Unique" ? "unique" : "multiple"
    if(meta.feature_types =~ /GEX/){
    """
    cellFile=\$(mktemp -p ./)
    zcat ${starsolo_filteredDir}/barcodes.tsv.gz > \$cellFile
    get_sequencing_saturation.sh ${whitelist} \$cellFile ${multiMapper} ${starsoloBAM} ${task.cpus} \\
                                 ${meta.id}_${meta.feature_types}.saturation_data.json \\
                                 ${meta.id}_${meta.feature_types}.UMI_hist.tsv \\
                                 ${meta.id}_${meta.feature_types}.gene_hist.tsv \\
                                 ${meta.id}_${meta.feature_types}.totalGeneCount.tsv
    rm \$cellFile
    combine_saturation_data.R ${meta.id}_${meta.feature_types} \\
                              ${meta.id}_${meta.feature_types}.saturation_data.json \\
                              ${meta.id}_${meta.feature_types}.UMI_hist.tsv \\
                              ${meta.id}_${meta.feature_types}.gene_hist.tsv \\
                              ${meta.id}_${meta.feature_types}.totalGeneCount.tsv \\
                              ${meta.id}_${meta.feature_types}.saturation_out.json
    """
    }else{
    """
    touch ${meta.id}_${meta.feature_types}.saturation_out.json
    """
    }
}
