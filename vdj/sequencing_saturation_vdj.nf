process CHECK_SATURATION_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_high'

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsolo_filteredDir)
    path(whitelist)

    output:
    tuple val(meta), path("${meta.id}_${meta.feature_types}.saturation_out.json"), emit: outJSON

    script:
    def multiMapper = params.soloMultiMappers == "Unique" ? "unique" : "multiple"
    """
    cellFile=\$(mktemp -p ./)
    zcat ${starsolo_filteredDir}/barcodes.tsv.gz > \$cellFile
    get_sequencing_saturation.sh ${whitelist} \$cellFile ${multiMapper} ${starsoloBAM} ${task.cpus} ${meta.id}_${meta.feature_types}.saturation_out.json
    rm \$cellFile
    """
}
