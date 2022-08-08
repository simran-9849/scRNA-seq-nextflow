process CHECK_SATURATION {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsolo_filteredDir)
    path(whitelist)

    output:
    tuple val(meta), path("${meta.id}.saturation_out.json"), emit: outJSON

    script:
    def multiMapper = params.soloMultiMappers == "Unique" ? "unique" : "multiple"
    """
    cellFile=\$(mktemp -p ./)
    zcat ${starsolo_filteredDir}/barcodes.tsv.gz > \$cellFile
    get_sequencing_saturation.sh ${whitelist} \$cellFile ${multiMapper} ${starsoloBAM} ${task.cpus} ${meta.id}.saturation_out.json
    rm \$cellFile
    """
}
