process CHECK_SATURATION {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/saturation/",
        mode: "${params.publish_dir_mode}",
        enabled: params.publishSaturation as boolean

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
    get_sequencing_saturation.sh ${whitelist} \$cellFile ${multiMapper} ${starsoloBAM} ${task.cpus} ${meta.id}.saturation_data.json ${meta.id}.UMI_hist.tsv ${meta.id}.gene_hist.tsv ${meta.id}.totalGeneCount.tsv
    rm \$cellFile
    combine_saturation_data.R ${meta.id} ${meta.id}.saturation_data.json ${meta.id}.UMI_hist.tsv ${meta.id}.gene_hist.tsv ${meta.id}.totalGeneCount.tsv ${meta.id}.saturation_out.json
    """
}
