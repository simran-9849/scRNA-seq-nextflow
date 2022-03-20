process CHECK_SATURATION {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsolo_outdir)
    path(whitelist)

    output:
    tuple val(meta), path("${meta.id}.saturation_out.json"), emit: outJSON

    shell:
    '''
    cellFile=$(mktemp -p ./)
    zcat !{starsolo_outdir}/filtered/barcodes.tsv.gz > $cellFile
    get_sequencing_saturation.sh !{whitelist} $cellFile !{starsoloBAM} !{task.cpus} !{meta.id}.saturation_out.json
    '''
}
