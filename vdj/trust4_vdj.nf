process TRUST4_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_high'
    publishDir "${params.outdir}/trust4/${meta.id}_${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsoloSummary)
    path(trust4_vdj_refGenome_fasta)
    path(trust4_vdj_imgt_fasta)
    path(whitelist)

    output:
    tuple val(meta), path("${meta.id}_${meta.feature_types}.vdj_metrics.json"), emit: metricsJSON
    tuple val(meta), path("${meta.id}_${meta.feature_types}.knee_input.tsv"), emit: kneeData
    tuple val(meta), path("${meta.id}_${meta.feature_types}.cloneType_out.tsv"), emit: cloneType

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
    
    trust4_metrics.sh ${meta.id}_${meta.feature_types} \\
    ${meta.feature_types} \\
    ${starsoloSummary} \\
    ${meta.id}_${meta.feature_types}.vdj_metrics.json \\
    ${meta.id}_${meta.feature_types}.knee_input.tsv \\
    ${meta.id}_${meta.feature_types}.cloneType_out.tsv
    """
    }else{
    """
    touch ${meta.id}_${meta.feature_types}.vdj_metrics.json
    touch ${meta.id}_${meta.feature_types}.knee_input.tsv
    touch ${meta.id}_${meta.feature_types}.cloneType_out.tsv
    """
    }
}