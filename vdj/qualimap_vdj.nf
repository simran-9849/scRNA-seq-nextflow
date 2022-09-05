process QUALIMAP_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_medium'
    publishDir "${params.outdir}/qualimap",
        mode: params.publish_dir_mode,
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("${prefix}"), emit: results

    script:
    prefix         = "${meta.id}_${meta.feature_types}"
    def paired_end = meta.single_end ? '' : '-pe'
    def memory     = task.memory.toGiga() + "G"

    def strandedness = 'non-strand-specific'
    if (meta.strandedness == 'forward') {
        strandedness = 'strand-specific-forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'strand-specific-reverse'
    }
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        rnaseq \\
        -bam $bam \\
        -gtf $gtf \\
        -p $strandedness \\
        $paired_end \\
        -outdir $prefix
    """
}