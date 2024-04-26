process CAT_FASTQ {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_low'

    input:
    tuple val(meta), path(read1_list), path(read2_list)

    output:
    tuple val(meta), path("${meta.id}*_1.merged.fq.gz"), emit: read1
    tuple val(meta), path("${meta.id}*_2.merged.fq.gz"), emit: read2

    script:
    def prefix   = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    def read1 = read1_list.collect{ it.toString() }
    def read2 = read2_list.collect{ it.toString() }
    if(read1.size == 1 && read2.size == 1){
        """
        ln -s ${read1.sort().join(' ')} ${prefix}_1.merged.fq.gz
        ln -s ${read2.sort().join(' ')} ${prefix}_2.merged.fq.gz
        """
    }else{
        """
        cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fq.gz
        cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fq.gz
        """
    }
}
