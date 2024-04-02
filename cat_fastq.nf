process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1.merged.fq.gz"), emit: read1
    tuple val(meta), path("*_2.merged.fq.gz"), emit: read2

    script:
    def prefix   = "${meta.id}"
    def readList = reads.collect{ it.toString() }
    if (readList.size >= 2) {
        def read1 = []
        def read2 = []
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
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
    }else{
        exit 1, 'Please provide both the read1 and the read2'
    }
}
