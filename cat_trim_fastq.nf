process CAT_TRIM_FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/cutqc/${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
        if(filename=~/versions.yml/){
            return null
        }else if(!params.save_merged_fastq && filename=~/\.(fq|fastq)\.gz/){
            return null // don't publish fastq file
        }else{
            return filename
        }
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1.merged.trimmed.fq.gz"), emit: read1
    tuple val(meta), path("*_2.merged.trimmed.fq.gz"), emit: read2
    tuple val(meta), path("*_cutqc_report.html"), emit: cutqc_report

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
            
            cutqc.sh ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz \\
            ${prefix}_cutqc_report.html \\
            $baseDir/bin/fastqc_report.Rmd \\
            -m $params.trimLength \\
            -j $task.cpus \\
            -q 0 \\
            -Q 30,30 \\
            $params.trimOpt
            
            ## remove merged fq.gz
            rm ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz
            
            """
        }else{
            """
            cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fq.gz
            cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fq.gz
            
            cutqc.sh ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz \\
            ${prefix}_cutqc_report.html \\
            $baseDir/bin/fastqc_report.Rmd \\
            -m $params.trimLength \\
            -j $task.cpus \\
            -q 0 \\
            -Q 30,30 \\
            $params.trimOpt
            
            ## remove merged fq.gz
            rm ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz
            
            """
        }
    }else{
        exit 1, 'Please provide both the read1 and the read2'
    }
}
