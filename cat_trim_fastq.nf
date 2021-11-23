// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName; getPathFromList } from './modules/nf-core_rnaseq/functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_TRIM_FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
        if(filename=~/versions.yml/){
            return null
        }else if(!params.save_merged_fastq && filename=~/\.(fq|fastq)\.gz/){
            return null // don't publish fastq file
        }else{
            return "${getPathFromList(['cutqc'])}/$filename"
        }
    }

    // conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    // } else {
    //     container "biocontainers/biocontainers:v1.2.0_cv1"
    // }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_1.merged.trimmed.fq.gz"), emit: read1
    tuple val(meta), path("*_2.merged.trimmed.fq.gz"), emit: read2
    tuple val(meta), path("*_cutqc_report.html"), emit: cutqc_report
    path "versions.yml"                       , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def readList = reads.collect{ it.toString() }
    if (readList.size >= 2) {
        def read1 = []
        def read2 = []
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
        """
        cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fq.gz
        cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fq.gz

        cutqc.sh ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz \\
        ${prefix}_cutqc_report.html \\
        $baseDir/bin/fastqc_report.Rmd \\
        -m $params.trimLength \\
        -j $task.cpus \\
        -q 30 \\
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \\
        -g TCTTTCCCTACACGACGCTCTTCCGATCT \\
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGA  \\
        -A "T{30}" \\
        -G GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \\
        -G AAGCAGTGGTATCAACGCAGAGTACATGGG \\

        ## remove merged fq.gz
        rm ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
        END_VERSIONS
        """
    }else{
        exit 1, 'Please provide both the read1 and the read2'
    }
}
