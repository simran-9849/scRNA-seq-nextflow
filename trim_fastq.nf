process TRIM_FASTQ {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}/cutqc",
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
    tuple val(meta), path(bc_read), path(cDNA_read)

    output:
    tuple val(meta), path("*_bc.trimmed.fq.gz"),   emit: bc_read
    tuple val(meta), path("*_cDNA.trimmed.fq.gz"), emit: cDNA_read
    tuple val(meta), path("${meta.id}.cutadapt.json"), emit: report_JSON
    
    script:
    def prefix   = "${meta.id}"
    """
    cutadapt -m $params.trimLength \\
             -j $task.cpus \\
             -q 0 \\
             -Q 30,30 \\
             --json ${prefix}.cutadapt.json \\
             -o ${prefix}_bc.trimmed.fq.gz -p ${prefix}_cDNA.trimmed.fq.gz \\
             $params.trimOpt \\
             ${bc_read} ${cDNA_read}
    """
}
