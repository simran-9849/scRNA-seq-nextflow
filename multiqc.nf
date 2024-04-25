process MULTIQC {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}/multiqc",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/report.html/){
            return filename
        }else{
            return null
        }
    }

    input:
    tuple val(meta),
          path(merged_read1),
          path(merged_read2),
          path(trimmed_bc_read),
          path(trimmed_cDNA_read),
          path(cutadapt_report)

    output:
    tuple val(meta), path("${meta.id}_multiqc_report.html"), emit: report    

    script:
    def prefix = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    def fastqcThreads = Math.min(4, task.cpus)
    """
    ## rename input files
    ln -s ${merged_read1} ${prefix}_raw_R1.fq.gz
    ln -s ${merged_read2} ${prefix}_raw_R2.fq.gz
    ln -s ${trimmed_bc_read} ${prefix}_bc_trimmed.fq.gz
    ln -s ${trimmed_cDNA_read} ${prefix}_cDNA_trimmed.fq.gz
    fastqc -t $fastqcThreads --nogroup ${prefix}_raw_R1.fq.gz ${prefix}_raw_R2.fq.gz ${prefix}_bc_trimmed.fq.gz ${prefix}_cDNA_trimmed.fq.gz
    multiqc --title ${prefix} .
    """
}