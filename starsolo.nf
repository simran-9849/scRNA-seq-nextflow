process STARSOLO {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/starsolo/${meta.id}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/sortedByCoord.out.bam/){
            return null
        }else{
            return filename
        }
    }
    input:
    tuple val(meta), path(cDNA_read)
    tuple val(meta), path(bc_read)
    path index
    path gtf
    path whitelist

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*d.out.bam.bai')   , emit: bai
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    tuple val(meta), path('*Solo.out/Gene*')   , emit: solo_out

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab


    script:
    def prefix     = "${meta.id}"
    //def barcodeMate = params.bc_read == "fastq_1" ? 1 : 2
    // Since starsolo default use single-end mode, activate this label for qualimap option
    meta.single_end = true

    """
    STAR --runThreadN $task.cpus \\
    --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --soloBarcodeMate 0 \\
    --readFilesIn $cDNA_read $bc_read \\
    --outFileNamePrefix ${prefix}. \\
    --soloStrand $params.soloStrand \\
    --soloCBstart $params.soloCBstart \\
    --soloCBlen $params.soloCBlen \\
    --soloUMIstart $params.soloUMIstart \\
    --soloUMIlen $params.soloUMIlen \\
    --soloCBwhitelist $whitelist \\
    --soloBarcodeReadLength 0 \\
    --readFilesCommand zcat \\
    --soloType $params.soloType \\
    --clipAdapterType $params.clipAdapterType \\
    --outFilterScoreMin $params.outFilterScoreMin \\
    --soloCBmatchWLtype $params.soloCBmatchWLtype \\
    --soloUMIfiltering $params.soloUMIfiltering \\
    --soloUMIdedup $params.soloUMIdedup \\
    --soloMultiMappers $params.soloMultiMappers \\
    --soloFeatures $params.soloFeatures \\
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
    --outSAMtype BAM SortedByCoordinate

    samtools index ${prefix}.Aligned.sortedByCoord.out.bam

    pigz -p $task.cpus ${prefix}.Solo.out/*/raw/*
    pigz -p $task.cpus ${prefix}.Solo.out/*/filtered/*
    """
}

process STARSOLO_COMPLEX {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/starsolo/${meta.id}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(cDNA_read)
    tuple val(meta), path(bc_read)
    path index
    path gtf
    path whitelist
    path whitelist2

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*d.out.bam.bai')   , emit: bai
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    tuple val(meta), path('*Solo.out/Gene')   , emit: solo_out
    //path "versions.yml"                       , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    script:
    def prefix     = "${meta.id}"
    //def barcodeMate = params.bc_read == "fastq_1" ? 1 : 2
    // Since starsolo default use single-end mode, activate this label for qualimap option
    meta.single_end = true
    """
    STAR --runThreadN $task.cpus \\
    --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --soloBarcodeMate 0 \\
    --readFilesIn $cDNA_read $bc_read \\
    --soloBarcodeReadLength 0 \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${prefix}. \\
    --soloStrand $params.soloStrand \\
    --soloType $params.soloType \\
    --soloCBposition $params.complexCBposition \\
    --soloUMIposition $params.complexUMIposition \\
    --soloAdapterSequence $params.complexAdapterSequence \\
    --soloAdapterMismatchesNmax $params.complexAdapterMismatchesNmax \\
    --soloCBwhitelist $whitelist $whitelist2 \\
    --clipAdapterType $params.clipAdapterType \\
    --outFilterScoreMin $params.outFilterScoreMin \\
    --soloCBmatchWLtype $params.soloCBmatchWLtype \\
    --soloUMIfiltering $params.soloUMIfiltering \\
    --soloUMIdedup $params.soloUMIdedup \\
    --soloMultiMappers $params.soloMultiMappers \\
    --soloFeatures $params.soloFeatures \\
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
    --outSAMtype BAM SortedByCoordinate

    samtools index *.bam

    pigz -p $task.cpus ${prefix}.Solo.out/Gene/raw/*
    pigz -p $task.cpus ${prefix}.Solo.out/Gene/filtered/*
    """
}

process STAR_MKREF {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: "move",
        enabled: params.outdir as boolean
    input:
    path genomeFasta
    path genomeGTF
    output:
    path("${params.refoutDir}")       , emit: refPath
    script:
    """
    STAR --runMode genomeGenerate \\
    --runThreadN $task.cpus \\
    --genomeDir $params.refoutDir \\
    --genomeFastaFiles $genomeFasta \\
    --sjdbGTFfile $genomeGTF $params.mkrefOpt
    """
}
