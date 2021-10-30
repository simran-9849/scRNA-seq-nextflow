// Code modified from nf-core ranseq pipeline
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './modules/nf-core_rnaseq/functions'

params.options = [:]
options        = initOptions(params.options)

process STARSOLO {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    path whitelist

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*d.out.bam.bai')   , emit: bai
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    tuple val(meta), path('*Solo.out/Gene')   , emit: solo_out
    path "versions.yml"                       , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab


    script:
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    //def barcodeMate = params.bc_read == "fastq_1" ? 1 : 2
    // Since starsolo default use single-end mode, activate this label for qualimap option
    meta.single_end = true
    // Assign read1 read2 to different list, code from cat_fastq
    def readList = reads.collect{ it.toString() }
    def read1 = []
    def read2 = []
    def cDNA_read = ""
    def bc_read = ""
    if (readList.size >= 2) {
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
    }else{
        exit 1, 'Please provide both the read1 and the read2'
    }

    if ( params.bc_read == "fastq_1" ){
        cDNA_read = read2.sort().join(",")
        bc_read = read1.sort().join(",")
    }else{
        cDNA_read = read1.sort().join(",")
        bc_read = read2.sort().join(",")
    }
    """
    STAR --runThreadN $task.cpus \\
         --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --soloBarcodeMate 0 \\
         --readFilesIn $cDNA_read $bc_read \\
         --outFileNamePrefix ${prefix}. \\
         --soloCBstart $params.soloCBstart \\
         --soloCBlen $params.soloCBlen \\
         --soloUMIstart $params.soloUMIstart \\
         --soloUMIlen $params.soloUMIlen \\
         --soloCBwhitelist $params.whitelist \\
         --soloBarcodeReadLength 0 \\
         --readFilesCommand zcat \\
         --soloType $params.soloType \\
         --clipAdapterType $params.clipAdapterType \\
         --outFilterScoreMin $params.outFilterScoreMin \\
         --soloCBmatchWLtype $params.soloCBmatchWLtype \\
         --soloUMIfiltering $params.soloUMIfiltering \\
         --soloUMIdedup $params.soloUMIdedup \\
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
         --outSAMtype BAM SortedByCoordinate \\

    samtools index *.bam

    ## ntfs fuseblk, gzip operation permission issue, use gzip -c instead
    for i in ${prefix}.Solo.out/Gene/raw/*; do gzip -c \$i > \${i}.gz; done
    for i in ${prefix}.Solo.out/Gene/filtered/*; do gzip -c \$i > \${i}.gz; done

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
    ${getSoftwareName(task.process)}: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
