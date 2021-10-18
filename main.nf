#!/usr/bin/env nextflow

// A lot of codes were directly inherited from https://github.com/nf-core/rnaseq

nextflow.enable.dsl=2

// NF-CORE MODULES

def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

include { CAT_FASTQ } from './modules/nf-core/modules/cat/fastq/main' addParams( options: cat_fastq_options )
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './modules/nf-core_rnaseq/functions'
include { QUALIMAP_RNASEQ } from './modules/nf-core/modules/qualimap/rnaseq/main' addParams( options: modules['qualimap_rnaseq'] )

// check mandatory params
if (!params.input) { exit 1, 'Input samplesheet not specified!' }
if (!params.genomeDir) { exit 1, 'Genome index DIR not specified!' }
if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    if (!file(row.bc_read).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Barcode read FastQ file does not exist!\n${row.bc_read}"
    }

    if (!file(row.cDNA_read).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> cDNA read FastQ file does not exist!\n${row.cDNA_read}"
    }
    // Ensure bc_read is the first fastq
    array = [ meta, [ file(row.bc_read), file(row.cDNA_read) ] ]

    return array
}

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
    """
    STAR -- runThreadN $task.cpus \\
         --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $reads \\
         --outFileNamePrefix ${prefix}. \\
         --soloCBstart $params.soloCBstart \\
         --soloCBlen $params.soloCBlen \\
         --soloUMIstart $params.soloUMIstart \\
         --soloUMIlen $params.soloUMIlen \\
         --soloCBwhitelist $params.whitelist \\
         --soloBarcodeReadLength $params.soloBarcodeReadLength \\
         --readFilesCommand zcat \\
         --soloType $params.soloType \\
         --clipAdapterType $params.clipAdapterType \\
         --outFilterScoreMin $params.outFilterScoreMin \\
         --soloCBmatchWLtype $params.soloCBmatchWLtype \\
         --soloUMIfiltering $params.soloUMIfiltering \\
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
         --outSAMtype BAM SortedByCoordinate \\

    samtools index *.bam
    """
}

workflow {
    Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map{ create_fastq_channel(it) }
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    // MODULE: Concatenate FastQ files from same sample if required

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    ch_genomeDir = Channel.fromPath(params.genomeDir)
    ch_genomeGTF = Channel.fromPath(params.genomeGTF)
    ch_whitelist = Channel.fromPath(params.whitelist)

    if(params.trimReads){}

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_star_multiqc               = Channel.empty()

    STARSOLO(
        ch_cat_fastq
        ch_genomeDir
        ch_genomeGTF
        path(params.whitelist)
    )

    ch_genome_bam                 = STARSOLO.out.bam
    ch_genome_bam_index           = STARSOLO.out.bai

    ch_qualimap_multiqc           = Channel.empty()
    QUALIMAP_RNASEQ(
        ch_genome_bam,
        ch_genomeGTF
    )
    ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results
}