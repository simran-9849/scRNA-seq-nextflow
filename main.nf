#!/usr/bin/env nextflow

// A lot of codes were directly inherited from https://github.com/nf-core/rnaseq

nextflow.enable.dsl=2

// NF-CORE MODULES

include { CAT_FASTQ } from './cat_fastq' addParams( options: ['publish_files': false] )
include { STARSOLO } from "./starsolo"
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './modules/nf-core_rnaseq/functions'
include { QUALIMAP_RNASEQ } from './modules/nf-core/modules/qualimap/rnaseq/main'
include { REPORT } from "./report"

// check mandatory params
if (!params.input) { exit 1, 'Input samplesheet not specified!' }
if (!params.genomeDir) { exit 1, 'Genome index DIR not specified!' }
if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.bc_read}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.cDNA_read}"
    }

    array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]

    return array
}


workflow {
    Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map{ create_fastq_channel(it) }
    .map {
        meta, fastq ->
            // meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .map {
        meta, fastq ->
            return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    // MODULE: Concatenate FastQ files from same sample if required

    ch_read1 = Channel.empty()
    ch_read2 = Channel.empty()

    CAT_FASTQ (
        ch_fastq
    )
    ch_read1 = CAT_FASTQ.out.read1
    ch_read2 = CAT_FASTQ.out.read2

    ch_genomeDir = Channel.fromPath(params.genomeDir)
    ch_genomeGTF = Channel.fromPath(params.genomeGTF)
    ch_whitelist = Channel.fromPath(params.whitelist)

    // if(params.trimReads){}

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_starsolo_out               = Channel.empty()
    ch_star_multiqc               = Channel.empty()

    if(params.bc_read == "fastq_1"){
        STARSOLO(
            ch_read2,
            ch_read1,
            ch_genomeDir,
            ch_genomeGTF,
            ch_whitelist
        )
    }else{
        STARSOLO(
            ch_read1,
            ch_read2,
            ch_genomeDir,
            ch_genomeGTF,
            ch_whitelist
        )
    }

    ch_genome_bam                 = STARSOLO.out.bam
    ch_genome_bam_index           = STARSOLO.out.bai
    ch_starsolo_out               = STARSOLO.out.solo_out
    ch_qualimap_multiqc           = Channel.empty()
    QUALIMAP_RNASEQ(
        ch_genome_bam,
        ch_genomeGTF
    )
    ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results

    REPORT(
        ch_starsolo_out,
        ch_qualimap_multiqc
    )
}