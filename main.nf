#!/usr/bin/env nextflow

// A lot of codes were directly inherited from https://github.com/nf-core/rnaseq

nextflow.enable.dsl=2

// NF-CORE MODULES

include { CAT_FASTQ } from './cat_fastq'
include { TRIM_FASTQ } from './trim_fastq'
include { MULTIQC } from './multiqc'
include { STARSOLO; STAR_MKREF; } from "./starsolo"
include { GENECOVERAGE; FEATURESTATS } from "./gene_coverage"
include { CHECK_SATURATION } from "./sequencing_saturation"
include { GET_VERSIONS } from "./present_version"
include { REPORT } from "./report"

def entryWorkflow() {
    int idx = -1
    cmd = workflow.commandLine.split(" ")
    if(cmd.contains("-entry")){
        idx = Arrays.asList(cmd).indexOf("-entry") + 1
        return cmd[idx]
    }else{
        return "scRNAseq"
    }
}

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
    
    if(!row.expected_cells){
        exit 1, "ERROR: Please check input samplesheet -> Please specify expected_cells column"
    }
    meta.expected_cells = row.expected_cells

    if(entryWorkflow() == "vdj"){
        if(!row.feature_types){
            exit 1, "ERROR: Please check input samplesheet -> feature_types column does not exist!\n"
        }else{
            meta.feature_types = row.feature_types
        }
    }

    array = [ meta, file(row.fastq_1), file(row.fastq_2) ]

    return array
}


workflow {
    // default workflow scRNA-seq
    scRNAseq()
}

workflow scRNAseq {
    // check mandatory params
    if (!params.input) { exit 1, 'Input samplesheet not specified!' }
    if (!params.genomeDir) { exit 1, 'Genome index DIR not specified!' }
    if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
    // only support unique-gene reads for now
    if (params.soloMultiMappers != "Unique") { exit 1, 'Only support unique-gene reads for now, please use soloMultiMappers = "Unique" instead!' }

    Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map{ create_fastq_channel(it) }
    .groupTuple(by: [0])
    .set { ch_fastq }

    // MODULE: Concatenate FastQ files from the same sample if required
    CAT_FASTQ( ch_fastq )

    ch_bc_read = Channel.empty()
    ch_cDNA_read = Channel.empty()

    if ( params.bc_read == "fastq_1" ){
        ch_bc_read = CAT_FASTQ.out.read1
        ch_cDNA_read = CAT_FASTQ.out.read2
    }else{
        ch_bc_read = CAT_FASTQ.out.read2
        ch_cDNA_read = CAT_FASTQ.out.read1
    }
    ch_merged_fastq = ch_bc_read.join(ch_cDNA_read, by:[0])

    TRIM_FASTQ( ch_merged_fastq )
    CAT_FASTQ.out.read1
    .join(CAT_FASTQ.out.read2, by:[0])
    .join(TRIM_FASTQ.out.bc_read, by:[0])
    .join(TRIM_FASTQ.out.cDNA_read, by:[0])
    .join(TRIM_FASTQ.out.report_JSON, by:[0])
    .set { ch_multiqc_input }

    MULTIQC( ch_multiqc_input )

    ch_genomeDir = file(params.genomeDir)
    ch_genomeGTF = file(params.genomeGTF)
    
    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_starsolo_out               = Channel.empty()
    ch_star_multiqc               = Channel.empty()

    ch_whitelist = Channel.fromPath(params.whitelist.split(" ").toList())

    // concatenate reads as single input channel
    // after concatenation, the first one is bc read, the second is cDNA read
    ch_starsolo_inputReads = TRIM_FASTQ.out.bc_read.join(TRIM_FASTQ.out.cDNA_read, by:[0])
    STARSOLO(
        ch_starsolo_inputReads,
        ch_genomeDir,
        ch_genomeGTF,
        ch_whitelist.toList()
    )
    ch_genome_bam       = STARSOLO.out.bam
    ch_genome_bam_index = STARSOLO.out.bai
    ch_filteredDir      = STARSOLO.out.filteredDir
    ch_rawDir           = STARSOLO.out.rawDir
    ch_starsolo_summary = STARSOLO.out.summary_unique
    ch_starsolo_UMI     = STARSOLO.out.UMI_file_unique

    ch_saturation_input = ch_genome_bam.join(ch_filteredDir, by:[0])
    CHECK_SATURATION(
        ch_saturation_input,
        ch_whitelist.toList()
    )
    GET_VERSIONS()

    ch_featureStats   = Channel.empty()
    ch_geneCoverage   = Channel.empty()

    ch_geneCoverage_input = ch_genome_bam.join(ch_genome_bam_index, by:[0])
    FEATURESTATS(
        ch_geneCoverage_input,
        ch_genomeGTF
    )
    
    GENECOVERAGE(
        ch_geneCoverage_input,
        ch_genomeGTF
    )

    ch_featureStats  = FEATURESTATS.out.stats
    ch_geneCoverage  = GENECOVERAGE.out.matrix


    ch_starsolo_summary
    .join(ch_starsolo_UMI, by:[0])
    .join(ch_rawDir, by:[0])
    .join(ch_filteredDir, by:[0])
    .join(ch_featureStats, by:[0])
    .join(ch_geneCoverage, by:[0])
    .join(CHECK_SATURATION.out.outJSON, by:[0])
    .set{ ch_report_input }

    REPORT(
        ch_report_input,
        GET_VERSIONS.out.json
    )

}

workflow mkref {
    // check mandatory params
    if (!params.genomeFasta) { exit 1, 'Genome Fasta not specified!' }
    if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
    if (!params.refoutDir) { exit 1, 'Reference output directory not specified!'}

    ch_genomeFasta = file(params.genomeFasta)
    ch_genomeGTF = file(params.genomeGTF)
    STAR_MKREF(
        ch_genomeFasta,
        ch_genomeGTF
    )
}

// sub-workflow for vdj analysis

include { VDJ_CELLCALLING_WITHOUTGEX; VDJ_CELLCALLING_WITHGEX; VDJ_ASSEMBLY; VDJ_METRICS } from "./vdj/trust4_vdj"
include { GET_VERSIONS_VDJ } from "./vdj/present_version_vdj"
include { REPORT_VDJ } from "./vdj/report_vdj"

workflow vdj {
    vdj_process()
    vdj_report(
        vdj_process.out.starsolo_summary,
        vdj_process.out.starsolo_umi,
        vdj_process.out.starsolo_filteredDir,
        vdj_process.out.featureStats,
        vdj_process.out.geneCoverage,
        vdj_process.out.saturation_json,
        vdj_process.out.trust4_report,
        vdj_process.out.trust4_airr,
        vdj_process.out.trust4_readsAssign,
        vdj_process.out.trust4_barcode,
        vdj_process.out.trust4_finalOut,
        vdj_process.out.trust4_kneeOut,
        vdj_process.out.trust4_cellOut,
        vdj_process.out.trust4_metrics,
        vdj_process.out.trust4_cloneType
    )
}

workflow vdj_process {    
    main:
    // check mandatory params
    if (!params.input) { exit 1, 'Input samplesheet not specified!' }
    if (!params.genomeDir) { exit 1, 'Genome index DIR not specified!' }
    if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
    if (!params.trust4_vdj_refGenome_fasta) { exit 1, 'TRUST4 refGenome fasta not specified!' }
    if (!params.trust4_vdj_imgt_fasta) { exit 1, 'TRUST4 imgt fasta not specified!' }
    
    // use different sampleList for vdj pipeline
    // one more column: feature_types
    Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map{ create_fastq_channel(it) }
    .groupTuple(by: [0])
    .set { ch_fastq }

    // MODULE: Concatenate FastQ files from the same sample if required
    CAT_FASTQ( ch_fastq )

    // process vdj first
    ch_bc_read = Channel.empty()
    ch_cDNA_read = Channel.empty()

    if ( params.bc_read == "fastq_1" ){
        ch_bc_read = CAT_FASTQ.out.read1
        ch_cDNA_read = CAT_FASTQ.out.read2
    }else{
        ch_bc_read = CAT_FASTQ.out.read2
        ch_cDNA_read = CAT_FASTQ.out.read1
    }

    ch_merged_fastq = ch_bc_read.join(ch_cDNA_read, by:[0])

    TRIM_FASTQ( ch_merged_fastq )

    CAT_FASTQ.out.read1
    .join(CAT_FASTQ.out.read2, by:[0])
    .join(TRIM_FASTQ.out.bc_read, by:[0])
    .join(TRIM_FASTQ.out.cDNA_read, by:[0])
    .join(TRIM_FASTQ.out.report_JSON, by:[0])
    .set { ch_multiqc_input }

    MULTIQC( ch_multiqc_input )
    
    ch_genomeDir = file(params.genomeDir)
    ch_genomeGTF = file(params.genomeGTF)
    ch_vdj_refGenome_fasta = file(params.trust4_vdj_refGenome_fasta)
    ch_vdj_imgt_fasta = file(params.trust4_vdj_imgt_fasta)

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_starsolo_summary           = Channel.empty()
    ch_starsolo_umi               = Channel.empty()
    ch_starsolo_filteredDir       = Channel.empty()
    ch_starsolo_rawDir            = Channel.empty()
    ch_saturation_json            = Channel.empty()
    // force using parameters for 5'-RNAseq
    //params.soloStrand = "Reverse"

    ch_whitelist = Channel.fromPath(params.whitelist.split(" ").toList())
    // concatenate reads as single input channel
    // after concatenation, the first one is bc read, the second is cDNA read
    ch_starsolo_inputReads = TRIM_FASTQ.out.bc_read.join(TRIM_FASTQ.out.cDNA_read, by:[0])

    STARSOLO(
        ch_starsolo_inputReads,
        ch_genomeDir,
        ch_genomeGTF,
        ch_whitelist.toList()
    )
    ch_genome_bam           = STARSOLO.out.bam
    ch_genome_bam_index     = STARSOLO.out.bai
    ch_starsolo_filteredDir = STARSOLO.out.filteredDir
    ch_starsolo_rawDir      = STARSOLO.out.rawDir
    ch_starsolo_summary     = STARSOLO.out.summary_unique
    ch_starsolo_umi         = STARSOLO.out.UMI_file_unique
    
    ch_genome_bam
    .join(ch_starsolo_filteredDir, by:[0])
    .filter{ it[0].feature_types == "GEX" }
    .set{ ch_saturation_input }

    CHECK_SATURATION(
        ch_saturation_input,
        ch_whitelist.toList()
    )
    ch_saturation_json = CHECK_SATURATION.out.outJSON

    
    ch_featureStats   = Channel.empty()
    ch_geneCoverage   = Channel.empty()

    ch_genome_bam
        .join(ch_genome_bam_index, by:[0])
        .filter{ it[0].feature_types == "GEX" }
        .set{ ch_geneCoverage_input }

    FEATURESTATS(
        ch_geneCoverage_input,
        ch_genomeGTF
    )
    
    GENECOVERAGE(
        ch_geneCoverage_input,
        ch_genomeGTF
    )

    ch_featureStats  = FEATURESTATS.out.stats
    ch_geneCoverage  = GENECOVERAGE.out.matrix
    
    // aggregate BAM files to check if GEX library exists
    TRIM_FASTQ.out.cDNA_read
    .join(TRIM_FASTQ.out.bc_read, by:[0])
    .join(ch_genome_bam, by:[0])
    .set{ ch_cellCalling_input }

    // 
    ch_cellCalling_input
    .map{
        meta, cDNAread, bcRead, bam ->
        tuple(meta.id, meta.expected_cells, meta.feature_types)
    }
    .groupTuple(by: [0])
    .filter{ !it[2].contains("GEX") }
    .transpose()
    .filter{ it[2] != "GEX" }
    .map{ id, expected_cells, feature_types -> tuple([id: id, expected_cells: expected_cells, feature_types: feature_types])}
    .join(ch_cellCalling_input, by:[0])
    .set{ ch_cellCalling_withoutGEX }

    ch_cellCalling_input
    .map{
        meta, cDNAread, bcRead, bam ->
        tuple(meta.id, meta.expected_cells, meta.feature_types)
    }
    .groupTuple(by: [0])
    .filter{ it[2].contains("GEX") }
    .transpose()
    .filter{ it[2] != "GEX" }
    .map{ id, expected_cells, feature_types -> tuple([id: id, expected_cells: expected_cells, feature_types: feature_types])}
    .join(ch_cellCalling_input, by:[0])
    .set{ ch_cellCalling_withGEX_noFilteredDir }

    ch_starsolo_filteredDir
    .filter{ it[0].feature_types == "GEX" }
    .map{ meta, filteredDir -> tuple(meta.id, filteredDir)}
    .set{ ch_cellCalling_withGEX_barcodes }

    ch_cellCalling_withGEX_noFilteredDir
    .map{
        meta, cDNAread, bcRead, bam ->
        tuple(meta.id, meta.expected_cells, meta.feature_types, cDNAread, bcRead, bam)
    }
    .combine(ch_cellCalling_withGEX_barcodes, by:[0])
    .map{
        id, expected_cells, feature_types, cDNAread, bcRead, bam, filteredDir ->
        tuple([id: id, expected_cells: expected_cells, feature_types: feature_types], cDNAread, bcRead, bam, filteredDir)
    }
    .set{ ch_cellCalling_withGEX }
    
    VDJ_CELLCALLING_WITHOUTGEX(
        ch_cellCalling_withoutGEX
    )

    VDJ_CELLCALLING_WITHGEX(
        ch_cellCalling_withGEX
    )

    ch_trust4_cDNAread = VDJ_CELLCALLING_WITHOUTGEX.out.cDNAread
        .concat(VDJ_CELLCALLING_WITHGEX.out.cDNAread)
    ch_trust4_bcRead = VDJ_CELLCALLING_WITHOUTGEX.out.bcRead
        .concat(VDJ_CELLCALLING_WITHGEX.out.bcRead)
    ch_trust4_kneeOut = VDJ_CELLCALLING_WITHOUTGEX.out.kneeOut
        .concat(VDJ_CELLCALLING_WITHGEX.out.kneeOut)
    ch_trust4_raw_cellOut = VDJ_CELLCALLING_WITHOUTGEX.out.cellOut
        .concat(VDJ_CELLCALLING_WITHGEX.out.cellOut)
    ch_trust4_readList = VDJ_CELLCALLING_WITHOUTGEX.out.readList
        .concat(VDJ_CELLCALLING_WITHGEX.out.readList)
    ch_trust4_barcode = VDJ_CELLCALLING_WITHOUTGEX.out.barcode
        .concat(VDJ_CELLCALLING_WITHGEX.out.barcode)
    ch_trust4_umi = VDJ_CELLCALLING_WITHOUTGEX.out.umi
        .concat(VDJ_CELLCALLING_WITHGEX.out.umi)

    ch_trust4_cDNAread
    .join(ch_trust4_bcRead, by:[0])
    .join(ch_trust4_readList, by:[0])
    .join(ch_trust4_barcode, by:[0])
    .join(ch_trust4_umi, by:[0])
    .join(ch_trust4_kneeOut, by:[0])
    .join(ch_trust4_raw_cellOut, by:[0])
    .set{ ch_vdj_assembly_input }

    VDJ_ASSEMBLY(
        ch_vdj_assembly_input,
        ch_vdj_refGenome_fasta,
        ch_vdj_imgt_fasta
    )

    ch_trust4_report        = VDJ_ASSEMBLY.out.report
    ch_trust4_airr          = VDJ_ASSEMBLY.out.airr
    ch_trust4_readsAssign   = VDJ_ASSEMBLY.out.readsAssign
    ch_trust4_finalOut      = VDJ_ASSEMBLY.out.finalOut

    ch_trust4_report
    .join(ch_trust4_airr, by:[0])
    .join(ch_trust4_readsAssign, by:[0])
    .join(ch_trust4_barcode, by:[0])
    .join(ch_trust4_kneeOut, by:[0])
    .join(ch_starsolo_summary.filter{ it[0].feature_types != "GEX" }, by:[0])
    .set{ ch_vdj_metrics_input }

    VDJ_METRICS(
        ch_vdj_metrics_input
    )

    ch_trust4_cellOut       = VDJ_METRICS.out.cellOut
    ch_trust4_metrics       = VDJ_METRICS.out.metricsJSON
    ch_trust4_cloneType     = VDJ_METRICS.out.cloneType


    emit:
        starsolo_summary     = ch_starsolo_summary
        starsolo_bam         = ch_genome_bam
        starsolo_umi         = ch_starsolo_umi
        starsolo_filteredDir = ch_starsolo_filteredDir
        featureStats         = ch_featureStats
        geneCoverage         = ch_geneCoverage
        saturation_json      = ch_saturation_json
        trust4_report        = ch_trust4_report
        trust4_airr          = ch_trust4_airr
        trust4_readsAssign   = ch_trust4_readsAssign
        trust4_barcode       = ch_trust4_barcode
        trust4_finalOut      = ch_trust4_finalOut
        trust4_kneeOut       = ch_trust4_kneeOut
        trust4_cellOut       = ch_trust4_cellOut
        trust4_metrics       = ch_trust4_metrics
        trust4_cloneType     = ch_trust4_cloneType
}

def collapse_vdj_ch ( ch_input ) {
    //ch_input.view()
    ch_expanded = ch_input
        .map{
            meta, file ->
            tuple(meta.id, meta.feature_types, file)
        }
        .groupTuple(by:[0], sort:true)
        .map{
            id, feature_types, file ->
            // feature_types: ["VDJ-B", "VDJ-T"]
            // file: [BCR_file, TCR_file]
            def bcr_file = file.any{ it.getName() =~ /VDJ-B/ } ? file.find{ it.getName() =~/VDJ-B/ } : ""
            def tcr_file = file.any{ it.getName() =~ /VDJ-T/ }? file.find{ it.getName() =~ /VDJ-T/ } : ""
            tuple([id: id], [bcr_file, tcr_file])
        }

    return ch_expanded
}

def select_gex_ch ( ch_input ) {
    ch_selected = ch_input
        .map {
            meta, file ->
            if( meta.feature_types == "GEX" )
                tuple([id: meta.id], file)
        }
        .groupTuple(by:[0])

    return ch_selected
}


workflow vdj_report {
    take:
    starsolo_summary
    starsolo_umi
    starsolo_filteredDir
    featureStats
    geneCoverage
    saturation_json
    trust4_report
    trust4_airr
    trust4_readsAssign
    trust4_barcode
    trust4_finalOut
    trust4_kneeOut
    trust4_cellOut
    trust4_metrics
    trust4_cloneType

    main:

    ch_vdj_metrics = collapse_vdj_ch(trust4_metrics)
    ch_vdj_cloneType = collapse_vdj_ch(trust4_cloneType)
    ch_vdj_cellOut = collapse_vdj_ch(trust4_cellOut)
    ch_vdj_report = collapse_vdj_ch(trust4_report)
    ch_vdj_airr = collapse_vdj_ch(trust4_airr)
    ch_vdj_kneeOut = collapse_vdj_ch(trust4_kneeOut)
    ch_vdj_finalOut = collapse_vdj_ch(trust4_finalOut)

    ch_starsolo_summary_gex = select_gex_ch(starsolo_summary)
    ch_featureStats_gex = select_gex_ch(featureStats)
    ch_geneCoverage_gex = select_gex_ch(geneCoverage)
    ch_starsolo_umi_gex = select_gex_ch(starsolo_umi)
    ch_starsolo_filteredDir_gex = select_gex_ch(starsolo_filteredDir)
    ch_saturation_json_gex = select_gex_ch(saturation_json)

    //ch_vdj_report contains all samples, including ones without GEX library
    ch_vdj_report
        .map{
            meta, file ->
            tuple(meta)
        }
        .join(ch_starsolo_summary_gex, remainder: true, by:[0])
        .join(ch_starsolo_umi_gex, remainder: true, by:[0])
        .join(ch_starsolo_filteredDir_gex, remainder: true, by:[0])
        .join(ch_featureStats_gex, remainder: true, by:[0])
        .join(ch_geneCoverage_gex, remainder: true, by:[0])
        .join(ch_saturation_json_gex, remainder: true, by:[0])
        .join(ch_vdj_report, remainder: true, by:[0])
        .join(ch_vdj_airr, remainder: true, by:[0])
        .join(ch_vdj_kneeOut, remainder: true, by:[0])
        .join(ch_vdj_finalOut, remainder: true, by:[0])
        .join(ch_vdj_cellOut, remainder: true, by:[0])
        .join(ch_vdj_metrics, remainder: true, by:[0])
        .join(ch_vdj_cloneType, remainder: true, by:[0])
        .set{ ch_report_input }

    GET_VERSIONS_VDJ()

    REPORT_VDJ(
        ch_report_input,
        GET_VERSIONS_VDJ.out.json
    )
}
