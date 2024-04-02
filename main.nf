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

    array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]

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
    CAT_FASTQ.out.rea1
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

    //if(params.soloMultiMappers != "Unique"){
    //    STARSOLO_MULTIPLE(
    //        STARSOLO.out.rawDir
    //    )
    //    CHECK_SATURATION(
    //        ch_genome_bam,
    //        STARSOLO_MULTIPLE.out.filteredDir,
    //        ch_whitelist.toList()
    //    )
    //    STARSOLO_MULT_SUMMARY(
    //        STARSOLO.out.cellReads_stats,
    //        STARSOLO_MULTIPLE.out.filteredDir,
    //        STARSOLO.out.summary_unique,
    //        CHECK_SATURATION.out.outJSON
    //    )
    //    STARSOLO_MULT_UMI(
    //        STARSOLO.out.cellReads_stats   
    //    )
    //    GET_VERSIONS(
    //        CHECK_SATURATION.out.outJSON
    //    )
    //}else{

        ch_saturation_input = ch_genome_bam.join(ch_filteredDir, by:[0])
        CHECK_SATURATION(
            ch_saturation_input,
            ch_whitelist.toList()
        )
        GET_VERSIONS(
            CHECK_SATURATION.out.outJSON
        )
    //}

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

    //if(params.soloMultiMappers != "Unique"){
    //    REPORT(
    //        STARSOLO_MULT_SUMMARY.out.summary_multiple,
    //        STARSOLO_MULT_UMI.out.UMI_file_multiple,
    //        STARSOLO_MULTIPLE.out.filteredDir,
    //        ch_featureStats,
    //        ch_geneCoverage,
    //        CHECK_SATURATION.out.outJSON,
    //        GET_VERSIONS.out.versions
    //    )
    //}else{
        ch_starsolo_summary
        .join(ch_starsolo_UMI, by:[0])
        .join(ch_rawDir, by:[0])
        .join(ch_filteredDir, by:[0])
        .join(ch_featureStats, by:[0])
        .join(ch_geneCoverage, by:[0])
        .join(CHECK_SATURATION.out.outJSON, by:[0])
        .join(GET_VERSIONS.out.versions, by:[0])
        .set{ ch_report_input }
        REPORT(
            ch_report_input
        )
    //}
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

include { CAT_TRIM_FASTQ_VDJ } from "./vdj/cat_trim_fastq_vdj"
include { STARSOLO_VDJ; STARSOLO_MULTIPLE_VDJ; STARSOLO_MULT_SUMMARY_VDJ; STARSOLO_MULT_UMI_VDJ; } from "./vdj/starsolo_vdj"
include { QUALIMAP_VDJ } from "./vdj/qualimap_vdj"
include { CHECK_SATURATION_VDJ } from "./vdj/sequencing_saturation_vdj.nf"
include { TRUST4_VDJ; VDJ_METRICS } from "./vdj/trust4_vdj"
include { GET_VERSIONS_VDJ } from "./vdj/present_version_vdj"
include { REPORT_VDJ } from "./vdj/report_vdj"

workflow vdj {
    vdj_process()
    vdj_report(
        vdj_process.out.starsolo_summary,
        vdj_process.out.starsolo_umi,
        vdj_process.out.starsolo_filteredDir,
        vdj_process.out.qualimap_outDir,
        vdj_process.out.saturation_json,
        vdj_process.out.trust4_report,
        vdj_process.out.trust4_airr,
        vdj_process.out.trust4_readsAssign,
        vdj_process.out.trust4_barcode,
        vdj_process.out.trust4_finalOut,
        vdj_process.out.trust4_kneeOut,
        vdj_process.out.trust4_cellOut
    )
}

workflow vdj_process {    
    main:
    // check mandatory params
    if (!params.input) { exit 1, 'Input samplesheet not specified!' }
    if (!params.genomeDir) { exit 1, 'Genome index DIR not specified!' }
    if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
    
    // use different sampleList for vdj pipeline
    // one more column: feature_types
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

    // process vdj first
    ch_bc_read = Channel.empty()
    ch_cDNA_read = Channel.empty()
    CAT_TRIM_FASTQ_VDJ( ch_fastq )
    if ( params.bc_read == "fastq_1" ){
        ch_bc_read = CAT_TRIM_FASTQ_VDJ.out.read1
        ch_cDNA_read = CAT_TRIM_FASTQ_VDJ.out.read2
    }else{
        ch_bc_read = CAT_TRIM_FASTQ_VDJ.out.read2
        ch_cDNA_read = CAT_TRIM_FASTQ_VDJ.out.read1
    }

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
    ch_qualimap_outDir            = Channel.empty()
    ch_saturation_json            = Channel.empty()
    // force using parameters for 5'-RNAseq
    //params.soloStrand = "Reverse"

    ch_whitelist = Channel.fromPath(params.whitelist.split(" ").toList())

    STARSOLO_VDJ(
        ch_cDNA_read,
        ch_bc_read,
        ch_genomeDir,
        ch_genomeGTF,
        ch_whitelist.toList(),
    )
    ch_genome_bam           = STARSOLO_VDJ.out.bam_sorted
    ch_genome_bam_index     = STARSOLO_VDJ.out.bai
    ch_starsolo_filteredDir = STARSOLO_VDJ.out.filteredDir
    ch_starsolo_summary     = STARSOLO_VDJ.out.summary_unique
    ch_starsolo_umi         = STARSOLO_VDJ.out.UMI_file_unique

    
    //if(params.soloMultiMappers != "Unique"){
    //    STARSOLO_MULTIPLE_VDJ(
    //        STARSOLO_VDJ.out.rawDir
    //    )
    //    CHECK_SATURATION_VDJ(
    //        ch_genome_bam,
    //        STARSOLO_MULTIPLE_VDJ.out.filteredDir,
    //        ch_whitelist
    //    )
    //    STARSOLO_MULT_SUMMARY_VDJ(
    //        STARSOLO_VDJ.out.cellReads_stats,
    //        STARSOLO_MULTIPLE_VDJ.out.filteredDir,
    //        STARSOLO_VDJ.out.summary_unique,
    //        CHECK_SATURATION_VDJ.out.outJSON
    //    )
    //    STARSOLO_MULT_UMI_VDJ(
    //        STARSOLO_VDJ.out.cellReads_stats   
    //    )
    //    //GET_VERSIONS(
    //    //    CHECK_SATURATION_VDJ.out.outJSON
    //    //)
    //}else{
    //    CHECK_SATURATION_VDJ(
    //        STARSOLO_VDJ.out.bam,
    //        STARSOLO_VDJ.out.filteredDir,
    //        ch_whitelist
    //    )
    //    //GET_VERSIONS(
    //    //    CHECK_SATURATION_VDJ.out.outJSON
    //    //)
    //}

    CHECK_SATURATION_VDJ(
        ch_genome_bam,
        ch_starsolo_filteredDir,
        ch_whitelist.toList()
    )
    ch_saturation_json = CHECK_SATURATION_VDJ.out.outJSON

    QUALIMAP_VDJ(
        ch_genome_bam,
        ch_genomeGTF
    )
    ch_qualimap_outDir = QUALIMAP_VDJ.out.results

    // aggregate BAM files to check if GEX library exists
    ch_genome_bam
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, meta.expected_cells, file ]
    }
    .groupTuple(by:[0])
    .set{ ch_bam_grouped }

    ch_bc_read
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, meta.expected_cells, file ]
    }
    .groupTuple(by:[0])
    .set{ ch_bcRead_grouped }

    ch_cDNA_read
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, meta.expected_cells, file ]
    }
    .groupTuple(by:[0])
    .set{ ch_cDNAread_grouped }

    ch_starsolo_filteredDir
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, meta.expected_cells, file ]
    }
    .groupTuple(by:[0])
    .set{ ch_filteredDir_grouped }

    TRUST4_VDJ(
        ch_cDNAread_grouped,
        ch_bcRead_grouped,
        ch_bam_grouped,
        ch_filteredDir_grouped,
        ch_vdj_refGenome_fasta,
        ch_vdj_imgt_fasta
    )

    ch_trust4_report        = TRUST4_VDJ.out.report
    ch_trust4_airr          = TRUST4_VDJ.out.airr
    ch_trust4_readsAssign   = TRUST4_VDJ.out.readsAssign
    ch_trust4_barcode       = TRUST4_VDJ.out.barcode
    ch_trust4_finalOut      = TRUST4_VDJ.out.finalOut
    ch_trust4_kneeOut       = TRUST4_VDJ.out.kneeOut
    ch_trust4_cellOut       = TRUST4_VDJ.out.cellOut

    emit:
        starsolo_summary     = ch_starsolo_summary
        starsolo_bam         = ch_genome_bam
        starsolo_umi         = ch_starsolo_umi
        starsolo_filteredDir = ch_starsolo_filteredDir
        qualimap_outDir      = ch_qualimap_outDir
        saturation_json      = ch_saturation_json
        trust4_report        = ch_trust4_report
        trust4_airr          = ch_trust4_airr
        trust4_readsAssign   = ch_trust4_readsAssign
        trust4_barcode       = ch_trust4_barcode
        trust4_finalOut      = ch_trust4_finalOut
        trust4_kneeOut       = ch_trust4_kneeOut
        trust4_cellOut       = ch_trust4_cellOut
}

workflow vdj_report {
    take:
    starsolo_summary
    starsolo_umi
    starsolo_filteredDir
    qualimap_outDir
    saturation_json
    trust4_report
    trust4_airr
    trust4_readsAssign
    trust4_barcode
    trust4_finalOut
    trust4_kneeOut
    trust4_cellOut

    main:

    starsolo_summary
    .map {
        meta, file ->
            //if( meta.feature_types == "GEX" )
                [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ starsolo_summary_collapsed }

    starsolo_umi
    .map {
        meta, file ->
            //if( meta.feature_types == "GEX" )
                [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ starsolo_umi_collapsed }

    starsolo_filteredDir
    .map {
        meta, file ->
            //if( meta.feature_types == "GEX" )
                [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ starsolo_filteredDir_collapsed }

    trust4_report
    .map {
         meta, file ->
            def tmp = []
            def feature = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ vdj_report }

    trust4_airr
    .map {
         meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ vdj_airr }

    trust4_readsAssign
    .map {
         meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ vdj_readsAssign }

    trust4_barcode
    .map {
         meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ vdj_barcode }

    trust4_finalOut
    .map {
         meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ vdj_finalOut }

    trust4_kneeOut
    .map {
         meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ vdj_kneeOut }

    VDJ_METRICS(
        vdj_report,
        vdj_airr,
        vdj_readsAssign,
        vdj_barcode,
        vdj_kneeOut,
        starsolo_summary_collapsed
    )
    
    qualimap_outDir
    .map {
        meta, file ->
            //if( meta.feature_types == "GEX" )
                [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ qualimap_outDir_collapsed }

     saturation_json
    .map {
        meta, file ->
            //if( meta.feature_types == "GEX" )
                [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ saturation_json_collapsed }

    VDJ_METRICS.out.cellOut
    .map {
        meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set { trust4_cells_collapsed }
    //VDJ_METRICS.out.metricsJSON.view()
    
    VDJ_METRICS.out.metricsJSON
    .map {
        meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set { trust4_metrics_collapsed }
    //trust4_metrics_collapsed.view()

    //VDJ_METRICS.out.cloneType.view()
    VDJ_METRICS.out.cloneType
    .map {
        meta, file ->
            def tmp = []
            if(file instanceof List){
                tmp = file
            }else{
                tmp = [ file ]
            }
            tmp.findAll { it =~ /VDJ-[BT]/ }
            .collect {
                if(it =~ /VDJ-T/){
                    ["VDJ-T", it]
                }else if(it =~ /VDJ-B/){
                    ["VDJ-B", it]
                }
            }
            .transpose()
            .plus(0, [id:meta.id])
    }
    .set{ trust4_cloneType_collapsed }
    //trust4_cloneType_collapsed.view()
    
    GET_VERSIONS_VDJ()

    REPORT_VDJ(
        starsolo_summary_collapsed,
        starsolo_umi_collapsed,
        starsolo_filteredDir_collapsed,
        qualimap_outDir_collapsed,
        saturation_json_collapsed,
        vdj_report,
        vdj_airr,
        vdj_kneeOut,
        vdj_finalOut,
        trust4_cells_collapsed,
        trust4_metrics_collapsed,
        trust4_cloneType_collapsed,
        GET_VERSIONS_VDJ.out.versions
    )
}