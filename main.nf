#!/usr/bin/env nextflow

// A lot of codes were directly inherited from https://github.com/nf-core/rnaseq

nextflow.enable.dsl=2

// NF-CORE MODULES

include { CAT_TRIM_FASTQ } from './cat_trim_fastq' addParams( options: ['publish_files': false] )
include { STARSOLO; STARSOLO_COMPLEX; STAR_MKREF; STARSOLO_MULTIPLE; STARSOLO_MULT_SUMMARY; STARSOLO_MULT_UMI } from "./starsolo"
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './modules/nf-core_rnaseq/functions'
include { QUALIMAP_RNASEQ } from './modules/nf-core/modules/qualimap/rnaseq/main'
include { CHECK_SATURATION } from "./sequencing_saturation"
include { GET_VERSIONS } from "./present_version"
include { REPORT } from "./report"


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
    // default workflow scRNA-seq
    scRNAseq()
}

workflow scRNAseq {
    // check mandatory params
    if (!params.input) { exit 1, 'Input samplesheet not specified!' }
    if (!params.genomeDir) { exit 1, 'Genome index DIR not specified!' }
    if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }

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
    ch_bc_read = Channel.empty()
    ch_cDNA_read = Channel.empty()
    CAT_TRIM_FASTQ( ch_fastq )
    if ( params.bc_read == "fastq_1" ){
        ch_bc_read = CAT_TRIM_FASTQ.out.read1
        ch_cDNA_read = CAT_TRIM_FASTQ.out.read2
    }else{
        ch_bc_read = CAT_TRIM_FASTQ.out.read2
        ch_cDNA_read = CAT_TRIM_FASTQ.out.read1
    }

    ch_genomeDir = file(params.genomeDir)
    ch_genomeGTF = file(params.genomeGTF)
    ch_whitelist = file(params.whitelist)
    ch_barcodelist = Channel.fromPath(params.barcodelist.split(" ").toList())

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_starsolo_out               = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    if(params.soloType == "CB_UMI_Complex"){
        STARSOLO_COMPLEX(
            ch_cDNA_read,
            ch_bc_read,
            ch_genomeDir,
            ch_genomeGTF,
            ch_barcodelist
        )
        ch_genome_bam       = STARSOLO_COMPLEX.out.bam
        ch_genome_bam_index = STARSOLO_COMPLEX.out.bai
        ch_filteredDir      = STARSOLO_COMPLEX.out.filteredDir
        ch_starsolo_summary = STARSOLO_COMPLEX.out.summary_unique
        ch_starsolo_UMI     = STARSOLO_COMPLEX.out.UMI_file_unique
    }else{
        STARSOLO(
            ch_cDNA_read,
            ch_bc_read,
            ch_genomeDir,
            ch_genomeGTF,
            ch_whitelist
        )
        ch_genome_bam       = STARSOLO.out.bam
        ch_genome_bam_index = STARSOLO.out.bai
        ch_filteredDir      = STARSOLO.out.filteredDir
        ch_starsolo_summary = STARSOLO.out.summary_unique
        ch_starsolo_UMI     = STARSOLO.out.UMI_file_unique
    }

    if(params.soloMultiMappers != "Unique"){
        STARSOLO_MULTIPLE(
            STARSOLO.out.rawDir
        )
        CHECK_SATURATION(
            ch_genome_bam,
            STARSOLO_MULTIPLE.out.filteredDir,
            ch_whitelist
        )
        STARSOLO_MULT_SUMMARY(
            STARSOLO.out.cellReads_stats,
            STARSOLO_MULTIPLE.out.filteredDir,
            STARSOLO.out.summary_unique,
            CHECK_SATURATION.out.outJSON
        )
        STARSOLO_MULT_UMI(
            STARSOLO.out.cellReads_stats   
        )
        GET_VERSIONS(
            CHECK_SATURATION.out.outJSON
        )
    }else{
        CHECK_SATURATION(
            ch_genome_bam,
            ch_filteredDir,
            ch_whitelist
        )
        GET_VERSIONS(
            CHECK_SATURATION.out.outJSON
        )
    }

    ch_qualimap_multiqc           = Channel.empty()
    QUALIMAP_RNASEQ(
        ch_genome_bam,
        ch_genomeGTF
    )
    ch_qualimap_multiqc = QUALIMAP_RNASEQ.out.results

    if(params.soloMultiMappers != "Unique"){
        REPORT(
            STARSOLO_MULT_SUMMARY.out.summary_multiple,
            STARSOLO_MULT_UMI.out.UMI_file_multiple,
            STARSOLO_MULTIPLE.out.filteredDir,
            ch_qualimap_multiqc,
            CHECK_SATURATION.out.outJSON,
            GET_VERSIONS.out.versions
        )
    }else{
        REPORT(
            ch_starsolo_summary,
            ch_starsolo_UMI,
            ch_filteredDir,
            ch_qualimap_multiqc,
            CHECK_SATURATION.out.outJSON,
            GET_VERSIONS.out.versions
        )
    }
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
include { STARSOLO_VDJ; STARSOLO_MULTIPLE_VDJ; STARSOLO_MULT_SUMMARY_VDJ; STARSOLO_MULT_UMI_VDJ; STARSOLO_COMPLEX_VDJ } from "./vdj/starsolo_vdj"
include { QUALIMAP_VDJ } from "./vdj/qualimap_vdj"
include { CHECK_SATURATION_VDJ } from "./vdj/sequencing_saturation_vdj.nf"
include { TRUST4_VDJ; VDJ_METRICS } from "./vdj/trust4_vdj"
include { GET_VERSIONS_VDJ } from "./vdj/present_version_vdj"
include { REPORT_VDJ } from "./vdj/report_vdj"

def create_fastq_channel_featureTypes(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.bc_read}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.cDNA_read}"
    }

    if(!row.feature_types){
        exit 1, "ERROR: Please check input samplesheet -> feature_types column does not exist!\n"
    }

    meta.feature_types = row.feature_types
    array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]

    return array
}

workflow vdj {
    vdj_process()
    vdj_report(
        vdj_process.out.starsolo_summary,
        vdj_process.out.starsolo_bam,
        vdj_process.out.starsolo_umi,
        vdj_process.out.starsolo_filteredDir,
        vdj_process.out.qualimap_outDir,
        vdj_process.out.saturation_json
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
        .map{ create_fastq_channel_featureTypes(it) }
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
    ch_whitelist = file(params.whitelist)

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

    if(params.soloType == "CB_UMI_Complex"){
        // barcodelist will only be required for complex mode
        ch_barcodelist = Channel.fromPath(params.barcodelist.split(" ").toList())
        STARSOLO_COMPLEX_VDJ(
            ch_cDNA_read,
            ch_bc_read,
            ch_genomeDir,
            ch_genomeGTF,
            ch_barcodelist.toList(),
        )
        ch_genome_bam           = STARSOLO_COMPLEX_VDJ.out.bam
        ch_genome_bam_index     = STARSOLO_COMPLEX_VDJ.out.bai
        ch_starsolo_filteredDir = STARSOLO_COMPLEX_VDJ.out.filteredDir
        ch_starsolo_summary     = STARSOLO_COMPLEX_VDJ.out.summary_unique
        ch_starsolo_umi         = STARSOLO_COMPLEX_VDJ.out.UMI_file_unique

    }else{
        STARSOLO_VDJ(
            ch_cDNA_read,
            ch_bc_read,
            ch_genomeDir,
            ch_genomeGTF,
            ch_whitelist,
        )
        ch_genome_bam           = STARSOLO_VDJ.out.bam
        ch_genome_bam_index     = STARSOLO_VDJ.out.bai
        ch_starsolo_filteredDir = STARSOLO_VDJ.out.filteredDir
        ch_starsolo_summary     = STARSOLO_VDJ.out.summary_unique
        ch_starsolo_umi         = STARSOLO_VDJ.out.UMI_file_unique

    }
    
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
        ch_whitelist
    )
    ch_saturation_json = CHECK_SATURATION_VDJ.out.outJSON

    QUALIMAP_VDJ(
        ch_genome_bam,
        ch_genomeGTF
    )
    ch_qualimap_outDir = QUALIMAP_VDJ.out.results

    emit:
        starsolo_summary     = ch_starsolo_summary
        starsolo_bam         = ch_genome_bam
        starsolo_umi         = ch_starsolo_umi
        starsolo_filteredDir = ch_starsolo_filteredDir
        qualimap_outDir      = ch_qualimap_outDir
        saturation_json      = ch_saturation_json

}

workflow vdj_report {
    take:
    starsolo_summary
    starsolo_bam
    starsolo_umi
    starsolo_filteredDir
    qualimap_outDir
    saturation_json

    main:
    ch_vdj_refGenome_fasta = file(params.trust4_vdj_refGenome_fasta)
    ch_vdj_imgt_fasta = file(params.trust4_vdj_imgt_fasta)
    ch_whitelist = file(params.whitelist)
    
    //starsolo_bam.view()
    TRUST4_VDJ(
        starsolo_bam,
        starsolo_summary,
        ch_vdj_refGenome_fasta,
        ch_vdj_imgt_fasta,
        ch_whitelist
    )

    ch_trust4_report = TRUST4_VDJ.out.trust4_report
    ch_trust4_airr   = TRUST4_VDJ.out.trust4_airr
    ch_toassemble_bc = TRUST4_VDJ.out.toassemble_bc

    ch_trust4_report
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ vdj_report }

    ch_trust4_airr
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ vdj_airr }

    ch_toassemble_bc
    .map {
        meta, file ->
        [ [id:meta.id], meta.feature_types, file ]
    }
    .groupTuple(by:[0])
    .set{ vdj_toassemble_bc }

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

    VDJ_METRICS(
        vdj_report,
        vdj_airr,
        vdj_toassemble_bc,
        starsolo_summary_collapsed,
        starsolo_filteredDir_collapsed
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

    VDJ_METRICS.out.metricsJSON.view()
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
    trust4_metrics_collapsed.view()

    //VDJ_METRICS.out.kneeData.view()
    VDJ_METRICS.out.kneeData
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
    .set{ trust4_kneeData_collapsed }
    //trust4_kneeData_collapsed.view()

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

    //starsolo_summary_collapsed.view()
    REPORT_VDJ(
        starsolo_summary_collapsed,
        starsolo_umi_collapsed,
        starsolo_filteredDir_collapsed,
        qualimap_outDir_collapsed,
        saturation_json_collapsed,
        trust4_metrics_collapsed,
        trust4_kneeData_collapsed,
        trust4_cloneType_collapsed,
        GET_VERSIONS_VDJ.out.versions
    )
}