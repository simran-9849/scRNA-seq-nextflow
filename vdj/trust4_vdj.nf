process VDJ_CELLCALLING_WITHOUTGEX {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_high'
    cache 'lenient'
    fair true

    input:
    tuple val(meta), path(cDNAread), path(bcRead), path(starsoloBAM)

    output:
    tuple val(meta), path(cDNAread),                                         emit: cDNAread
    tuple val(meta), path(bcRead),                                           emit: bcRead
    tuple val(meta), path("${meta.id}*.kneeOut.tsv"),                        emit: kneeOut
    tuple val(meta), path("${meta.id}*.rawCellOut.tsv"),                     emit: cellOut
    tuple val(meta), path("${meta.id}*_readID.lst"),                         emit: readList
    tuple val(meta), path("${meta.id}*_barcode.fa"),                         emit: barcode
    tuple val(meta), path("${meta.id}*_umi.fa"),                             emit: umi

    script:
    def prefix = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    def CBtag = params.whitelist == "None" ? "CR" : "CB"
    def UMItag = params.soloType == "CB_samTagOut" ? "None" : "UB"

    """
    ## process bam and generate input reads file
    vdj_cellCalling.sh --inputBAM ${starsoloBAM} \\
    --gexBarcode None \\
    --expectedCells ${meta.expected_cells} \\
    --percentile 0.95 \\
    --umi_fold 10 \\
    --kneeOut ${prefix}.kneeOut.tsv \\
    --cellOut ${prefix}.rawCellOut.tsv \\
    --readIDout ${prefix}_readID.lst \\
    --threads ${task.cpus} \\
    --barcode_fasta ${prefix}_barcode.fa \\
    --umi_fasta ${prefix}_umi.fa \\
    --CBtag ${CBtag} \\
    --UMItag ${UMItag} \\
    --downSample ${params.trust4_downSample}
    """
}

process VDJ_CELLCALLING_WITHGEX {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_high'
    cache 'lenient'
    fair true

    input:
    tuple val(meta), path(cDNAread), path(bcRead), path(starsoloBAM), path(filteredDir)

    output:
    tuple val(meta), path(cDNAread),                                         emit: cDNAread
    tuple val(meta), path(bcRead),                                           emit: bcRead
    tuple val(meta), path("${meta.id}*.kneeOut.tsv"),                        emit: kneeOut
    tuple val(meta), path("${meta.id}*.rawCellOut.tsv"),                     emit: cellOut
    tuple val(meta), path("${meta.id}*_readID.lst"),                         emit: readList
    tuple val(meta), path("${meta.id}*_barcode.fa"),                         emit: barcode
    tuple val(meta), path("${meta.id}*_umi.fa"),                             emit: umi

    script:
    def prefix = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    def CBtag = params.whitelist == "None" ? "CR" : "CB"
    def UMItag = params.soloType == "CB_samTagOut" ? "None" : "UB"
    def use_UMI = UMItag == "None" ? "false" : "true"
    def use_cDNAread_only = params.trust4_cDNAread_only ? "true" : "false"

    """
    ## process bam and generate input reads file
    gex_cells=\$(mktemp -p ./)
    gzip -cd ${filteredDir}/barcodes.tsv.gz > \$gex_cells
    
    vdj_cellCalling.sh --inputBAM ${starsoloBAM} \\
    --gexBarcode \$gex_cells \\
    --kneeOut ${prefix}.kneeOut.tsv \\
    --cellOut ${prefix}.rawCellOut.tsv \\
    --readIDout ${prefix}_readID.lst \\
    --threads ${task.cpus} \\
    --barcode_fasta ${prefix}_barcode.fa \\
    --umi_fasta ${prefix}_umi.fa \\
    --CBtag ${CBtag} \\
    --UMItag ${UMItag} \\
    --downSample ${params.trust4_downSample}
    """
}

process VDJ_ASSEMBLY {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_high'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/${meta.id}/trust4/${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta),
          path(cDNAread),
          path(bcRead),
          path(readList),
          path(barcodeFASTA),
          path(umiFASTA),
          path(kneeOut),
          path(cellOut)
    path(trust4_vdj_refGenome_fasta)
    path(trust4_vdj_imgt_fasta)

    output:
    tuple val(meta), path("${meta.id}*_readsAssign.out"),                    emit: readsAssign
    tuple val(meta), path("${meta.id}*_barcode_report.filterDiffusion.tsv"), emit: report
    tuple val(meta), path("${meta.id}*_barcode_airr.tsv"),                   emit: airr
    tuple val(meta), path("${meta.id}*_final.out"),                          emit: finalOut

    script:
    def prefix = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    def UMItag = params.soloType == "CB_samTagOut" ? "None" : "UB"
    def use_UMI = UMItag == "None" ? "false" : "true"
    def use_cDNAread_only = params.trust4_cDNAread_only ? "true" : "false"

    """
    ## extract trust4 input reads
    seqtk subseq ${bcRead} ${readList} | pigz -p 6 > trust4_${prefix}_input.R1.fq.gz
    seqtk subseq ${cDNAread} ${readList} | pigz -p 6 > trust4_${prefix}_input.R2.fq.gz

    ## extraction will fail if readsID tails with "/1" or "/2"
    extracted_reads=\$(zcat trust4_${prefix}_input.R1.fq.gz | wc -l)
    if [[ \$extracted_reads -eq 0 ]]
    then
        zcat ${bcRead} |
            awk 'NR%4==1{sub("/1", "", \$1); print \$1; getline; print; getline; print; getline; print}' |
            seqtk subseq - ${readList} | pigz -p 6 > trust4_${prefix}_input.R1.fq.gz
    fi
    extracted_reads=\$(zcat trust4_${prefix}_input.R2.fq.gz | wc -l)
    if [[ \$extracted_reads -eq 0 ]]
    then
        zcat ${cDNAread} |
            awk 'NR%4==1{sub("/2", "", \$1); print \$1; getline; print; getline; print; getline; print}' |
            seqtk subseq - ${readList} | pigz -p 6 > trust4_${prefix}_input.R2.fq.gz
    fi

    use_cDNAread_only=${use_cDNAread_only}
    use_UMI=${use_UMI}
    if [[ \$use_UMI == true ]]
    then
        barcode_umi_opt="--barcode ${barcodeFASTA} --UMI ${umiFASTA}"
    else
        barcode_umi_opt="--barcode ${barcodeFASTA}"
    fi

    if [[ \$use_cDNAread_only == false ]]
    then
        trust4_input_opt="--inputR1 trust4_${prefix}_input.R1.fq.gz --inputR2 trust4_${prefix}_input.R2.fq.gz"
    else
        trust4_input_opt="--inputR2 trust4_${prefix}_input.R2.fq.gz --R2_Only"
    fi

    singleBC_trust4.sh --genomeRef ${trust4_vdj_refGenome_fasta} \\
           --imgtRef  ${trust4_vdj_imgt_fasta} \\
           \$barcode_umi_opt \\
           \$trust4_input_opt \\
           --readFormat ${params.trust4_readFormat} \\
           --threads ${task.cpus} \\
           --reportOut ${prefix}_barcode_report.tsv \\
           --airrOut ${prefix}_barcode_airr.tsv \\
           --readsAssign ${prefix}_readsAssign.out \\
           --annotFasta ${prefix}_annot.fa \\
           --finalOut ${prefix}_final.out

    ## filter barcode by diffusion clonetypes
    barcoderep-filter.py -b ${prefix}_barcode_report.tsv \\
                         -a ${prefix}_annot.fa > ${prefix}_barcode_report.filterDiffusion.tsv
    """
}

process VDJ_METRICS {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/${meta.id}/trust4/${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta),
          path(report),
          path(airr),
          path(readsAssign),
          path(barcode),
          path(kneeOut),
          path(starsoloSummary)

    output:
    tuple val(meta), path("${meta.id}_*.vdj_cellOut.tsv"),       emit: cellOut
    tuple val(meta), path("${meta.id}_*.vdj_metrics.json"),      emit: metricsJSON
    tuple val(meta), path("${meta.id}_*.cloneType_out.tsv"),     emit: cloneType

    script:
    def prefix   = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    """
    trust4_metrics.sh ${report} \\
                      ${airr} \\
                      ${readsAssign} \\
                      ${barcode} \\
                      ${kneeOut} \\
                      ${meta.feature_types} \\
                      ${starsoloSummary} \\
                      ${prefix}.vdj_cellOut.tsv \\
                      ${prefix}.vdj_metrics.json \\
                      ${prefix}.cloneType_out.tsv
    """
}
