process TRUST4_VDJ {
    tag "${meta.id}"
    label 'process_high'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/trust4/${meta.id}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/_GEX/){
            return null
        }else if(filename=~/Solo.out/){
            return filename.split("/")[-1]
        }else{
            return filename
        }
    }

    input:
    tuple val(meta), val(cDNAread_featureTypes),    val(cDNAread_expectedCells),    path(cDNAread)
    tuple val(meta), val(bcRead_featureTypes),      val(bcRead_expectedCells),      path(bcRead)
    tuple val(meta), val(starsoloBAM_featureTypes), val(starsoloBAM_expectedCells), path(starsoloBAM)
    tuple val(meta), val(filteredDir_featureTypes), val(filteredDir_expectedCells), path(filteredDir)
    path(trust4_vdj_refGenome_fasta)
    path(trust4_vdj_imgt_fasta)

    output:
    tuple val(meta), path("${meta.id}*_readsAssign.out"),                    emit: readsAssign
    tuple val(meta), path("${meta.id}*_barcode_report.filterDiffusion.tsv"), emit: report
    tuple val(meta), path("${meta.id}*_barcode_airr.tsv"),                   emit: airr
    tuple val(meta), path("${meta.id}*_final.out"),                          emit: finalOut
    tuple val(meta), path("${meta.id}*_barcode.fa"),                         emit: barcode
    tuple val(meta), path("${meta.id}*.kneeOut.tsv"),                        emit: kneeOut
    tuple val(meta), path("${meta.id}*.rawCellOut.tsv"),                     emit: cellOut

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        def dList = []
        if(feature_types.size() == 1 && data_list.getClass() == nextflow.processor.TaskPath){
            dList = [data_list]
        }else{
            dList = data_list
        }
        map = [feature_types, dList].transpose().collectEntries()
        //if(feature_types.size() > 1){
        //  map = [feature_types, data_list].transpose().collectEntries()
        //}else if(feature_types.size() == 1){
        //    map = [feature_types, [data_list]].transpose().collectEntries()
        //}
        return map
    }

    def cDNAread_map      = associate_feature_type(cDNAread_featureTypes,     cDNAread)
    def bcRead_map        = associate_feature_type(bcRead_featureTypes,       bcRead)
    def starsoloBAM_map   = associate_feature_type(starsoloBAM_featureTypes,  starsoloBAM)
    def filteredDir_map   = associate_feature_type(filteredDir_featureTypes,  filteredDir)
    // extract expectedCells from BAM input channel                           
    def expectedCells_map = associate_feature_type(starsoloBAM_featureTypes,  starsoloBAM_expectedCells)
    
    def CBtag = params.whitelist == "None" ? "CR" : "CB"
    def UMItag = params.soloType == "CB_samTagOut" ? "None" : "UB"
    def use_UMI = UMItag == "None" ? "false" : "true"
    def use_cDNAread_only = params.trust4_cDNAread_only ? "true" : "false"
    def scriptString = []
    def vdj_featureTypes = starsoloBAM_featureTypes.collect()
    vdj_featureTypes.remove("GEX")
    vdj_featureTypes.forEach{
        if(starsoloBAM_featureTypes.contains("GEX")){
            scriptString.push(
            """
            ## process bam and generate input reads file
            gex_cells=\$(mktemp -p ./)
            gzip -cd ${filteredDir_map["GEX"]}/barcodes.tsv.gz > \$gex_cells
            
            vdj_cellCalling.sh --inputBAM ${starsoloBAM_map[it]} \\
            --gexBarcode \$gex_cells \\
            --kneeOut ${meta.id}_${it}.kneeOut.tsv \\
            --cellOut ${meta.id}_${it}.rawCellOut.tsv \\
            --readIDout ${it}_readID.lst \\
            --threads ${task.cpus} \\
            --barcode_fasta ${it}_barcode.fa \\
            --umi_fasta ${it}_umi.fa \\
            --CBtag ${CBtag} \\
            --UMItag ${UMItag} \\
            --downSample ${params.trust4_downSample}
            """.stripIndent()
            )
        }else{
            scriptString.push(
            """
            ## process bam and generate input reads file
            vdj_cellCalling.sh --inputBAM ${starsoloBAM_map[it]} \\
            --gexBarcode None \\
            --expectedCells ${expectedCells_map[it]} \\
            --percentile 0.95 \\
            --umi_fold 10 \\
            --kneeOut ${meta.id}_${it}.kneeOut.tsv \\
            --cellOut ${meta.id}_${it}.rawCellOut.tsv \\
            --readIDout ${it}_readID.lst \\
            --threads ${task.cpus} \\
            --barcode_fasta ${it}_barcode.fa \\
            --umi_fasta ${it}_umi.fa \\
            --CBtag ${CBtag} \\
            --UMItag ${UMItag} \\
            --downSample ${params.trust4_downSample}
            
            """.stripIndent()
            )
        }

        scriptString.push(
        """
        ## output barcode fasta
        cp ${it}_barcode.fa ${meta.id}_${it}_barcode.fa 

        ## extract trust4 input reads
        seqtk subseq ${bcRead_map[it]} ${it}_readID.lst | pigz -p 6 > trust4_${it}_input.R1.fq.gz
        seqtk subseq ${cDNAread_map[it]} ${it}_readID.lst | pigz -p 6 > trust4_${it}_input.R2.fq.gz

        use_cDNAread_only=${use_cDNAread_only}
        use_UMI=${use_UMI}
        if [[ \$use_UMI == true ]]
        then
            barcode_umi_opt="--barcode ${it}_barcode.fa --UMI ${it}_umi.fa"
        else
            barcode_umi_opt="--barcode ${it}_barcode.fa"
        fi

        if [[ \$use_cDNAread_only == false ]]
        then
            trust4_input_opt="--inputR1 trust4_${it}_input.R1.fq.gz --inputR2 trust4_${it}_input.R2.fq.gz"
        else
            trust4_input_opt="--inputR2 trust4_${it}_input.R2.fq.gz --R2_Only"
        fi

        singleBC_trust4.sh --genomeRef ${trust4_vdj_refGenome_fasta} \\
               --imgtRef  ${trust4_vdj_imgt_fasta} \\
               \$barcode_umi_opt \\
               \$trust4_input_opt \\
               --readFormat ${params.trust4_readFormat} \\
               --threads ${task.cpus} \\
               --reportOut ${meta.id}_${it}_barcode_report.tsv \\
               --airrOut ${meta.id}_${it}_barcode_airr.tsv \\
               --readsAssign ${meta.id}_${it}_readsAssign.out \\
               --annotFasta ${meta.id}_${it}_annot.fa \\
               --finalOut ${meta.id}_${it}_final.out

        ## filter barcode by diffusion clonetypes
        barcoderep-filter.py -b ${meta.id}_${it}_barcode_report.tsv \\
        -a ${meta.id}_${it}_annot.fa > ${meta.id}_${it}_barcode_report.filterDiffusion.tsv
        """.stripIndent()
        )
    }

    scriptString.reverse().join("\n")
}

process VDJ_METRICS {
    tag "${meta.id}"
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/trust4/${meta.id}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/_GEX/){
            return null
        }else if(filename=~/Solo.out/){
            return filename.split("/")[-1]
        }else{
            return filename
        }
    }

    input:
    tuple val(meta), val(report_featureTypes),          path(report)
    tuple val(meta), val(airr_featureTypes),            path(airr)
    tuple val(meta), val(readsAssign_featureTypes),     path(readsAssign)
    tuple val(meta), val(barcode_featureTypes),         path(barcode)
    tuple val(meta), val(kneeOut_featureTypes),         path(kneeOut)
    tuple val(meta), val(starsoloSummary_featureTypes), path(starsoloSummary)

    output:
    tuple val(meta), path("${meta.id}_*.vdj_cellOut.tsv"),       emit: cellOut
    tuple val(meta), path("${meta.id}_*.vdj_metrics.json"),  emit: metricsJSON
    tuple val(meta), path("${meta.id}_*.cloneType_out.tsv"), emit: cloneType

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        def dList = []
        if(feature_types.size() == 1 && data_list.getClass() == nextflow.processor.TaskPath){
            dList = [data_list]
        }else{
            dList = data_list
        }
        map = [feature_types, dList].transpose().collectEntries()
        //if(feature_types.size() > 1){
        //    map = [feature_types, data_list].transpose().collectEntries()
        //}else if(feature_types.size() == 1){
        //    map = [feature_types, [data_list]].transpose().collectEntries()
        //}//else{
            // return empty, the program will stop and throw the warnning
            //def map = [:]
        //}
        return map
    }
    def report_map          = associate_feature_type(report_featureTypes, report)
    def airr_map            = associate_feature_type(airr_featureTypes, airr)
    def readsAssign_map     = associate_feature_type(readsAssign_featureTypes, readsAssign)
    def barcode_map         = associate_feature_type(barcode_featureTypes, barcode)
    def kneeOut_map         = associate_feature_type(kneeOut_featureTypes, kneeOut)
    def starsoloSummary_map = associate_feature_type(starsoloSummary_featureTypes, starsoloSummary)

    def scriptString = []
    def vdj_featureTypes = report_featureTypes.collect()
    vdj_featureTypes.forEach{
        scriptString.push(
        """
        trust4_metrics.sh ${report_map[it]} \\
                          ${airr_map[it]} \\
                          ${readsAssign_map[it]} \\
                          ${barcode_map[it]} \\
                          ${kneeOut_map[it]} \\
                          ${it} \\
                          ${starsoloSummary_map[it]} \\
                          ${meta.id}_${it}.vdj_cellOut.tsv \\
                          ${meta.id}_${it}.vdj_metrics.json \\
                          ${meta.id}_${it}.cloneType_out.tsv
        """.stripIndent()
        )
    }
    scriptString.reverse().join("\n")
}
