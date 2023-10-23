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
    tuple val(meta), path("TRUST4_OUT/${meta.id}*_toassemble_bc.fa"), emit: toassemble_bc
    tuple val(meta), path("TRUST4_OUT/${meta.id}*_barcode_report.tsv"), emit: report
    tuple val(meta), path("TRUST4_OUT/${meta.id}*_barcode_airr.tsv"), emit: airr
    tuple val(meta), path("TRUST4_OUT/${meta.id}*_final.out"), emit: finalOut
    tuple val(meta), path("${meta.id}*_kneeOut.tsv"), emit: kneeOut
    tuple val(meta), path("${meta.id}*_cellOut.tsv"), emit: cellOut

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        if(feature_types.size() > 1){
            map = [feature_types, data_list].transpose().collectEntries()
        }else if(feature_types.size() == 1){
            map = [feature_types, [data_list]].transpose().collectEntries()
        }//else{
            // return empty, them program will stop and throw the warnning
            //def map = [:]
        //}
        return map
    }

    def cDNAread_map      = associate_feature_type(cDNAread_featureTypes,     cDNAread)
    def bcRead_map        = associate_feature_type(bcRead_featureTypes,       bcRead)
    def starsoloBAM_map   = associate_feature_type(starsoloBAM_featureTypes,  starsoloBAM)
    def filteredDir_map   = associate_feature_type(filteredDir_featureTypes,  filteredDir)
    // extract expectedCells from BAM input channel                           
    def expectedCells_map = associate_feature_type(starsoloBAM_featureTypes,  starsoloBAM_expectedCells)

    def striptString = []
    def vdj_featureTypes = starsoloBAM_featureTypes.collect()
    vdj_featureTypes.remove("GEX")
    if(starsoloBAM_featureTypes.contains("GEX")){
        vdj_featureTypes.forEach{
            striptString.push(
            """
            ## process bam and generate input reads file
            gex_cells=\$(mktemp -p ./)
            gzip -cd ${filteredDir_map["GEX"]}/barcodes.tsv.gz > \$gex_cells
            
            vdj_cellCalling.sh --inputBAM ${starsoloBAM_map[it]} \\
            --gexBarcode \$gex_cells \\
            --kneeOut ${meta.id}_${it}_kneeOut.tsv \\
            --cellOut ${meta.id}_${it}_cellOut.tsv \\
            --readIDout ${it}_readID.lst \\
            --barcode_fasta ${it}_barcode.fa \\
            --umi_fasta ${it}_umi.fa \\
            --threads ${task.cpus}
            
            ## extract trust4 input reads
            seqtk subseq ${bcRead_map[it]} ${it}_readID.lst | pigz -p 6 > trust4_${it}_input.R1.fq.gz
            seqtk subseq ${cDNAread_map[it]} ${it}_readID.lst | pigz -p 6 > trust4_${it}_input.R2.fq.gz
            
            run-trust4 -f ${trust4_vdj_refGenome_fasta} \\
            --ref ${trust4_vdj_imgt_fasta} \\
            -o ${meta.id}_${it} \\
            --od TRUST4_OUT \\
            -1 trust4_${it}_input.R1.fq.gz \\
            -2 trust4_${it}_input.R2.fq.gz \\
            --barcode ${it}_barcode.fa \\
            --UMI ${it}_umi.fa \\
            --readFormat r1:59:-1 \\
            --outputReadAssignment \\
            -t ${task.cpus}

            ## filter barcode by diffusion clonetypes
            barcoderep-filter.py -b TRUST4_OUT/${meta.id}_${it}_barcode_report.tsv \\
            -a TRUST4_OUT/${meta.id}_${it}_annot.fa > TRUST4_OUT/${meta.id}_${it}_barcode_report.filterDiffusion.tsv
            """.stripIndent()
            )
        }
    }else{
        vdj_featureTypes.forEach{
            striptString.push(
            """
            ## process bam and generate input reads file
            vdj_cellCalling.sh --inputBAM ${starsoloBAM_map[it]} \\
            --gexBarcode None \\
            --expectedCells ${expectedCells_map[it]} \\
            --percentile 0.95 \\
            --umi_fold 10 \\
            --kneeOut ${meta.id}_${it}_kneeOut.tsv \\
            --cellOut ${meta.id}_${it}_cellOut.tsv \\
            --readIDout ${it}_readID.lst \\
            --barcode_fasta ${it}_barcode.fa \\
            --umi_fasta ${it}_umi.fa \\
            --threads ${task.cpus}
            
            ## extract trust4 input reads
            seqtk subseq ${bcRead_map[it]} ${it}_readID.lst | pigz -p 6 > trust4_${it}_input.R1.fq.gz
            seqtk subseq ${cDNAread_map[it]} ${it}_readID.lst | pigz -p 6 > trust4_${it}_input.R2.fq.gz
            
            run-trust4 -f ${trust4_vdj_refGenome_fasta} \\
            --ref ${trust4_vdj_imgt_fasta} \\
            -o ${meta.id}_${it} \\
            --od TRUST4_OUT \\
            -1 trust4_${it}_input.R1.fq.gz \\
            -2 trust4_${it}_input.R2.fq.gz \\
            --barcode ${it}_barcode.fa \\
            --UMI ${it}_umi.fa \\
            --readFormat r1:59:-1 \\
            --outputReadAssignment \\
            -t ${task.cpus}

            ## filter barcode by diffusion clonetypes
            barcoderep-filter.py -b TRUST4_OUT/${meta.id}_${it}_barcode_report.tsv \\
            -a TRUST4_OUT/${meta.id}_${it}_annot.fa > TRUST4_OUT/${meta.id}_${it}_barcode_report.filterDiffusion.tsv
            """.stripIndent()
            )
        }
    }
    striptString.join("\n")
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
    tuple val(meta), val(toassemble_bc_featureTypes),   path(toassemble_bc)
    tuple val(meta), val(starsoloSummary_featureTypes), path(starsoloSummary)
    tuple val(meta), val(filteredDir_featureTypes),     path(filteredDir)

    output:
    tuple val(meta), path("${meta.id}_*.vdj_metrics.json"), emit: metricsJSON
    tuple val(meta), path("${meta.id}_*.knee_input.tsv"), emit: kneeData
    tuple val(meta), path("${meta.id}_*.cloneType_out.tsv"), emit: cloneType
    tuple val(meta), path("TRUST_${meta.id}_*_barcode_report.filtered.tsv"), emit: trust4_report
    tuple val(meta), path("TRUST_${meta.id}_*_barcode_airr.filtered.tsv"), emit: trust4_airr

    script:
    // https://stackoverflow.com/questions/49114850/create-a-map-in-groovy-having-two-collections-with-keys-and-values
    def associate_feature_type = { feature_types, data_list ->
        def map = [:]
        if(feature_types.size() > 1){
            map = [feature_types, data_list].transpose().collectEntries()
        }else if(feature_types.size() == 1){
            map = [feature_types, [data_list]].transpose().collectEntries()
        }//else{
            // return empty, the program will stop and throw the warnning
            //def map = [:]
        //}
        return map
    }
    def report_map          = associate_feature_type(report_featureTypes, report)
    def airr_map            = associate_feature_type(airr_featureTypes, airr)
    def toassemble_bc_map   = associate_feature_type(toassemble_bc_featureTypes, toassemble_bc)
    def starsoloSummary_map = associate_feature_type(starsoloSummary_featureTypes, starsoloSummary)
    def filteredDir_map     = associate_feature_type(filteredDir_featureTypes, filteredDir)

    def gex_cells           = filteredDir_featureTypes.contains("GEX") ? "${filteredDir_map['GEX']}/barcodes.tsv.gz" : "None"
    def featureArgument     = report_map.findAll{it.key != "GEX"}.keySet()
    def summaryArgument     = featureArgument.collect{ starsoloSummary_map[it] }.join(",")
    featureArgument         = featureArgument.join(",")
    """
    
    ## VDJ-B
    trust4_metrics_filter.sh ${meta.id} \\
                             ${featureArgument} \\
                             ${gex_cells} \\
                             ${summaryArgument}

    """
}