process STARSOLO_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_high'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/starsolo/${meta.id}_${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/sortedByCoord.out.bam/){
            return null
        }else if(filename=~/Solo.out/){
            return filename.split("/")[-1]
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
    tuple val(meta), path('*d.out.bam')                                                     , emit: bam
    tuple val(meta), path('*d.out.bam.bai')                                                 , emit: bamIndex
    tuple val(meta), path('*Log.final.out')                                                 , emit: log_final
    tuple val(meta), path('*Log.out')                                                       , emit: log_out
    tuple val(meta), path('*Log.progress.out')                                              , emit: log_progress
    tuple val(meta), path("${meta.id}_${meta.feature_types}.matrix_filtered")               , emit: filteredDir,     optional:true
    tuple val(meta), path("${meta.id}_${meta.feature_types}.matrix_raw")                    , emit: rawDir,          optional:true
    tuple val(meta), path("${meta.id}_${meta.feature_types}.Summary.unique.csv")            , emit: summary_unique,  optional:true
    tuple val(meta), path("${meta.id}_${meta.feature_types}.UMIperCellSorted.unique.txt")   , emit: UMI_file_unique, optional:true
    tuple val(meta), path("${meta.id}_${meta.feature_types}.CellReads.stats")               , emit: cellReads_stats, optional:true

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab


    script:
    // check CB_UMI_Complex params
    if(params.soloType == "CB_UMI_Complex"){
        if(!params.soloCBposition){
            exit 1, "ERROR: soloType = 'CB_UMI_Complex' -> Please provide soloCBposition parameter!\n"
        }
        if(!params.soloUMIposition){
            exit 1, "ERROR: soloType = 'CB_UMI_Complex' -> Please provide soloUMIposition parameter!\n"
        }
        if(params.soloCBmatchWLtype!="1MM" && params.soloCBmatchWLtype!="Exact"){
            exit 1, "ERROR: soloType = 'CB_UMI_Complex' -> Please use 1MM or Exact for soloCBmatchWLtype!\n"
        }
    }
    def prefix     = "${meta.id}_${meta.feature_types}"
    //def barcodeMate = params.bc_read == "fastq_1" ? 1 : 2
    // Since starsolo default use single-end mode, activate this label for qualimap option
    meta.single_end = true

    def bamOutTags = ""
    if(params.soloType == "CB_samTagOut"){
        bamOutTags = "NH HI nM AS CR UR CB GX GN gx gn sS sQ sM"
    }else{
        bamOutTags = "NH HI nM AS CR UR CB UB GX GN gx gn sS sQ sM"
    }
    def CBtag = params.whitelist == "None" ? "CR" : "CB"

    def scriptString = []
    scriptString.push(
    """
    STAR --runThreadN $task.cpus \\
    --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --soloBarcodeMate 0 \\
    --readFilesIn $cDNA_read $bc_read \\
    --outFileNamePrefix ${prefix}. \\
    --soloStrand $params.soloStrand \\
    --soloCBwhitelist $whitelist \\
    --soloBarcodeReadLength 0 \\
    --readFilesCommand zcat \\
    --clipAdapterType $params.clipAdapterType \\
    --outFilterScoreMin $params.outFilterScoreMin \\
    --soloCBmatchWLtype $params.soloCBmatchWLtype \\
    --soloUMIfiltering $params.soloUMIfiltering \\
    --soloUMIdedup $params.soloUMIdedup \\
    --soloMultiMappers $params.soloMultiMappers \\
    --soloFeatures $params.soloFeatures \\
    --soloCellReadStats Standard \\
    --outSAMattributes ${bamOutTags} \\
    --outSAMunmapped ${params.outSAMunmapped} \\
    --outSAMtype ${params.outSAMtype} \\
    --outBAMsortingBinsN 300 \\
    --soloCellFilter $params.soloCellFilter \\
    """.stripIndent()
    )
    if(params.soloType == "CB_UMI_Complex"){
    scriptString.push(
    """\
    --soloAdapterSequence $params.soloAdapterSequence \\
    --soloType $params.soloType \\
    --soloCBposition $params.soloCBposition \\
    --soloUMIposition $params.soloUMIposition
    
    """.stripIndent()
    )
    }else{
    scriptString.push(
    """\
    --soloCBstart $params.soloCBstart \\
    --soloCBlen $params.soloCBlen \\
    --soloUMIstart $params.soloUMIstart \\
    --soloUMIlen $params.soloUMIlen \\
    --soloType $params.soloType

    """.stripIndent()
    )
    }

    scriptString.push(
    """
    samtools index -c ${prefix}.Aligned.sortedByCoord.out.bam
    """.stripIndent()
    )
    if(meta.feature_types == "GEX"){
        scriptString.push(
        """
        pigz -p $task.cpus ${prefix}.Solo.out/*/raw/*
        pigz -p $task.cpus ${prefix}.Solo.out/*/filtered/*
        
        ## Rename outputs to prevent file name collision
        cp ${prefix}.Solo.out/*/Summary.csv ${meta.id}_${meta.feature_types}.Summary.unique.csv
        cp ${prefix}.Solo.out/*/UMIperCellSorted.txt ${meta.id}_${meta.feature_types}.UMIperCellSorted.unique.txt
        cp ${prefix}.Solo.out/Gene*/CellReads.stats ${meta.id}_${meta.feature_types}.CellReads.stats
        cp -r ${prefix}.Solo.out/*/filtered ${meta.id}_${meta.feature_types}.matrix_filtered
        cp -r ${prefix}.Solo.out/*/raw ${meta.id}_${meta.feature_types}.matrix_raw
        """.stripIndent()
        )
    }else{
        scriptString.push(
        """
        mkdir ${meta.id}_${meta.feature_types}.matrix_filtered
        summaryFile=\$(find ./ -name "Summary*.csv" | wc -l | awk '{print \$1}')
        if [[ \$summaryFile == 0 ]]
        then
            ## Create fake summary csv file and UMI file
            touch ${meta.id}_${meta.feature_types}.Summary.unique.csv
            touch ${meta.id}_${meta.feature_types}.UMIperCellSorted.unique.txt
            touch ${meta.id}_${meta.feature_types}.CellReads.stats

            ## calculate total reads and valid barcodes
            totalReads=\$(samtools view -@ ${task.cpus} -F 0x100 -F 0x800 ${prefix}.Aligned.sortedByCoord.out.bam | wc -l)
            validBCReads=\$(samtools view -@ ${task.cpus} -F 0x100 -F 0x800 ${prefix}.Aligned.sortedByCoord.out.bam | awk '{cb=\$0; gsub(/.*${CBtag}:Z:/, "", cb); gsub(/\\t.*\$/, "", cb); print \$1"\\t"cb}' | awk '\$2!="-"' | wc -l)
            validBCRatio=\$(awk -v x=\$validBCReads -v y=\$totalReads 'BEGIN{print x/y}')
            echo "Number of Reads,\$totalReads" >>  ${meta.id}_${meta.feature_types}.Summary.unique.csv
            echo "Reads With Valid Barcodes,\$validBCRatio" >>  ${meta.id}_${meta.feature_types}.Summary.unique.csv
        else
            cp ${prefix}.Solo.out/*/Summary.csv ${meta.id}_${meta.feature_types}.Summary.unique.csv
            cp ${prefix}.Solo.out/*/UMIperCellSorted.txt ${meta.id}_${meta.feature_types}.UMIperCellSorted.unique.txt
            cp ${prefix}.Solo.out/Gene*/CellReads.stats ${meta.id}_${meta.feature_types}.CellReads.stats
        fi
        """.stripIndent()
        )
    }
    scriptString.reverse().join()
}

process STARSOLO_MULTIPLE_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/starsolo/${meta.id}_${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(starsolo_rawDir)

    output:
    tuple val(meta), path('raw_mult')     , emit: rawDir
    tuple val(meta), path('filtered_mult'), emit: filteredDir

    script:
    """
    mkdir -p raw_mult
    mkdir -p filtered_mult
    cp ${starsolo_rawDir}/*.tsv.gz raw_mult/
    cp ${starsolo_rawDir}/UniqueAndMult*.mtx.gz raw_mult/matrix.mtx.gz
    gunzip -f raw_mult/*.gz
    STAR --runMode soloCellFiltering raw_mult/ filtered_mult/ \\
    --soloCellFilter $params.soloCellFilter

    pigz -p $task.cpus raw_mult/*
    pigz -p $task.cpus filtered_mult/*
    """
}

process STARSOLO_MULT_SUMMARY_VDJ {
    // Update starsolo summary file for multiple gene reads
    tag "${meta.id}:${meta.feature_types}"
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/starsolo/${meta.id}_${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(cellReads_stats)
    tuple val(meta), path(starsolo_filteredMultiDir)
    tuple val(meta), path(starsolo_unique_summary)
    tuple val(meta), path(saturation_json)

    output:
    tuple val(meta), path('Summary.multiple.csv')   , emit: summary_multiple
    
    script:
    """
    update_starsolo_summary.R ${cellReads_stats} \\
    ${starsolo_filteredMultiDir}/barcodes.tsv.gz \\
    ${starsolo_filteredMultiDir} \\
    ${saturation_json} \\
    ${starsolo_unique_summary} \\
    Summary.multiple.csv
    """
}

process STARSOLO_MULT_UMI_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/starsolo/${meta.id}_${meta.feature_types}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(meta), path(cellReads_stats)

    output:
    tuple val(meta), path('UMIperCellSorted.multiple.txt')   , emit: UMI_file_multiple
    
    script:
    """
    update_umi_file.R ${cellReads_stats} UMIperCellSorted.multiple.txt
    """
}

