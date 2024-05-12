process STARSOLO {
    tag { meta.feature_types ? "${meta.id}:${meta.feature_types}" : "${meta.id}" }
    label 'process_high'
    publishDir { meta.feature_types ? "${params.outdir}/${meta.id}/starsolo/${meta.feature_types}" : "${params.outdir}/${meta.id}/starsolo"},
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(!params.publishBAM && filename=~/sortedByCoord.out.bam/){
            return null
        }else if(filename=~/Solo.out/){
            return filename.split("/")[-1]
        }else{
            return filename
        }
    }

    input:
    tuple val(meta), path(bc_read), path(cDNA_read)
    path index
    path gtf
    path whitelist

    output:
    tuple val(meta), path('*d.out.bam')                                            , emit: bam
    tuple val(meta), path('*d.out.bam.bai')                                        , emit: bai
    tuple val(meta), path('*Log.final.out')                                        , emit: log_final
    tuple val(meta), path('*Log.out')                                              , emit: log_out
    tuple val(meta), path('*Log.progress.out')                                     , emit: log_progress
    tuple val(meta), path('*.matrix_filtered')                      , optional:true, emit: filteredDir
    tuple val(meta), path('*.matrix_raw')                           , optional:true, emit: rawDir
    tuple val(meta), path('*_summary.unique.csv')                   , optional:true, emit: summary_unique
    tuple val(meta), path('*_UMIperCellSorted.unique.txt')          , optional:true, emit: UMI_file_unique
    tuple val(meta), path('*.CellReads.stats')                      , optional:true, emit: cellReads_stats

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab


    script:
    def prefix     = meta.feature_types ? "${meta.id}_${meta.feature_types}" : "${meta.id}"
    //def barcodeMate = params.bc_read == "fastq_1" ? 1 : 2
    // Since starsolo default use single-end mode, activate this label for qualimap option
    //meta.single_end = true

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

    def bamOutTags = "NH HI nM AS CR UR CB UB GX GN gx gn sS sQ sM"
    if(params.soloType == "CB_samTagOut"){
        bamOutTags = "NH HI nM AS CR UR CB GX GN gx gn sS sQ sM"
    }
    def CBtag = params.whitelist == "None" ? "CR" : "CB"

    def soloCellFilter = "$params.soloCellFilter"
    if(params.soloCellFilter=="EmptyDrops_CR"){
        soloCellFilter = "EmptyDrops_CR ${meta.expected_cells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000"
    }else if(params.soloCellFilter=="CellRanger2.2"){
        soloCellFilter = "CellRanger2.2 ${meta.expected_cells} 0.99 10"
    }

    scriptString = []
    // common starsolo options
    scriptString.push(
    """
    ## Set ulimit to avoid number of opening files reaching to limit
    ulimit -n 40960
    ## Added "--outBAMsortingBinsN 300" option to solve sorting RAM issue when BAM is too large
    ## Refer to: https://github.com/alexdobin/STAR/issues/870
    STAR --runThreadN $task.cpus \\
    --genomeDir $index \\
    --sjdbGTFfile $gtf \\
    --soloBarcodeMate 0 \\
    --readFilesIn $cDNA_read $bc_read \\
    --soloCBwhitelist $whitelist \\
    --soloBarcodeReadLength 0 \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${prefix}. \\
    --soloStrand $params.soloStrand \\
    --clipAdapterType $params.clipAdapterType \\
    --outFilterScoreMin $params.outFilterScoreMin \\
    --soloCBmatchWLtype $params.soloCBmatchWLtype \\
    --soloUMIfiltering $params.soloUMIfiltering \\
    --soloUMIdedup $params.soloUMIdedup \\
    --soloMultiMappers $params.soloMultiMappers \\
    --soloFeatures $params.soloFeatures \\
    --soloCellFilter $soloCellFilter \\
    --soloCellReadStats Standard \\
    --outSAMattributes ${bamOutTags} \\
    --outSAMtype ${params.outSAMtype} \\
    --outSAMunmapped ${params.outSAMunmapped} \\
    --limitBAMsortRAM ${params.limitBAMsortRAM} \\
    --outBAMsortingBinsN 300 \\
    """.stripIndent()
    )

    if(params.soloType == "CB_UMI_Complex"){
    scriptString.push(
    """\
    --soloType $params.soloType \\
    --soloCBposition $params.soloCBposition \\
    --soloUMIposition $params.soloUMIposition \\
    --soloAdapterSequence $params.soloAdapterSequence \\
    --soloAdapterMismatchesNmax $params.soloAdapterMismatchesNmax
    
    """.stripIndent()
    )
    }else if(params.soloType == "CB_UMI_Simple"){
    scriptString.push(
    """\
    --soloType $params.soloType \\
    --soloCBstart $params.soloCBstart \\
    --soloCBlen $params.soloCBlen \\
    --soloUMIstart $params.soloUMIstart \\
    --soloUMIlen $params.soloUMIlen
    
    """.stripIndent()
    )
    }else{
        exit 1, "soloType only support CB_UMI_Simple and CB_UMI_Complex for now."
    }

    scriptString.push(
    """
    samtools index ${prefix}.Aligned.sortedByCoord.out.bam
    """.stripIndent()
    )

    if(!meta.feature_types || meta.feature_types == "GEX"){
        scriptString.push(
        """
        pigz -p $task.cpus ${prefix}.Solo.out/${params.soloFeatures}/raw/*
        pigz -p $task.cpus ${prefix}.Solo.out/${params.soloFeatures}/filtered/*
        cp ${prefix}.Solo.out/${params.soloFeatures}/Summary.csv ${prefix}_summary.unique.csv
        cp ${prefix}.Solo.out/${params.soloFeatures}/UMIperCellSorted.txt ${prefix}_UMIperCellSorted.unique.txt
        cp ${prefix}.Solo.out/${params.soloFeatures}/CellReads.stats ${prefix}.CellReads.stats
        cp -r ${prefix}.Solo.out/${params.soloFeatures}/filtered ${prefix}.matrix_filtered
        cp -r ${prefix}.Solo.out/${params.soloFeatures}/raw ${prefix}.matrix_raw
        """.stripIndent()
        )
    }else{
        scriptString.push(
        """
        mkdir ${prefix}.matrix_filtered
        summaryFile=\$(find ./ -name "Summary*.csv" | wc -l | awk '{print \$1}')
        if [[ \$summaryFile == 0 ]]
        then
            ## Create fake summary csv file and UMI file
            touch ${prefix}_summary.unique.csv
            touch ${prefix}_UMIperCellSorted.unique.txt
            touch ${prefix}.CellReads.stats

            ## calculate total reads and valid barcodes
            totalReads=\$(samtools view -@ ${task.cpus} -F 0x100 -F 0x800 ${prefix}.Aligned.sortedByCoord.out.bam | wc -l)
            validBCReads=\$(samtools view -@ ${task.cpus} -F 0x100 -F 0x800 ${prefix}.Aligned.sortedByCoord.out.bam | awk '{cb=\$0; gsub(/.*${CBtag}:Z:/, "", cb); gsub(/\\t.*\$/, "", cb); print \$1"\\t"cb}' | awk '\$2!="-"' | wc -l)
            validBCRatio=\$(awk -v x=\$validBCReads -v y=\$totalReads 'BEGIN{print x/y}')
            echo "Number of Reads,\$totalReads" >>  ${prefix}_summary.unique.csv
            echo "Reads With Valid Barcodes,\$validBCRatio" >>  ${prefix}_summary.unique.csv
        else
            cp ${prefix}.Solo.out/${params.soloFeatures}/Summary.csv ${prefix}_summary.unique.csv
            cp ${prefix}.Solo.out/${params.soloFeatures}/UMIperCellSorted.txt ${prefix}_UMIperCellSorted.unique.txt
            cp ${prefix}.Solo.out/${params.soloFeatures}/CellReads.stats ${prefix}.CellReads.stats
        fi
        """.stripIndent()
        )
    }
    scriptString.reverse().join()
}

process STAR_MKREF {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: "move",
        enabled: params.outdir as boolean
    input:
    path genomeFasta
    path genomeGTF
    output:
    path("${params.refoutDir}")       , emit: refPath
    script:
    """
    STAR --runMode genomeGenerate \\
    --runThreadN $task.cpus \\
    --genomeDir $params.refoutDir \\
    --genomeFastaFiles $genomeFasta \\
    --sjdbGTFfile $genomeGTF $params.mkrefOpt
    """
}
