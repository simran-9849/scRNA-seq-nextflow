process GENECOVERAGE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path gtf

    output:
    tuple val(meta), path('*.scaled.tab')        ,emit: matrix

    script:
    def prefix     = "${meta.id}"
    """
    ## Convert bam result to bigwiggle
    bamCoverage --bam $bam -o ${prefix}.bw -p $task.cpus
    ## Compute transcript coverage with 100 bins
    computeMatrix scale-regions -S ${prefix}.bw \\
    -R $gtf \\
    --regionBodyLength 1000 \\
    -o ${prefix}.matrix.mat.gz \\
    --outFileNameMatrix ${prefix}.scaled.tab \\
    -p $task.cpus \\
    -q
    """
}

process FEATURESTATS {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path gtf

    output:
    tuple val(meta), path('*.featureCoverage_stats.json')  , emit: stats

    script:
    """
    ## summarized transcripts
    featureCounts -T $task.cpus -O -a $gtf -o ${meta.id}_featureCounts_transcript.tsv -t transcript $bam
    featureCounts -T $task.cpus -O -a $gtf -o ${meta.id}_featureCounts_exon.tsv -t exon $bam

    intergenicReads=\$(awk '\$1=="Unassigned_NoFeatures"{print \$2}' ${meta.id}_featureCounts_transcript.tsv.summary)
    geneReads=\$(awk '\$1=="Assigned"{print \$2}' ${meta.id}_featureCounts_transcript.tsv.summary)
    exonReads=\$(awk '\$1=="Assigned"{print \$2}' ${meta.id}_featureCounts_exon.tsv.summary)
    intronReads=\$(awk -v geneReads=\$geneReads -v exonReads=\$exonReads 'BEGIN{print geneReads-exonReads}')
    totalMappedReads=\$(awk -v geneReads=\$geneReads -v intergenicReads=\$intergenicReads 'BEGIN{print geneReads+intergenicReads}')
    exonRatio=\$(awk -v totalMappedReads=\$totalMappedReads -v exonReads=\$exonReads 'BEGIN{print exonReads/totalMappedReads}')
    intronRatio=\$(awk -v totalMappedReads=\$totalMappedReads -v intronReads=\$intronReads 'BEGIN{print intronReads/totalMappedReads}')
    intergenicRatio=\$(awk -v totalMappedReads=\$totalMappedReads -v intergenicReads=\$intergenicReads 'BEGIN{print intergenicReads/totalMappedReads}')

    jq -n \\
    --arg totalMappedReads "\$totalMappedReads" \\
    --arg exonRatio "\$exonRatio" \\
    --arg intronRatio "\$intronRatio" \\
    --arg intergenicRatio "\$intergenicRatio" \\
    '{sampleName: \"${meta.id}\", totalReads: \$totalMappedReads, exonRatio: \$exonRatio, intronRatio: \$intronRatio, intergenicRatio: \$intergenicRatio}' > ${meta.id}.featureCoverage_stats.json
    """
}
