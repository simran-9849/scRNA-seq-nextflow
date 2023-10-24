process GET_VERSIONS_VDJ {
    tag "GET_VERSIONS_VDJ"
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    output:
    path("versions.json"), emit: versions

    script:
    def includeIntron = params.soloFeatures == "Gene" ? "FALSE" : "TRUE"
    def includeMultiReads = params.soloMultiMappers == "Unique" ? "FALSE" : "TRUE"
    """
    ## fastqc version
    fastqc_version=\$(fastqc --version | awk '{print \$2}')
    ## cutadapt version
    cutadapt_version=\$(cutadapt --version | sed '/^\$/d')
    ## STAR version
    star_version=\$(STAR --version)
    ## samtools version
    samtools_version=\$(samtools --version | head -1 |awk '{print \$2}')
    ## bedtools version
    bedtools_version=\$(bedtools --version | awk '{print \$2}')
    ## qualimap version
    qualimap_version=\$(qualimap -h | awk '\$1=="QualiMap"{print \$2}')
    cat<<-EOF > versions.json
	{
	  "pipeline_version": "$workflow.manifest.version",
	  "referenceDir": "${params.genomeDir}",
	  "referenceGTF": "${params.genomeGTF}",
	  "STAR_version": "\$star_version",
      "soloCBmatchWLtype": "${params.soloCBmatchWLtype}",
      "soloUMIfiltering": "${params.soloUMIfiltering}",
      "soloUMIdedup": "${params.soloUMIdedup}",
	  "soloCellFilter": "${params.soloCellFilter}",
	  "includeIntron": "${includeIntron}",
	  "includeMultiReads": "${includeMultiReads}",
	  "samtools_version": "\$samtools_version",
	  "bedtools_version": "\$bedtools_version",
	  "qualimap_version": "\$qualimap_version"
	}
	EOF
    """
}