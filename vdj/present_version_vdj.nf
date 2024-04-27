process GET_VERSIONS_VDJ {
    tag "GET_VERSIONS_VDJ"
    label 'process_low'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/report",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    output:
    path("versions.json"), emit: json

    script:
    def includeIntron = params.soloFeatures == "Gene" ? "FALSE" : "TRUE"
    def includeMultiReads = params.soloMultiMappers == "Unique" ? "FALSE" : "TRUE"
    def cDNAreadOnly = params.trust4_cDNAread_only ? "TRUE" : "FALSE"
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
    ## trust4 version
    trust4_version=\$(run-trust4 2>&1 > /dev/null | head -1 | awk '{print \$2}')
    cat<<-EOF > versions.json
	{
	  "pipeline_version": "$workflow.manifest.version",
	  "referenceDir": "${params.genomeDir}",
	  "referenceGTF": "${params.genomeGTF}",
      "whitelist": "${params.whitelist}",
	  "STAR_version": "\$star_version",
      "assembleWithcDNAreadOnly": "${cDNAreadOnly}",
      "soloCBmatchWLtype": "${params.soloCBmatchWLtype}",
      "soloUMIfiltering": "${params.soloUMIfiltering}",
      "soloUMIdedup": "${params.soloUMIdedup}",
	  "soloCellFilter": "${params.soloCellFilter}",
	  "includeIntron": "${includeIntron}",
	  "includeMultiReads": "${includeMultiReads}",
	  "samtools_version": "\$samtools_version",
	  "bedtools_version": "\$bedtools_version",
      "trust4_version": "\$trust4_version"
	}
	EOF
    """
}