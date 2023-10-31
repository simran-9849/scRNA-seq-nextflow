process CHECK_SATURATION_VDJ {
    tag "${meta.id}:${meta.feature_types}"
    label 'process_high'
    cache 'lenient'
    fair true
    publishDir "${params.outdir}/saturation/",
        mode: "${params.publish_dir_mode}",
        enabled: params.publishSaturation as boolean,
        saveAs: { filename ->
          if(filename=~/_VDJ/){
            return null
          }else{
            return filename
          }
        }

    input:
    tuple val(meta), path(starsoloBAM)
    tuple val(meta), path(starsolo_filteredDir)
    path(whitelist)

    output:
    tuple val(meta), path("${meta.id}_${meta.feature_types}.saturation_out.json"), emit: outJSON

    script:
    def multiMapper = params.soloMultiMappers == "Unique" ? "unique" : "multiple"
    def whitelist_files = whitelist.join(",")
    if(meta.feature_types =~ /GEX/){
    """
    whitelist_files="${whitelist_files}"
    if [[ \$whitelist_files =~ "," ]]
    then
        whitelist_combined=\$(mktemp -p ./)
        ##https://stackoverflow.com/questions/10586153/how-to-split-a-string-into-an-array-in-bash
        IFS=', ' read -r -a array <<< "\$whitelist_files"
        if [[ \${#array[@]} -eq 2 ]]
        then
            awk 'ARGIND==1{a[\$1]}ARGIND==2{b[\$1]}END{for(m in a){for(n in b){print m"_"n}}}' \${array[0]} \${array[1]} > \$whitelist_combined
        elif [[ \${#array[@]} -eq 3 ]]
        then
            awk 'ARGIND==1{a[\$1]}ARGIND==2{b[\$1]}ARGIND==3{c[\$1]}END{for(m in a){for(n in b){for(k in c){print m"_"n"_"k}}}}' \${array[0]} \${array[1]} \${array[2]} > \$whitelist_combined
        fi
    else
        whitelist_combined="\$whitelist_files"
    fi
    cellFile=\$(mktemp -p ./)
    zcat ${starsolo_filteredDir}/barcodes.tsv.gz > \$cellFile
    get_sequencing_saturation.sh \$whitelist_combined \$cellFile ${multiMapper} ${starsoloBAM} ${task.cpus} \\
                                 ${meta.id}_${meta.feature_types}.saturation_data.json \\
                                 ${meta.id}_${meta.feature_types}.UMI_hist.tsv \\
                                 ${meta.id}_${meta.feature_types}.gene_hist.tsv \\
                                 ${meta.id}_${meta.feature_types}.totalGeneCount.tsv
    rm \$cellFile
    combine_saturation_data.R ${meta.id}_${meta.feature_types} \\
                              ${meta.id}_${meta.feature_types}.saturation_data.json \\
                              ${meta.id}_${meta.feature_types}.UMI_hist.tsv \\
                              ${meta.id}_${meta.feature_types}.gene_hist.tsv \\
                              ${meta.id}_${meta.feature_types}.totalGeneCount.tsv \\
                              ${meta.id}_${meta.feature_types}.saturation_out.json
    """
    }else{
    """
    touch ${meta.id}_${meta.feature_types}.saturation_out.json
    """
    }
}
