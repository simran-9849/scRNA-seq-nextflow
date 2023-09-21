#! /usr/bin/env bash

sampleID=$1
## featureType="VDJ-B,VDJ-T"
featureType=$2
gex_cells=$3
## starsoloSummary="VDJ-B_summary.tsv,VDJ-T_summary.tsv"
starsoloSummary=$4

filter_airr_by_report(){

    local report_input=$1
    local airr_input=$2
    local airr_output=$3
    awk 'ARGIND==1{split($3, tmp, ","); a[tmp[8]]; split($4, tmp, ","); a[tmp[8]]}ARGIND==2{if(FNR==1){print}else if($1 in a){print}}' $report_input $airr_input > $airr_output
}


filter_VDJ(){

    local gex_cells=$1
    local vdj_type=$2
    local report_input=$3
    local airr_input=$4
    local toassembly_input=$5
    local starsolo_summary=$6

    local report_gexFiltered_out=$(basename $report_input)".gexFiltered"
    local toassembly_gexFiltered_out=$(basename $toassembly_input)".gexFiltered"
    local airr_gexFiltered_out=$(basename $airr_input)".gexFiltered"
    local report_finalFiltered_out=$(basename $report_input)".finalFiltered"
    local gex_barcodes=$(mktemp -p ./)

    if [[  -f $gex_cells ]] && [[ $gex_cells =~ \.gz$ ]]
    then
        zcat $gex_cells > $gex_barcodes
    elif [[  -f $gex_cells ]] && [[ $gex =~ .tsv$ ]]
    then
        cp $gex_cells $gex_barcodes
    fi

    if [[  $(wc -l $gex_barcodes | awk '{print $1}') > 0 ]]
    then
        ## filter report.tsv, airr,tsv and toassemble.fa
        awk 'ARGIND==1{a[$1]}ARGIND==2{if(FNR==1){print}else if($1 in a){print}}' $gex_barcodes $report_input > $report_gexFiltered_out
        awk 'ARGIND==1{a[$1]}ARGIND==2{if($1~/^>/){readID=substr($1,2); getline; CB=$1; if(CB in a){print ">"readID"\n"CB}}}' $gex_barcodes $toassembly_input > $toassembly_gexFiltered_out
        filter_airr_by_report $report_gexFiltered_out $airr_input $airr_gexFiltered_out
        ## samples with GEX library will perform no additional filtration
        ln -s $report_gexFiltered_out $report_finalFiltered_out
    else
        ## empty gex_barcodes indicates no GEX dataset
        ln -s $report_input $report_gexFiltered_out
        ln -s $toassembly_input $toassembly_gexFiltered_out
        ln -s $airr_input $airr_gexFiltered_out
        ## Filter report and airr by multiple criterias
        filter_trust4_result.R -r $report_gexFiltered_out -a $airr_gexFiltered_out -o $report_finalFiltered_out
    fi

    ## Collecting metrics
    trust4_metrics.sh $report_finalFiltered_out \
                      $toassembly_gexFiltered_out \
                      $vdj_type \
                      $starsolo_summary \
                      ${vdj_type}.vdj_metrics.json \
                      ${vdj_type}.knee_input.tsv \
                      ${vdj_type}.cloneType_out.tsv

}

## https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash
IFS=',' read -ra ADDR <<< "$featureType"
IFS=',' read -ra SUADDR <<< "$starsoloSummary"

for i in "${!ADDR[@]}"
do
    vdj_type=${ADDR[i]}
    starsoloSummary=${SUADDR[i]}
    report_input="TRUST_${sampleID}_${vdj_type}_barcode_report.tsv"
    airr_input="TRUST_${sampleID}_${vdj_type}_barcode_airr.tsv"
    toassembly_input="TRUST_${sampleID}_${vdj_type}_toassemble_bc.fa"

    filter_VDJ $gex_cells $vdj_type $report_input $airr_input $toassembly_input $starsoloSummary
    ## rename output
    mv ${vdj_type}.vdj_metrics.json ${sampleID}_${vdj_type}.vdj_metrics.json
    mv ${vdj_type}.knee_input.tsv ${sampleID}_${vdj_type}.knee_input.tsv
    mv ${vdj_type}.cloneType_out.tsv ${sampleID}_${vdj_type}.cloneType_out.tsv
    mv ${report_input}.finalFiltered TRUST_${sampleID}_${vdj_type}_barcode_report.filtered.tsv
    filter_airr_by_report TRUST_${sampleID}_${vdj_type}_barcode_report.filtered.tsv $airr_input TRUST_${sampleID}_${vdj_type}_barcode_airr.filtered.tsv
done

