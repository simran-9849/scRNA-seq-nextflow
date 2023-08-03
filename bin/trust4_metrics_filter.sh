#! /usr/bin/env bash

gex_cells=$1
report_input=$2
airr_input=$3
toassemble_bc_input=$4
starsoloSummary=$5
vdj_type=$6

report_out=$(basename $report_input)
report_out=${report_out%%.tsv}".filtered.tsv"

airr_out=$(basename $airr_input)
airr_out=${airr_out%%.tsv}".filtered.tsv"

toassemble_bc_out=$(basename $toassemble_bc_input)
toassemble_bc_out=${toassemble_bc_out%%.fa}".filtered.fa"

awk 'ARGIND==1{a[$1]}ARGIND==2{if(FNR==1){print}else if($1 in a){print}}' $gex_cells $report_input > $report_out
awk 'ARGIND==1{a[$1]}ARGIND==2{if($1~/^>/){readID=substr($1,2); getline; CB=$1; if(CB in a){print ">"readID"\n"CB}}}' $gex_cells $toassemble_bc_input > $toassemble_bc_out
awk 'ARGIND==1{a[$1]}ARGIND==2{if(FNR==1){print}else if($1 in a){print}}' $gex_cells $airr_input > $airr_out

trust4_metrics.sh $report_out \
                  $toassemble_bc_out \
                  $vdj_type \
                  $starsoloSummary \
                  ${vdj_type}.vdj_metrics.json \
                  ${vdj_type}.knee_input.tsv \
                  ${vdj_type}.cloneType_out.tsv
