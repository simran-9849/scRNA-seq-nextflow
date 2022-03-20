#! /usr/bin/env bash

## This script is to generate data for
## sequencing saturation curve

## The sequencing saturation inrespects cells (consistent with STARsolo)
## The median genes per cell considers only cells

## Input is STARsolo mapping BAM file
## dependency: samtools, bedtools, gawk, sort

set -euo pipefail

usage="Usage: $(basename $0) whitelist cellList inputBAM thread outputJSON

Please set whitelist as None when it's not available.
"

if [[ $# -eq 0 ]] || [[ $1 == "-h" ]] || [[ $1 == "-help" ]] || [[ $1 == "--help" ]]
then
    echo "$usage"
    exit 0
fi

whitelist=$1
cellList=$2
inputBAM=$3
thread=$4
outputJSON=$5

if [[ ! -f $whitelist ]] && [[ $whitelist != "None" ]]
then
    echo "$whitelist not found."
    exit 1
fi

if [[ ! -f $inputBAM ]]
then
    echo "$inputBAM not found."
    exit 1
fi

function calc_saturation {
    local readInfo=$(mktemp -p ./)
    ## percentage is  0.0 ≤ FLOAT ≤ 1.0
    local percentage=$1

    if [[ $whitelist != "None" ]]
    then
        samtools view -@ $thread -u -D CB:$whitelist $inputBAM |
            samtools view -@ $thread -d GX --subsample $percentage --subsample-seed 1324 |
            awk '
            {
                readID=$1;
                for(i=12;i<=NF;i++){
                    if($i~/^CB/){split($i, tmp, ":");
                        CB=tmp[3]
                    }else if($i~/^UB/){
                        split($i, tmp, ":");
                        UMI=tmp[3]
                    }else if($i~/^GX/){
                        split($i, tmp, ":");
                        gene=tmp[3]
                    }
                };
                print readID"\t"CB"\t"UMI"\t"gene
            }
            ' |
            sort -k 2,2 -k 3,3 -k 4,4 -k 1,1 -u --parallel $thread -S 10% -T ./ > $readInfo
        cat $readInfo | bedtools groupby -g 2,3 -c 1,4 -o distinct |
            awk '
            BEGIN{
                m=0;n=0
            }
            NR==FNR{
                cell[$1]
            }
            NR>FNR{
                if($2!="-"){m+=1};
                split($3, readArray, ",");
                n+=length(readArray);
                if($1 in cell){
                    readCount[$1]+=length(readArray)
                    if($2!="-"){barocodeGene[$1","$4]}
                }
            }
            END{
                asort(readCount)
                totalRead=0
                for(i=1;i<=length(readCount);i++){
                    totalRead+=readCount[i]
                }
                meanRead=totalRead/length(readCount)
                for(i in barocodeGene){
                    split(i, tmp, ",")
                    geneCount[tmp[1]]+=1
                }
                asort(geneCount)
                if(length(geneCount)%2==1){
                    medianGene=geneCount[(length(geneCount)+1)/2]
                }else{
                    medianGene=(geneCount[length(geneCount)/2] + geneCount[length(geneCount)/2+1])/2
                }
                print "'$percentage'\t"meanRead"\t"medianGene"\t"n"\t"1-m/n
            }
            ' $cellList -
    else
        samtools view -@ $thread -d GX --subsample $percentage --subsample-seed 1324 $inputBAM |
            awk '
            {
                readID=$1;
                for(i=12;i<=NF;i++){
                    if($i~/^CB/){split($i, tmp, ":");
                        CB=tmp[3]
                    }else if($i~/^UB/){
                        split($i, tmp, ":");
                        UMI=tmp[3]
                    }else if($i~/^GX/){
                        split($i, tmp, ":");
                        gene=tmp[3]
                    }
                };
                print readID"\t"CB"\t"UMI"\t"gene
            }
            ' |
            sort -k 2,2 -k 3,3 -k 4,4 -k 1,1 -u --parallel $thread -S 10% -T ./ > $readInfo
        cat $readInfo | bedtools groupby -g 2,3 -c 1,4 -o distinct |
            awk '
            BEGIN{
                m=0;n=0
            }
            NR==FNR{
                cell[$1]
            }
            NR>FNR{
                if($2!="-"){m+=1};
                split($3, readArray, ",");
                n+=length(readArray);
                if($1 in cell){
                    readCount[$1]+=length(readArray)
                    if($2!="-"){barocodeGene[$1","$4]}
                }
            }
            END{
                asort(readCount)
                totalRead=0
                for(i=1;i<=length(readCount);i++){
                    totalRead+=readCount[i]
                }
                meanRead=totalRead/length(readCount)
                for(i in barocodeGene){
                    split(i, tmp, ",")
                    geneCount[tmp[1]]+=1
                }
                asort(geneCount)
                if(length(geneCount)%2==1){
                    medianGene=geneCount[(length(geneCount)+1)/2]
                }else{
                    medianGene=(geneCount[length(geneCount)/2] + geneCount[length(geneCount)/2+1])/2
                }
                print "'$percentage'\t"meanRead"\t"medianGene"\t"n"\t"1-m/n
            }
            ' $cellList -
    fi
    rm $readInfo
}

## Start time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation started at %s %s\n" $currentTime

for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1};
do
    calc_saturation $i
done |
    jq --raw-input --slurp 'split("\n") |map(split("\t")) | .[0:-1] | map( { "percentage": .[0], "meanReadPerCell": .[1], "medianGenePercell": .[2], "reads": .[3], "saturation": .[4] } )' > $outputJSON

## End time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation ended at %s %s\n" $currentTime
