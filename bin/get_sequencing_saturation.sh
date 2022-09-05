#! /usr/bin/env bash

## This script is to generate data for
## sequencing saturation curve

## The sequencing saturation inrespects cells (consistent with STARsolo)
## The median genes per cell considers only cells

## Input is STARsolo mapping BAM file
## dependency: samtools, bedtools, gawk, sort

set -euo pipefail

function usage {

    cat<<-EOF
	Usage:
	=========
	$(basename $0) whiteList cellList multiMapper inputBAM thread outputJSON
	
	All inputs are positional.
	
	whiteList        Input whiteList file name. Please set
	                 whiteList as "None" when it's not available.
	                 
	cellList         Input cellList file, e.g. barcodes.tsv from
	                 STARsolo's output directory. The cellList doesn't
	                 support <(zcat barcodes.tsv.gz) as input.
	                 
	multiMapper      "unique" or "multiple", indicates calculating
	                 saturation results from unique gene reads or
	                 multiple gene reads.
	                 
	inputBAM         Input BAM file, i.e. STARsolo output BAM.
	                 
	thread           Number of threads to use.
	                 
	outputJSON       Output json file name
EOF
}

if [[ $# -eq 0 ]] || [[ $1 == "-h" ]] || [[ $1 == "-help" ]] || [[ $1 == "--help" ]]
then
    usage
    exit 0
fi

## default options
multiMapper="unique"
thread=4

whiteList=$1
cellList=$2
multiMapper=$3
inputBAM=$4
thread=$5
outputJSON=$6

if [[ ! -f $whiteList ]] && [[ $whiteList != "None" ]]
then
    >&2 echo "$whiteList not found."
    exit 1
fi

if [[ ! -f $inputBAM ]]
then
    echo "$inputBAM not found." >&2
    exit 1
fi

if [[ $multiMapper == "unique" ]]
then
    geneTag="GX"
elif [[ $multiMapper == "multiple" ]]
then
    geneTag="gx"
else
    echo "multiMapper option only supports unique or multiple" >&2
    exit 1
fi

readInfo=$(mktemp -p ./)

function processBAM {
    ## percentage is  0.0 ≤ FLOAT ≤ 1.0
    local percentage=$1
    if [[ $whiteList != "None" ]]
    then
        samtools view -@ $thread --subsample $percentage --subsample-seed 1324 $inputBAM |
            awk '
            NR==FNR{
                wl[$1]
            }
            NR>FNR{
                readID=$1;
                CB="-"
                UMI="-"
                gene="-"
                for(i=12;i<=NF;i++){
                    if($i~/^CB/){
                        CB=substr($i, 6)
                    }else if($i~/^UB/){
                        UMI=substr($i, 6)
                    }else if($i~/^'$geneTag'/){
                        gene=substr($i, 6)
                    }
                };
                if((CB in wl) && (gene !="-")){
                    print readID"\t"CB"\t"UMI"\t"gene
                }
            }
            ' $whiteList - |
            sort -k 2,2 -k 3,3 -k 4,4 -k 1,1 -u --parallel $thread -S 10% -T ./ > $readInfo
    else
        samtools view -@ $thread --subsample $percentage --subsample-seed 1324 $inputBAM |
            awk '
            {
                readID=$1;
                for(i=12;i<=NF;i++){
                    if($i~/^CB/){
                        CB=substr($i, 6)
                    }else if($i~/^UB/){
                        UMI=substr($i, 6)
                    }else if($i~/^'$geneTag'/){
                        gene=substr($i, 6)
                    }
                };
                if(CB!="-" && gene !="-"){
                    print readID"\t"CB"\t"UMI"\t"gene
                }
            }
            ' |
            sort -k 2,2 -k 3,3 -k 4,4 -k 1,1 -u --parallel $thread -S 10% -T ./ > $readInfo
    fi
}

function calc_saturation {
    ## percentage is  0.0 ≤ FLOAT ≤ 1.0
    local percentage=$1
    #local inputReadInfo=$2
    bedtools groupby -g 2,3 -c 1,4 -o distinct -i $readInfo |
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
}

## Start time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation started at %s %s\n" $currentTime

##for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1};
for i in {0.05,0.1,0.15}
do
    printf "Subsampling: %s\n" $i >&2
    processBAM $i
    printf "Processing BAM finished...\n" >&2
    calc_saturation $i
    printf "calculation finished...\n" >&2
done |
    jq --raw-input --slurp 'split("\n") |map(split("\t")) | .[0:-1] | map( { "percentage": .[0], "meanReadPerCell": .[1], "medianGenePerCell": .[2], "reads": .[3], "saturation": .[4] } )' > $outputJSON

rm $readInfo

## End time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation ended at %s %s\n" $currentTime
