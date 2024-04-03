#! /usr/bin/env bash

## This script is to generate data for
## sequencing saturation curve

## The sequencing saturation inrespects cells (consistent with STARsolo)
## The median genes per cell considers only cells

## Input is STARsolo mapping BAM file
## dependency: samtools, bedtools, gawk, sort, shuf

set -euo pipefail

usage(){

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

## Inputs
whiteList=$1
cellList=$2
multiMapper=$3
inputBAM=$4
thread=$5
## Outputs
outputJSON=$6
preseqR_UMI_hist=$7
preseqR_gene_hist=$8
totalGeneCount=$9

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



## Process the BAM file and generate readInfo file
currentTime=$(date +"%F %T")
>&2 printf "Starting processing bam file at %s %s\n" $currentTime

tmpReadList=$(mktemp -p ./)
samtools view -@ $thread $inputBAM |
awk '
{
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
    print readID"\t"CB"\t"UMI"\t"gene
}
' > $tmpReadList

readInfo=$(mktemp -p ./)
if [[ $whiteList != "None" ]] && [[ -f $whiteList ]]
then
    awk '
    NR==FNR{
      wl[$1]
    }
    NR>FNR{
      CB=$2
      gene=$4
      if((CB in wl) && (gene !="-")){
        print
      }
    }
    ' $whiteList $tmpReadList > $readInfo
else
    awk '
    {
      CB=$2
      gene=$4
      if((CB !="-") && (gene !="-")){
        print
      }
    }
    ' $tmpReadList > $readInfo
fi
rm $tmpReadList

currentTime=$(date +"%F %T")
>&2 printf "Processing bam file ended at %s %s\n" $currentTime

calc_saturation(){
    ## percentage is  0.0 ≤ FLOAT ≤ 1.0
    local percentage=$1
    local nTotal=$(wc -l $readInfo | awk '{print $1}')
    local nRow=$(awk -v nTotal=$nTotal -v frac=$percentage 'BEGIN{printf "%i\n", nTotal*frac}')
    tmpDownFile=$(mktemp -p ./)
    ## use sort -R instead of shuf to save memory
    awk -v cutoff=$percentage 'BEGIN{srand()}{if(rand()<=cutoff){print}}' $readInfo  > $tmpDownFile
    tmpSortedDownFile=$(mktemp -p ./)
    sort -k 2,2 -k 3,3 -k 4,4 -k 1,1 -u --parallel -S 512M $thread -T ./ $tmpDownFile > $tmpSortedDownFile
    rm $tmpDownFile
    bedtools groupby -g 2,3 -c 1,4 -o distinct -i $tmpSortedDownFile |
        awk -v percentage=$percentage '
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
            print percentage"\t"meanRead"\t"medianGene"\t"n"\t"1-m/n
        }
        ' $cellList - && rm $tmpSortedDownFile
}

## Start time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation started at %s %s\n" $currentTime


for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1};
##for i in {0.05,1}
do
    printf "Subsampling: %s\n" $i >&2
    calc_saturation $i
    printf "calculation finished...\n" >&2
done |
    jq --raw-input --slurp 'split("\n") |map(split("\t")) | .[0:-1] | map( { "percentage": .[0], "meanReadPerCell": .[1], "medianGenePerCell": .[2], "reads": .[3], "saturation": .[4] } )' > $outputJSON

## Get UMI hist for preseqR
## to invoke parallel in sort, avoiding pipe
tmpReadFile1=$(mktemp -p ./)
tmpReadFile2=$(mktemp -p ./)
cut -f 2,3 $readInfo > $tmpReadFile1
sort --parallel $thread -T ./ $tmpReadFile1 > $tmpReadFile2
uniq -c $tmpReadFile2 |
    awk '{print $1}' > $tmpReadFile1
sort --parallel $thread -T ./ $tmpReadFile1 > $tmpReadFile2
uniq -c $tmpReadFile2 |
    awk '{print $2"\t"$1}' > $tmpReadFile1
sort --parallel $thread -T ./ -k 1,1n $tmpReadFile1 > $preseqR_UMI_hist
rm $tmpReadFile1 $tmpReadFile2
## Get gene hist for preseqR
## totalGeneCount summarizes all genes from different barcodes
tmpReadFile1=$(mktemp -p ./)
tmpReadFile2=$(mktemp -p ./)
awk 'ARGIND==1{cell[$1]}ARGIND==2&&$3!="-"{if($2 in cell){print $4"\t"$2}}' $cellList $readInfo > $tmpReadFile1
sort --parallel $thread -T ./ $tmpReadFile1 |
    uniq -c > $tmpReadFile2
wc -l $tmpReadFile2 | awk '{print $1}' > $totalGeneCount
awk '{print $1}' $tmpReadFile2 > $tmpReadFile1
sort --parallel $thread -T ./ $tmpReadFile1 > $tmpReadFile2
uniq -c $tmpReadFile2 |
    awk '{print $2"\t"$1}' > $tmpReadFile1
sort --parallel $thread -T ./ -k 1,1n $tmpReadFile1 > $preseqR_gene_hist

rm $tmpReadFile1
rm $tmpReadFile2

rm $readInfo
## End time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation ended at %s %s\n" $currentTime
