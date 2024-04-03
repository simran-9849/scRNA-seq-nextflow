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
cellListInput=$2
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

logTime(){
    local info=$1
    local currentTime=$(date +"%F %T")
    >&2 printf "$info at %s %s\n" $currentTime
}

## Process the BAM file and generate readInfo file
logTime "Starting processing bam file"

readInfo=$(mktemp -p ./)
if [[ $whiteList != "None" ]] && [[ -f $whiteList ]]
then
    samtools view -@ $thread $inputBAM |
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
    ' $whiteList - > $readInfo
else
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
        if((CB !="-") && (gene !="-")){
            print readID"\t"CB"\t"UMI"\t"gene
        }
    }
    ' > $readInfo
fi

logTime "Processing bam file ended"

prepare_downSample_read(){
    awk '
    BEGIN{
        srand()
        percentage[1]=0.05
        percentage[2]=0.1
        percentage[3]=0.15
        percentage[4]=0.2
        percentage[5]=0.25
        percentage[6]=0.3
        percentage[7]=0.4
        percentage[8]=0.6
        percentage[9]=0.8
        for(i=1;i<=length(percentage);i++){
            downFile[i]="sub_"percentage[i]"_readInfo"
        }
    }
    {
        randNumber=rand()
        for(i=1;i<=length(percentage);i++){
            if(randNumber<=percentage[i]){
                print > downFile[i]
            }
        }
    }
    ' $readInfo
    ln -s $readInfo sub_1_readInfo
}

calc_saturation(){
    ## percentage is  0.0 ≤ FLOAT ≤ 1.0
    local cellList=$1
    local readFile=$2
    local percentage=${readFile##sub_}
    percentage=${percentage%%_readInfo}
    sort -k 2,2 -k 3,3 -k 4,4 -k 1,1 -u -S 512M -T ./ $readFile |
        bedtools groupby -g 2,3 -c 1,4 -o distinct |
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
        ' $cellList -
}

## Start time
logTime "Start preparing downsample files"

prepare_downSample_read
## remove the original one to save space

logTime "Downsample files ready"

## Start time
logTime "Starting calculating"
export -f calc_saturation
ls sub_*_readInfo | awk -v cellListInput=$cellListInput '{print "calc_saturation "cellListInput" "$1}' |
    parallel -P $thread |
    jq --raw-input --slurp 'split("\n") |map(split("\t")) | .[0:-1] | map( { "percentage": .[0], "meanReadPerCell": .[1], "medianGenePerCell": .[2], "reads": .[3], "saturation": .[4] } )' > $outputJSON

## Get UMI hist for preseqR
## to invoke parallel in sort, avoiding pipe
tmpReadFile1=$(mktemp -p ./)
tmpReadFile2=$(mktemp -p ./)
cut -f 2,3 $readInfo > $tmpReadFile1
sort --parallel $thread -S 512M -T ./ $tmpReadFile1 > $tmpReadFile2
uniq -c $tmpReadFile2 |
    awk '{print $1}' > $tmpReadFile1
sort --parallel $thread -S 512M -T ./ $tmpReadFile1 > $tmpReadFile2
uniq -c $tmpReadFile2 |
    awk '{print $2"\t"$1}' > $tmpReadFile1
sort --parallel $thread -S 512M -T ./ -k 1,1n $tmpReadFile1 > $preseqR_UMI_hist
rm $tmpReadFile1 $tmpReadFile2
## Get gene hist for preseqR
## totalGeneCount summarizes all genes from different barcodes
tmpReadFile1=$(mktemp -p ./)
tmpReadFile2=$(mktemp -p ./)
awk 'ARGIND==1{cell[$1]}ARGIND==2&&$3!="-"{if($2 in cell){print $4"\t"$2}}' $cellListInput $readInfo > $tmpReadFile1
sort --parallel $thread -S 512M -T ./ $tmpReadFile1 |
    uniq -c > $tmpReadFile2
wc -l $tmpReadFile2 | awk '{print $1}' > $totalGeneCount
awk '{print $1}' $tmpReadFile2 > $tmpReadFile1
sort --parallel $thread -S 512M -T ./ $tmpReadFile1 > $tmpReadFile2
uniq -c $tmpReadFile2 |
    awk '{print $2"\t"$1}' > $tmpReadFile1
sort --parallel $thread -S 512M -T ./ -k 1,1n $tmpReadFile1 > $preseqR_gene_hist

rm $tmpReadFile1
rm $tmpReadFile2

rm $readInfo
## End time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation ended at %s %s\n" $currentTime
