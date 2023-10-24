#! /usr/bin/env bash

## This script is to extract Barcode and UMI information
## from starsolo BAM output, and perform traditional
## cell calling (e.g. 3000, 0.99, 10) to define cells

## Dependency: samtools, bedtools

set -eo pipefail

## Default values
inputBAM=""
gexBarcode="None"
expectedCells=3000
percentile=0.95
umi_fold=10
threads=10
kneeDataOut="kneeOut.tsv"
cellListOut="cellOut.tsv"
readIDList="readID.lst"
barcode_fasta="bc.fa"
umi_fasta="umi.fa"

usage(){
    cat <<-EOF
	Basic Usage:
	==========
	$(basename $0) --inputBAM bamFile \\
	  --gexBarcode barcodeFile.tsv/None \\
	  --expectedCells 3000 \\
	  --percentile 0.95 \\
	  --umi_fold 10 \\
	  --kneeOut kneeOut.tsv \\
	  --cellOut cellOut.tsv \\
	  --readIDout readID.lst \\
	  --barcode_fasta barcode.fa \\
	  --umi_fasta umi.fa \\
	  --threads 10

	Dependency:
	  samtools, bedtools
	EOF
}

if [[ -z $1 ]]
then
    usage
    exit 0
fi

while [[ ! -z $1 ]]
do
    case $1 in
        "" | "-h" | "-help" | "--help")
            usage
            exit 0
            ;;
        "--inputBAM")
            inputBAM=$2
            shift 2
            ;;
        "--gexBarcode")
            gexBarcode=$2
            shift 2
            ;;
        "--expectedCells")
            expectedCells=$2
            shift 2
            ;;
        "--percentile")
            percentile=$2
            shift 2
            ;;
        "--umi_fold")
            umi_fold=$2
            shift 2
            ;;
        "--kneeOut")
            kneeDataOut=$2
            shift 2
            ;;
        "--cellOut")
            cellListOut=$2
            shift 2
            ;;
        "--barcode_fasta")
            barcode_fasta=$2
            shift 2
            ;;
        "--umi_fasta")
            umi_fasta=$2
            shift 2
            ;;
        "--readIDout")
            readIDList=$2
            shift 2
            ;;
        "--threads")
            threads=$2
            shift 2
            ;;
        *)
            echo "$1 is not a valid argument"
            exit 1
    esac
done

## If gexBarcode is None, use traditional cell calling method
if [[  $gexBarcode == "None" ]]
then
    ## calculate cutoff rank
    umi_cutoff_rank=$(awk -v expectedCells=$expectedCells -v percentil=$percentile 'BEGIN{printf "%d\n", expectedCells*(1-percentile)}')

    CB_reads=$(mktemp -p ./)
    samtools view -@ $threads $inputBAM |
        awk '
    {
        for(i=12;i<=NF;i++){
            if($i~/CB:Z:/){
                split($i, tmp, ":");
                CB=tmp[3]
            }
            if($i~/UB:Z:/){
                split($i, tmp, ":");
                UB=tmp[3]
            }
        }
        if(CB!="-" && UB!="-"){
            print CB"\t"UB"\t"$1
        }
    }
    ' |
        sort -u -T ./ --parallel $threads |
        sort -k 1,1 -k 2,2 -T ./ --parallel $threads > $CB_reads

    cut -f 1,2 $CB_reads |
        sort -u -T ./ --parallel $threads |
        sort -k1,1 -T ./ --parallel $threads |
        bedtools groupby -g 1 -c 2 -o count |
        sort -k 2,2rn -T ./ --parallel $threads > $kneeDataOut

    awk -v umi_cutoff_rank=$umi_cutoff_rank 'NR<=umi_cutoff_rank{print $1}' $kneeDataOut > $cellListOut
else
    CB_reads=$(mktemp -p ./)
    samtools view -@ $threads -D CB:$gexBarcode $inputBAM |
        awk '
        {
            for(i=12;i<=NF;i++){
                if($i~/CB:Z:/){
                    split($i, tmp, ":");
                    CB=tmp[3]
                }
                if($i~/UB:Z:/){
                    split($i, tmp, ":");
                    UB=tmp[3]
                }
            }
            if(CB!="-" && UB!="-"){
                print CB"\t"UB"\t"$1
            }
        }
        ' |
        sort -u -T ./ --parallel $threads |
        sort -k 1,1 -k 2,2 -T ./ --parallel $threads > $CB_reads

    cut -f 1,2 $CB_reads |
        sort -u -T ./ --parallel $threads |
        sort -k1,1 -T ./ --parallel $threads |
        bedtools groupby -g 1 -c 2 -o count |
        sort -k 2,2rn -T ./ --parallel $threads > $kneeDataOut

    cat $gexBarcode > $cellListOut
fi

## Extract cell reads ID
## Downsample if there are more than 50,000 reads for a cell
readID_down=$(mktemp -p ./)
awk '
    ARGIND==1{
        bc[$1]
    }
    ARGIND==2{
        if($1!=CB){
            umi_idx=1
            CB=$1
            UMI=$2
        }else if($1==CB && $2!=UMI){
            umi_idx=umi_idx+1
            UMI=$2
        }
        if(($1 in bc) && umi_idx<=50000){
            print
        }
    }
    ' $cellListOut $CB_reads > $readID_down

awk '{print $3}' $readID_down > $readIDList
awk '{print ">"$3"\n"$1}' $readID_down > $barcode_fasta
awk '{print ">"$3"\n"$2}' $readID_down > $umi_fasta

rm $CB_reads $readID_down
