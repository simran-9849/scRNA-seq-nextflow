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
downSampleCount=50000
CBtag="CB" ## tag used to extract CB from BAM
UMItag="UB" ## tag used to extract UMI from BAM, if UMItag="None", then reads will be used directly
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
	  --CBtag CB \\
	  --UMItag UB/None \\
	  --kneeOut kneeOut.tsv \\
	  --cellOut cellOut.tsv \\
	  --readIDout readID.lst \\
	  --barcode_fasta barcode.fa \\
	  --umi_fasta umi.fa \\
	  --threads 10 \\
	  --downSample 50000

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
        "--CBtag")
            CBtag=$2
            shift 2
            ;;
        "--UMItag")
            UMItag=$2
            shift 2
            ;;
        "--downSample")
            downSampleCount=$2
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
    umi_cutoff_rank=$(awk -v expectedCells=$expectedCells -v percentile=$percentile 'BEGIN{printf "%d\n", expectedCells*(1-percentile)}')
    ##echo "Using umi_cutoff_rank: $umi_cutoff_rank"

    CB_reads=$(mktemp -p ./)
    samtools view -@ $threads $inputBAM |
        awk -v CBtag=$CBtag -v UMItag=$UMItag '
        BEGIN{
            CBpattern=CBtag":Z:"
            UMIpattern=UMItag":Z:"
        }
        {
            for(i=12;i<=NF;i++){
                if($i ~ CBpattern){
                    split($i, tmp, ":");
                    CB=tmp[3]
                }
                if($i ~ UMIpattern){
                    split($i, tmp, ":");
                    UB=tmp[3]
                }
            }
            if(UMItag=="None" && CB!="-"){
                print CB"\t"$1
            }else if(CB!="-" && UB!="-" && UMItag!="None"){
                print CB"\t"UB"\t"$1
            }
        }
        ' |
        sort -u -T ./ --parallel $threads |
        sort -k 1,1 -k 2,2 -T ./ --parallel $threads > $CB_reads

    if [[ $UMItag=="None" ]]
    then
        bedtools groupby -g 1 -c 2 -o count -i $CB_reads |
            sort -k 2,2rn -T ./ --parallel $threads > $kneeDataOut
    else
        cut -f 1,2 $CB_reads |
            sort -u -T ./ --parallel $threads |
            sort -k1,1 -T ./ --parallel $threads |
            bedtools groupby -g 1 -c 2 -o count |
            sort -k 2,2rn -T ./ --parallel $threads > $kneeDataOut
    fi
    awk -v umi_cutoff_rank=$umi_cutoff_rank -v umi_fold=$umi_fold '
    ARGIND==1 && FNR==umi_cutoff_rank{
        umi_cutoff=$2/umi_fold
    }
    ARGIND==2{
        if($2 >= umi_cutoff){
            print $1
        }
    }
    ' $kneeDataOut $kneeDataOut > $cellListOut
else
    CB_reads=$(mktemp -p ./)
    samtools view -@ $threads -D ${CBtag}:${gexBarcode} $inputBAM |
        awk -v CBtag=$CBtag -v UMItag=$UMItag '
        BEGIN{
            CBpattern=CBtag":Z:"
            UMIpattern=UMItag":Z:"
        }
        {
            for(i=12;i<=NF;i++){
                if($i ~ CBpattern){
                    split($i, tmp, ":");
                    CB=tmp[3]
                }
                if($i ~ UMIpattern){
                    split($i, tmp, ":");
                    UB=tmp[3]
                }
            }
            if(UMItag=="None" && CB!="-"){
                print CB"\t"$1
            }else if(CB!="-" && UB!="-" && UMItag!="None"){
                print CB"\t"UB"\t"$1
            }
        }
        ' |
        sort -u -T ./ --parallel $threads |
        sort -k 1,1 -k 2,2 -T ./ --parallel $threads > $CB_reads

    if [[ $UMItag=="None" ]]
    then
        bedtools groupby -g 1 -c 2 -o count -i $CB_reads |
            sort -k 2,2rn -T ./ --parallel $threads > $kneeDataOut
    else
        cut -f 1,2 $CB_reads |
            sort -u -T ./ --parallel $threads |
            sort -k1,1 -T ./ --parallel $threads |
            bedtools groupby -g 1 -c 2 -o count |
            sort -k 2,2rn -T ./ --parallel $threads > $kneeDataOut
    fi
    cat $gexBarcode > $cellListOut
fi

## Extract cell reads ID
## Downsample if there are more than 50,000 reads for a cell
readID_down=$(mktemp -p ./)
awk -v downSampleCount=$downSampleCount '
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
        if(($1 in bc) && umi_idx<=downSampleCount){
            print
        }
    }
    ' $cellListOut $CB_reads > $readID_down

if [[ $UMItag != "None" ]]
then
    awk '{print $3}' $readID_down > $readIDList
    awk '{print ">"$3"\n"$1}' $readID_down > $barcode_fasta
    awk '{print ">"$3"\n"$2}' $readID_down > $umi_fasta
else
    awk '{print $2}' $readID_down > $readIDList
    awk '{print ">"$2"\n"$1}' $readID_down > $barcode_fasta
    touch $umi_fasta
fi

rm $CB_reads $readID_down
