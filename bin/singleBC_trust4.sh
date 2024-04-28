#! /usr/bin/env bash

set -eo pipefail

usage(){
    cat <<-EOF
Basic Usage:
==========
$(basename $0) --genomeRef genome_bcr.fa \\
               --imgtRef imgt_bcr.fa \\
               --barcode VDJ-B_barcode.fa \\
               --UMI VDJ-B_umi.fa \\
               --inputR1 trust4_input.R1.fq.gz \\
               --inputR2 trust4_input.R2.fq.gz \\
               --readFormat r1:69:-1 \\
               --threads 8 \\
               --reportOut  combined_barcode_report.tsv \\
               --airrOut combined_barcode_airr.tsv \\
               --readsAssign combined_reads_assign.tsv \\
               --annotFasta combined_annot.fa \\
               --finalOut combined_final.out

Dependency:
  fastx_toolkit, seqtk, trust4, parallel
EOF
}

if [[ -z $1 ]]
then
    usage
    exit 0
fi

genomeRef=""
imgtRef=""
inputBC=""
inputUMI=""
inputR1=""
inputR2=""
readFormat="r1:69:-1"
threads=8
reportOut=""
airrOut=""
readsAssign=""
annotFa=""
finalOut=""
R2_Only=false

while [[ ! -z $1 ]]
do
    case $1 in
        "" | "-h" | "-help" | "--help")
            usage
            exit 0
            ;;
        "--genomeRef")
            genomeRef=$2
            shift 2
            ;;
        "--imgtRef")
            imgtRef=$2
            shift 2
            ;;
        "--barcode")
            inputBC=$2
            shift 2
            ;;
        "--UMI")
            inputUMI=$2
            shift 2
            ;;
        "--inputR1")
            inputR1=$2
            shift 2
            ;;
        "--inputR2")
            inputR2=$2
            shift 2
            ;;
        "--readFormat")
            readFormat=$2
            shift 2
            ;;
        "--threads")
            threads=$2
            shift 2
            ;;
        "--reportOut")
            reportOut=$2
            shift 2
            ;;
        "--airrOut")
            airrOut=$2
            shift 2
            ;;
        "--readsAssign")
            readsAssign=$2
            shift 2
            ;;
        "--R2_Only")
            R2_Only=true
            shift 1
            ;;
        "--annotFasta")
            annotFa=$2
            shift 2
            ;;
        "--finalOut")
            finalOut=$2
            shift 2
            ;;
        *)
            echo "$1 is not a valid argument"
            exit 1
    esac
done

bc_tsv=$(mktemp -p ./)
fasta_formatter -t -i $inputBC | cut -f 2 | sort -u > $bc_tsv

awk 'BEGIN{print "#barcode\tcell_type\tchain1\tchain2\tsecondary_chain1\tsecondary_chain2"}' > $reportOut
awk 'BEGIN{print "sequence_id\tsequence\trev_comp\tproductive\tv_call\td_call\tj_call\tc_call\tsequence_alignment\tgermline_alignment\tcdr1\tcdr2\tjunction\tjunction_aa\tv_cigar\td_cigar\tj_cigar\tv_identity\tj_identity\tcell_id\tcomplete_vdj\tconsensus_count"}' > $airrOut
if [[ -f $readsAssign ]]
then
    rm $readsAssign
fi
touch $readsAssign
if [[ -f $annotFa ]]
then
    rm $annotFa
fi
touch $annotFa
if [[ -f $finalOut ]]
then
    rm $finalOut
fi
touch $finalOut


split_bc(){

    local wd=$1
    local barcodeFa=$2
    local umiFa=$3

    mkdir -p "$wd/barcode"
    mkdir -p "$wd/readsList"
    mkdir -p "$wd/umi"

    ## combine barcode, umi, readID
    local bc_umi_reads=$(mktemp -p ./)
    local umiCounts=$(grep -c ">" $umiFa)
    if [[ $umiCounts -gt 0 ]]
    then
        join -j 1 -t $'\t' \
             <(fasta_formatter -t -i $barcodeFa | sort -k 1,1 -k 2,2 --parallel $threads) \
             <(fasta_formatter -t -i $umiFa | sort -k 1,1 -k 2,2 --parallel $threads) > $bc_umi_reads
        sort -k 2,2 -k 3,3 --parallel $threads $bc_umi_reads |
            awk -v wd="$wd" '
            {
                bc_out=wd"/barcode/"$2".single_bc.fa";
                umi_out=wd"/umi/"$2".single_umi.fa";
                reads_out=wd"/readsList/"$2".single_reads.lst";
                print ">"$1"\n"$2 > bc_out;
                print ">"$1"\n"$3 > umi_out;
                print $1 > reads_out;
            }
            '
        rm $bc_umi_reads
    else
        fasta_formatter -t -i $barcodeFa |
            awk -v wd="$wd" '
            {
                bc_out=wd"/barcode/"$2".single_bc.fa";
                reads_out=wd"/readsList/"$2".single_reads.lst";
                print ">"$1"\n"$2 > bc_out;
                print $1 > reads_out;
            }
            '
        fasta_formatter -t -i $barcodeFa | cut -f 2 | sort -u --parallel $threads |
            parallel -j $threads touch $wd/umi/{}".single_umi.fa"
    fi
}

singleBC_trust4(){
    ## function cannot access global variable when parsing to GNU parallel
    local bc=$1
    local readsList=$2
    local barcodeFa=$3
    local umiFa=$4
    local trust4_wd=$5

    mkdir -p $trust4_wd
    ## extract read list
    ##readID=$(mktemp -p ./ readID.XXXXXX)
    ##fasta_formatter -t -i $inputBC | awk -v bc="$bc" '$2==bc{print $1}' > $readID
    ##singleBC=$(mktemp -p ./ BC.XXXXXX.fa)
    ##fasta_formatter -t -i $inputBC | awk -v bc="$bc" '$2==bc{print ">"$1"\n"$NF}' > $singleBC
    local umiCounts=$(grep -c ">" $umiFa)
    if [[ $umiCounts -eq 0 ]]
    then
        local barcode_umi_opt="--barcode $barcodeFa"
    else
        ##singleUMI=$(mktemp -p ./ UMI.XXXXXX.fa)
        ##fasta_formatter -t -i $inputUMI | awk 'ARGIND==1{a[$1]}ARGIND==2{if($1 in a){print ">"$1"\n"$NF}}' $readID -  > $singleUMI
        local barcode_umi_opt="--barcode $barcodeFa --UMI $umiFa"
    fi

    if [[ $R2_Only == true ]]
    then
        local tmpR2=$(mktemp -p "$trust4_wd" R2.XXXXXX.fq)
        seqtk subseq $inputR2 $readsList > $tmpR2
        local trust4_input_opt="-u $tmpR2"
    else
        local tmpR1=$(mktemp -p "$trust4_wd" R1.XXXXXX.fq)
        local tmpR2=$(mktemp -p "$trust4_wd" R2.XXXXXX.fq)
        seqtk subseq $inputR1 $readsList > $tmpR1
        seqtk subseq $inputR2 $readsList > $tmpR2
        local trust4_input_opt="-1 $tmpR1 -2 $tmpR2"
    fi
    local trust4_out=$(mktemp -d -p "$trust4_wd" "${bc}.XXXXXX")
    ## separate trust4 steps, since fastq-extractor will perform filtration
    ## and there may be no reads left
    fastq-extractor -t 1 \
                    -f $genomeRef \
                    -o $trust4_out/singleBC_toassemble \
                    --readFormat $readFormat \
                    $trust4_input_opt \
                    $barcode_umi_opt
    ## check reads number
    local readsNum=$(grep -c ">" $trust4_out/singleBC_toassemble_bc.fa)
    echo $readsNum
    if [[ $readsNum -gt 0 ]]
    then
        run-trust4 -f $genomeRef \
                   --ref $imgtRef \
                   -o singleBC \
                   --od $trust4_out \
                   $trust4_input_opt \
                   $barcode_umi_opt \
                   --readFormat $readFormat \
                   --outputReadAssignment \
                   -t 1 \
                   --stage 1

        ##awk 'NR>1{print}' $trust4_out/singleBC_barcode_report.tsv >> $reportOut
        ##awk 'NR>1{print}' $trust4_out/singleBC_barcode_airr.tsv >> $airrOut
        ##cat $trust4_out/singleBC_assign.out >> $readsAssign
        ##cat $trust4_out/singleBC_annot.fa >> $annotFa
        ##cat $trust4_out/singleBC_final.out >> $finalOut
    fi

    if [[ $R2_Only == true ]]
    then
        rm $tmpR2
    else
        rm $tmpR1 $tmpR2
    fi
    ##rm -rf $trust4_out
}


## split barcode
workDir=$(mktemp -d -p ./)
split_bc $workDir $inputBC $inputUMI

export -f singleBC_trust4
export inputR1
export inputR2
export genomeRef
export imgtRef
export R2_Only
export readFormat
##export reportOut
##export airrOut
##export readsAssign
##export annotFa
##export finalOut

cat $bc_tsv |
    parallel -j $threads singleBC_trust4 {} $workDir/readsList/{}.single_reads.lst $workDir/barcode/{}.single_bc.fa $workDir/umi/{}.single_umi.fa $workDir

awk 'FNR>1{print}' $workDir/*/singleBC_barcode_report.tsv >> $reportOut
awk 'FNR>1{print}' $workDir/*/singleBC_barcode_airr.tsv >> $airrOut
cat $workDir/*/singleBC_assign.out >> $readsAssign
cat $workDir/*/singleBC_annot.fa >> $annotFa
cat $workDir/*/singleBC_final.out >> $finalOut

##readarray -t Arr < $bc_tsv
##for i in "${Arr[@]}"
##do
##    singleBC_trust4 $i
##done
rm $bc_tsv
## Remove trust4 temp workDir
rm $workDir
