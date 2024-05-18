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
    mkdir -p "$wd/readsFastq"

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
                bc_out=wd"/barcode/"$2".sc_bc.fa";
                umi_out=wd"/umi/"$2".sc_umi.fa";
                reads_out=wd"/readsList/"$2".sc_reads.lst";
                print ">"$1"\n"$2 > bc_out;
                print ">"$1"\n"$3 > umi_out;
                print $1 > reads_out;
            }
            '
    else
        fasta_formatter -t -i $barcodeFa | awk '{print $1"\t"$2"\t-"}' > $bc_umi_reads
        sort -k 2,2 -k 3,3 --parallel $threads $bc_umi_reads |
            awk -v wd="$wd" '
            {
                bc_out=wd"/barcode/"$2".sc_bc.fa";
                reads_out=wd"/readsList/"$2".sc_reads.lst";
                print ">"$1"\n"$2 > bc_out;
                print $1 > reads_out;
            }
            '
        fasta_formatter -t -i $barcodeFa | cut -f 2 | sort -u --parallel $threads |
            parallel -j $threads touch $wd/umi/{}".sc_umi.fa"
    fi
    ## prepare reads
    tmpR1_tsv=$(mktemp -p ./ R1.XXXXXX.tsv)
    tmpR1_sorted=$(mktemp -p ./ R1.XXXXXX.sorted.tsv)
    tmpR2_tsv=$(mktemp -p ./ R2.XXXXXX.tsv)
    tmpR2_sorted=$(mktemp -p ./ R2.XXXXXX.sorted.tsv)
    tmpbc_sorted=$(mktemp -p ./ bc.XXXXXX.sorted.tsv)
    sort -k 1,1 --parallel $threads $bc_umi_reads > $tmpbc_sorted
    rm $bc_umi_reads

    zcat $inputR1 | paste - - - - | awk '{print substr($1, 2)"\t"$0}' > $tmpR1_tsv
    sort -k 1,1 --parallel $threads $tmpR1_tsv > $tmpR1_sorted
    rm $tmpR1_tsv
    join -j 1 -t $'\t' $tmpbc_sorted $tmpR1_sorted |
        awk -F"\t" -v wd="$wd" '
        {
            reads_out=wd"/readsFastq/"$2".sc_reads.R1.fq"
            print $4"\n"$5"\n"$6"\n"$7 > reads_out
        }
        '

    zcat $inputR2 | paste - - - - | awk '{print substr($1, 2)"\t"$0}' > $tmpR2_tsv
    sort -k 1,1 --parallel $threads $tmpR2_tsv > $tmpR2_sorted
    rm $tmpR2_tsv
    join -j 1 -t $'\t' $tmpbc_sorted $tmpR2_sorted |
        awk -F"\t" -v wd="$wd" '
        {
            reads_out=wd"/readsFastq/"$2".sc_reads.R2.fq"
            print $4"\n"$5"\n"$6"\n"$7 > reads_out
        }
        '
    ## remove temp files
    rm $tmpR1_sorted
    rm $tmpR2_sorted
    rm $tmpbc_sorted
}

singleBC_trust4(){
    ## function cannot access global variable when parsing to GNU parallel
    local bc=$1
    local trust4_wd=$2


    local umiCounts=$(grep -c ">" ${trust4_wd}/umi/${bc}.sc_umi.fa)
    if [[ $umiCounts -eq 0 ]]
    then
        local barcode_umi_opt="--barcode ${trust4_wd}/barcode/${bc}.sc_bc.fa"
    else
        local barcode_umi_opt="--barcode ${trust4_wd}/barcode/${bc}.sc_bc.fa --UMI ${trust4_wd}/umi/${bc}.sc_umi.fa"
    fi

    if [[ $R2_Only == true ]]
    then
        local trust4_input_opt="-u ${trust4_wd}/readsFastq/${bc}.sc_reads.R2.fq"
    else
        local trust4_input_opt="-1 ${trust4_wd}/readsFastq/${bc}.sc_reads.R1.fq -2 ${trust4_wd}/readsFastq/${bc}.sc_reads.R2.fq"
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
    else
        echo "No reads extracted for $bc"
    fi
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
    parallel -j $threads singleBC_trust4 {} $workDir

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
rm -rf $workDir
