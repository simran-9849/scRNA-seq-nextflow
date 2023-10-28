#! /usr/bin/env bash


usage () {
    cat<<-EOF
	$(basename $0) <trust4_report> <trust4_airr> <trust4_toassemble_bc> <cellType> <starsoloSummary> <metricsOutput> <kneeDataOutput> <cloneTypeResultFile>
	
	This script is to perform some statistics from TRUST4 result
	and generate a metrics json output. User will have to provide
	trust4_report, trust4_toassemble_bc.fa, kneeData.tsv, cellType
	and starsolo summary.csv as input, and define the names of
	three output files.
	
	cellType only support "VDJ-T" and "VDJ-B"
	
	example: $(basename $0) trust4_report.tsv \\
	                        trust4_airr.tsv \\
	                        trust4_toassemble_bc.fa \\
	                        sampleID_VDJ-T.kneeData.tsv \\
	                        VDJ-T \\
	                        sampleID_VDJ-T.Summary.unique.csv \\
	                        sampleID_VDJ-T.cellOut.tsv \\
	                        sampleID_VDJ-T.vdj_metrics.json \\
	                        sampleID_VDJ-T.cloneType_out.tsv
	EOF
}

if [[ $# -eq 0 ]] || [[ $1 == "-h" ]]
then
    usage
    exit 0
fi

trust4_report=$1
trust4_airr=$2
trust4_toassemble_bc=$3
kneeInput=$4
cellType=$5
starsolo_summary=$6
cellOut=$7
metricsOut=$8
cloneTypeResult=$9

sampleID=${trust4_report%%_VDJ*}

if [[ $cellType == "VDJ-T" ]]
then
    cellName="abT"
elif [[ $cellType == "VDJ-B" ]]
then
    cellName="B"
else
    echo "cellType only supports VDJ-T and VDJ-B" &>2
    exit 1
fi


## cells were defined as chain1 (IGH or TRB) and chain2 (IGK/L or TRA) has 3 UMI in total
cellBC=$(mktemp -p ./)
awk -v cellName=$cellName '
    BEGIN{
        print "CB\tChain1_UMI\tChain2_UMI"
    }
    {
        split($3, chain1, ",");
        split($4, chain2, ",");
        chain1_umi=chain1[7];
        chain2_umi=chain2[7];
        if($1!="-" && $2==cellName && chain1_umi+chain2_umi>=3){
            print $1"\t"chain1_umi"\t"chain2_umi
        }
    }
    ' $trust4_report > $cellOut
awk 'NR>1{print $1}' $cellOut > $cellBC
cellNum=$(wc -l $cellBC | awk '{print $1}')

## calculate UMI in cells and background barcodes
## include more information in the kneeData
## columns:
## 1.barcode 2.chain1_umi 3.chain1_genotype 4.chain1_cdr3 5.chain2_umi 6.chain2_genotype 7.chain2_cdr3
##awk 'NR>1{if($3!="*"){split($3, chain1, ",");  chain1_umi=chain1[7]; chain1_genotype=chain[1]"+"chain1[2]"+"chain1[3]"+"chain1[4]; chain1_cdr3=chain1[6];}else{chain1_umi=0; chain1_genotype="*"; chain1_cdr3="*";} if($4!="*"){split($4, chain2, ","); chain2_umi=chain2[7]; chain2_genotype=chain2[1]"+"chain2[2]"+"chain2[3]"+"chain2[4]; chain2_cdr3=chain2[6];}else{chain2_umi=0; chain2_genotype="*"; chain2_cdr3="*"};print $1"\t"chain1_umi"\t"chain1_genotype"\t"chain1_cdr3"\t"chain2_umi"\t"chain2_genotype"\t"chain2_cdr3}' $trust4_report |
##    awk 'BEGIN{print "CB\tchain1_UMI\tchain1_genotype\tchain1_cdr3\tchain2_UMI\tchain2_genotype\tchain2_cdr3"}{print}' > $kneeInput

## Extract read and umi in cells
readCBList=$(mktemp -p ./)
awk '
    $1~/^>/{
        readID=substr($1,2);
        getline;
        CB=$1;
        print readID"\t"CB
    }
    ' $trust4_toassemble_bc > $readCBList
readCellList=$(mktemp -p ./)
awk 'ARGIND==1{cell[$1]}ARGIND==2{if($2 in cell){print}}' $cellBC $readCBList > $readCellList

## calculate mean/median reads per cell
meanReadsPerCell=$(awk '{print $2}' $readCellList | sort | uniq -c | awk '{print $1}' | sort -k 1,1rn | awk 'BEGIN{t=0}{t+=$1}END{print t/NR}')
medianReadsPerCell=$(awk '{print $2}' $readCellList | sort | uniq -c | awk '{a[$1]}END{asorti(a); n=length(a); if(n%2==1){print a[n/2+0.5]}else{print (a[n/2]+a[n/2+1])/2}}')

## total reads of barcode
totalReadsInCB=$(wc -l $readCBList| awk '{print $1}')
totalReadsInCell=$(wc -l $readCellList| awk '{print $1}')
fractionReadsInCells=$(awk -v cell=$totalReadsInCell -v total=$totalReadsInCB 'BEGIN{print cell/total}')

## calculate mean/median UMIs per cell, UMI of chain1 (heavy) and chain2 (light)
meanUMIsPerCell=$(awk 'ARGIND==1{bc[$1]}ARGIND==2{if($1 in bc){print}}' $cellBC $kneeInput| sort -k 2,2rn | awk 'BEGIN{t=0}{t+=$2}END{print t/NR}')
medianUMIsPerCell=$(awk 'ARGIND==1{bc[$1]}ARGIND==2{if($1 in bc){print}}' $cellBC $kneeInput | sort -k 2,2rn | awk '{a[$2]}END{asorti(a); n=length(a); if(n%2==1){print a[n/2+0.5]}else{print (a[n/2]+a[n/2+1])/2}}')
medianUMIsChain1=$(awk 'ARGIND==1{bc[$1]}ARGIND==2{split($3, chain1, ","); split($4, chain2, ","); chain1_umi=chain1[7]; chain2_umi=chain2[7]; if($1 in bc){print $1"\t"chain1_umi}}' $cellBC $trust4_report | awk '{a[$2]}END{asorti(a); n=length(a); if(n%2==1){print a[n/2+0.5]}else{print (a[n/2]+a[n/2+1])/2}}')
medianUMIsChain2=$(awk 'ARGIND==1{bc[$1]}ARGIND==2{split($3, chain1, ","); split($4, chain2, ","); chain1_umi=chain1[7]; chain2_umi=chain2[7]; if($1 in bc){print $1"\t"chain2_umi}}' $cellBC $trust4_report | awk '{a[$2]}END{asorti(a); n=length(a); if(n%2==1){print a[n/2+0.5]}else{print (a[n/2]+a[n/2+1])/2}}')

## Generate clonetype result table
awk '
    ARGIND==1{
        bc[$1]
    }
    ARGIND==2{
        if($1 in bc){print}
    }
    ' $cellBC $trust4_report |
    awk '
    {
        if($3!="*"){
            split($3, chain1, ",");
            s1=chain1[6]
            gsub(/\*.*/, "", chain1[1])
            gsub(/\*.*/, "", chain1[2])
            gsub(/\*.*/, "", chain1[3])
        }else{
            s1="NA"
        };
        if($4!="*"){
            split($4, chain2, ",");
            s2=chain2[6]
            gsub(/\*.*/, "", chain2[1])
            gsub(/\*.*/, "", chain2[3])
        }else{
            s2="NA"
        }
        name1=substr($3,1,3);
        name2=substr($4,1,3);
        print name1":"s1";"name2":"s2"\t"chain1[1]"|"chain1[2]"|"chain1[3]"\t"chain2[1]"|"chain2[3]
    }
    ' | sort | uniq -c | sort -k 1,1rn |
    awk -v cellNum=$cellNum 'BEGIN{print "cloneType\tHchain_VDJ\tLchain_VDJ\tCellCount\tFrequency"}{print $2"\t"$3"\t"$4"\t"$1"\t"$1/cellNum}' > $cloneTypeResult

## Extract metrics from starsolo summary output
extractSummaryTerm(){
    local summaryFile=$1
    local term=$2
    local tmpCount=$(grep "$term" $summaryFile | wc -l)
    local outValue="NoInfor"
    if [[ $tmpCount -gt 0 ]]
    then
        outValue=$(awk -F"," -v term="$term" '$1==term{print $2}' $summaryFile)
    fi
    echo $outValue
}

## Number of Reads
totalRawReads=$(extractSummaryTerm $starsolo_summary "Number of Reads")
## Reads With Valid Barcodes
validBCreads=$(extractSummaryTerm $starsolo_summary "Reads With Valid Barcodes")
## Sequencing Saturation
saturation=$(extractSummaryTerm $starsolo_summary "Sequencing Saturation")
## Q30 Bases in CB+UMI
q30InCBandUMI=$(extractSummaryTerm $starsolo_summary "Q30 Bases in CB+UMI")
## Q30 Bases in RNA read
q30InRNA=$(extractSummaryTerm $starsolo_summary "Q30 Bases in RNA read")
## Reads mapped to genome (U+M)
totalMappedReads=$(extractSummaryTerm $starsolo_summary "Reads Mapped to Genome: Unique+Multiple")
## Reads mapped to genome (U)
uniquelyMappedReads=$(extractSummaryTerm $starsolo_summary "Reads Mapped to Genome: Unique")

## Total cloneTypes detected
totalCloneTypes=$(wc -l $cloneTypeResult| awk '{print $1-1}')

pairingCellNum=$(awk 'BEGIN{t=0}NR>1 && $1!~/*/{t+=$4}END{print t}' $cloneTypeResult)
pairingRate=$(awk -v pairingCellNum=$pairingCellNum -v cellNum=$cellNum 'BEGIN{print pairingCellNum/cellNum}')

countProductiveCells(){
    local cellBC=$1
    local airr=$2
    local report=$3
    local type=$4
    awk -v type=$type '
    ARGIND==1{
        cells[$1]
    }
    ARGIND==2&&FNR>1{
        productive[$1]=$4
    }
    ARGIND==3{
        if($1 in cells){
            split($3, chain1, ",");
            split($4, chain2, ",");
            name1=substr(chain1[1], 1, 3);
            name2=substr(chain2[1], 1, 3);
            chain1_id=chain1[8];
            chain2_id=chain2[8];
            if(type == "IGH" || type == "TRB" || type == "TRD"){
                if(name1==type && productive[chain1_id]=="T"){print $1}
            }else if(type == "IGK" || type == "IGL" || type == "TRA" || type == "TRG"){
                if(name2==type && productive[chain2_id]=="T"){print $1}
            }else if(type == "IGH_IGK" || type == "IGH_IGL" || type == "TRB_TRA" || type == "TRG_TRD"){
                split(type, typeArray, "_");
                if(name1==typeArray[1] && name2==typeArray[2] && productive[chain1_id]=="T" && productive[chain2_id]=="T"){print $1}
            }
        }
    }
    ' $cellBC $airr $report | wc -l
}

if [[ $cellType == "VDJ-B" ]]
then
    ## Calculate cells with productive V-J/V-D-J fragments (full length, with in-frame CDR3 amino acids)
    ## cells with productive IGH chain, no matter if complete VDJ was found, no need to have both chains
    cellsWithProductiveChainIGH=$(countProductiveCells $cellBC $trust4_airr $trust4_report "IGH")
    ## cells with productive IGK chain, no matter if complete pair was found
    cellsWithProductiveChainIGK=$(countProductiveCells $cellBC $trust4_airr $trust4_report "IGK")
    ## cells with productive IGL chain, no matter if complete pair was found
    cellsWithProductiveChainIGL=$(countProductiveCells $cellBC $trust4_airr $trust4_report "IGL")
    ## cells with both productive productive IGH and IGK
    cellsWithProductiveChainIGKIGH=$(countProductiveCells $cellBC $trust4_airr $trust4_report "IGH_IGK")
    ## cells with both productive productive IGH and IGL
    cellsWithProductiveChainIGLIGH=$(countProductiveCells $cellBC $trust4_airr $trust4_report "IGH_IGL")


    ## Generate json metrics
    jq -n \
       --arg sampleName "$sampleID" \
       --arg totalRawReads "$totalRawReads" \
       --arg validBCreads "$validBCreads" \
       --arg saturation "$saturation" \
       --arg q30InCBandUMI "$q30InCBandUMI" \
       --arg q30InRNA "$q30InRNA" \
       --arg totalMappedReads "$totalMappedReads" \
       --arg uniquelyMappedReads "$uniquelyMappedReads" \
       --arg cells "$cellNum" \
       --arg meanReadsPerCell "$meanReadsPerCell" \
       --arg medianReadsPerCell "$medianReadsPerCell" \
       --arg totalReadsInCell "$totalReadsInCell" \
       --arg fractionReadsInCells "$fractionReadsInCells" \
       --arg meanUMIsPerCell "$meanUMIsPerCell" \
       --arg medianUMIsPerCell "$medianUMIsPerCell" \
       --arg medianUMIsChain1 "$medianUMIsChain1" \
       --arg medianUMIsChain2 "$medianUMIsChain2" \
       --arg totalCloneTypes "$totalCloneTypes" \
       --arg pairingRate "$pairingRate" \
       --arg cellsWithProductiveChainIGH "$cellsWithProductiveChainIGH" \
       --arg cellsWithProductiveChainIGK "$cellsWithProductiveChainIGK" \
       --arg cellsWithProductiveChainIGL "$cellsWithProductiveChainIGL" \
       --arg cellsWithProductiveChainIGKIGH "$cellsWithProductiveChainIGKIGH" \
       --arg cellsWithProductiveChainIGLIGH "$cellsWithProductiveChainIGLIGH" \
       '{sampleName: $sampleName, totalRawReads: $totalRawReads, validBCreads: $validBCreads, saturation: $saturation, q30InCBandUMI: $q30InCBandUMI, q30InRNA: $q30InRNA, totalMappedReads: $totalMappedReads, uniquelyMappedReads: $uniquelyMappedReads, cells: $cells, totalReadsInCell: $totalReadsInCell, meanReadsPerCell: $meanReadsPerCell, medianReadsPerCell: $medianReadsPerCell, fractionReadsInCells: $fractionReadsInCells, meanUMIsPerCell: $meanUMIsPerCell, medianUMIsPerCell: $medianUMIsPerCell, medianUMIsChain1: $medianUMIsChain1, medianUMIsChain2: $medianUMIsChain2, totalCloneTypes: $totalCloneTypes, pairingRate: $pairingRate, cellsWithProductiveChainIGH: $cellsWithProductiveChainIGH, cellsWithProductiveChainIGK: $cellsWithProductiveChainIGK, cellsWithProductiveChainIGL: $cellsWithProductiveChainIGL, cellsWithProductiveChainIGKIGH: $cellsWithProductiveChainIGKIGH, cellsWithProductiveChainIGLIGH: $cellsWithProductiveChainIGLIGH}' > $metricsOut

elif [[ $cellType == "VDJ-T" ]]
then
    ## Calculate cells with productive V-J/V-D-J fragments (full length, with in-frame CDR3 amino acids)
    ## cells with productive IGH chain, no matter if complete pair was found
    cellsWithProductiveChainTRB=$(countProductiveCells $cellBC $trust4_airr $trust4_report "TRB")
    cellsWithProductiveChainTRD=$(countProductiveCells $cellBC $trust4_airr $trust4_report "TRD")
    ## cells with productive IGK chain, no matter if complete pair was found
    cellsWithProductiveChainTRA=$(countProductiveCells $cellBC $trust4_airr $trust4_report "TRA")
    ## cells with productive IGL chain, no matter if complete pair was found
    cellsWithProductiveChainTRG=$(countProductiveCells $cellBC $trust4_airr $trust4_report "TRG")
    ## cells with both productive TRA and TRB
    cellsWithProductiveChainTRATRB=$(countProductiveCells $cellBC $trust4_airr $trust4_report "TRB_TRA")
    ## cells with both productive TRD and TRG
    cellsWithProductiveChainTRGTRD=$(countProductiveCells $cellBC $trust4_airr $trust4_report "TRD_TRG")

    jq -n \
       --arg sampleName "$sampleID" \
       --arg totalRawReads "$totalRawReads" \
       --arg validBCreads "$validBCreads" \
       --arg saturation "$saturation" \
       --arg q30InCBandUMI "$q30InCBandUMI" \
       --arg q30InRNA "$q30InRNA" \
       --arg totalMappedReads "$totalMappedReads" \
       --arg uniquelyMappedReads "$uniquelyMappedReads" \
       --arg cells "$cellNum" \
       --arg meanReadsPerCell "$meanReadsPerCell" \
       --arg medianReadsPerCell "$medianReadsPerCell" \
       --arg totalReadsInCell "$totalReadsInCell" \
       --arg fractionReadsInCells "$fractionReadsInCells" \
       --arg meanUMIsPerCell "$meanUMIsPerCell" \
       --arg medianUMIsPerCell "$medianUMIsPerCell" \
       --arg medianUMIsChain1 "$medianUMIsChain1" \
       --arg medianUMIsChain2 "$medianUMIsChain2" \
       --arg totalCloneTypes "$totalCloneTypes" \
       --arg pairingRate "$pairingRate" \
       --arg cellsWithProductiveChainTRB "$cellsWithProductiveChainTRB" \
       --arg cellsWithProductiveChainTRD "$cellsWithProductiveChainTRD" \
       --arg cellsWithProductiveChainTRA "$cellsWithProductiveChainTRA" \
       --arg cellsWithProductiveChainTRG "$cellsWithProductiveChainTRG" \
       --arg cellsWithProductiveChainTRATRB "$cellsWithProductiveChainTRATRB" \
       --arg cellsWithProductiveChainTRGTRD "$cellsWithProductiveChainTRGTRD" \
       '{sampleName: $sampleName, totalRawReads: $totalRawReads, validBCreads: $validBCreads, saturation: $saturation, q30InCBandUMI: $q30InCBandUMI, q30InRNA: $q30InRNA, totalMappedReads: $totalMappedReads, uniquelyMappedReads: $uniquelyMappedReads, cells: $cells, totalReadsInCell: $totalReadsInCell, meanReadsPerCell: $meanReadsPerCell, medianReadsPerCell: $medianReadsPerCell, fractionReadsInCells: $fractionReadsInCells, meanUMIsPerCell: $meanUMIsPerCell, medianUMIsPerCell: $medianUMIsPerCell, medianUMIsChain1: $medianUMIsChain1, medianUMIsChain2: $medianUMIsChain2, totalCloneTypes: $totalCloneTypes, pairingRate: $pairingRate, cellsWithProductiveChainTRB: $cellsWithProductiveChainTRB, cellsWithProductiveChainTRD: $cellsWithProductiveChainTRD, cellsWithProductiveChainTRA: $cellsWithProductiveChainTRA, cellsWithProductiveChainTRG: $cellsWithProductiveChainTRG, cellsWithProductiveChainTRATRB: $cellsWithProductiveChainTRATRB, cellsWithProductiveChainTRGTRD: $cellsWithProductiveChainTRGTRD}' > $metricsOut

else
    echo "cellType only supports VDJ-T and VDJ-B" &>2
    exit 1
fi

## remove temp files
rm $cellBC $readCBList $readCellList
