#!/usr/bin/env Rscript
#  
#  
#  Author: Yuan Zhen, Mail:yuanzhen@thunder-bio.com
#  
# this scripts can filter your VDJ results from trust4
# need two file *barcode_airr.tsv and *barcode_report.tsv for input
# and a file as filtered results for output
# 
# the filtered parameters as flow:
# 1.Coding region has an open reading frame
# 2.No defect in the start codon, splicing sites or regulatory elements
# 3.No internal stop codons
# 4.An in-frame junction region
# 5.NOT USE,Sequence_alignment includes both the first V gene codon that encodes the mature polypeptide chain  and the last complete codon of the J gene,not use because of 3' scRNA-se
# 6.Number of reads contributing to the UMI consensus or contig assembly for this sequence need greater than 3
# 7.The sum of umi of chain one and chain two is greater than 6
# 8.No stop coden or N base in CDR3 base sequence
# 9.If barcode X's two chains CDR3s are identical to another barcode Y, and X's chain abundance is significantly lower(100 times) than Y's, filter X
#
#
#


# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
  make_option(c("-r", "--barcode_report"), type = "character",
              help="*barcode_report.tsv from TRUST4"),
  make_option(c("-a", "--barcode_airr"),type = "character", 
              help = "*barcode_airr.tsv from TRUST4"),
  make_option(c("-o", "--output"),type = "character",
              help = "output filename")
)




# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$barcode_airr) | is.null(opt$barcode_report)){
  stop("VDJ_trust4_filter wrong Use -h or --help for help")
}
if(!file.exists(opt$barcode_airr)){
  stop(sprintf("your *barcode_airr.tsv: ( %s ) does not exist,check your path and files", opt$barcode_airr))
}else{
  barcode_airr=opt$barcode_airr
}

if(!file.exists(opt$barcode_report)){
  stop("your *barcode_report.tsv: ( %s ) not exits ,check your path and files",opt$barcode_report)
}else{
  barcode_report=opt$barcode_report
}



suppressPackageStartupMessages(library('tidyverse'))

#suppressPackageStartupMessages(library('magrittr'))
suppressWarnings(airr <-read_delim(file=barcode_airr, 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE,show_col_types = FALSE))
barcode_star <- read_table(barcode_report,show_col_types = FALSE) %>% 
  filter(!`#barcode`=="-")
message("We have detect raw barcode ",nrow(barcode_star))

airr$sequence_id=str_remove_all(airr$sequence_id,"_[0-9]+")
airr=airr[!str_detect(airr$sequence_id,"-"),]



# 1.Coding region has an open reading frame
# 2.No defect in the start codon, splicing sites or regulatory elements
# 3.No internal stop codons
# 4.An in-frame junction region
airr=airr[airr$productive,]
message("We will filter productive contig:")
message("1.Coding region has an open reading frame")
message("2.No defect in the start codon, splicing sites or regulatory elements")
message("3.No internal stop codons")
message("4.An in-frame junction region")

# 5.NOT USE,Sequence_alignment includes both the first V gene codon that encodes the mature polypeptide chain  and the last complete codon of the J gene
#airr=airr[airr$complete_vdj,]

#6.Number of reads contributing to the UMI consensus or contig assembly for this sequence need greater than 3
airr=airr[airr$consensus_count>3,]
message("6.Number of reads contributing to the UMI consensus or contig assembly for this sequence need greater than 3")





# 7.The sum of umi of chain one and chain two is greater than 6
barcode_star=filter(barcode_star,`#barcode` %in% airr$sequence_id)
message("After step 1,2,3,4,6 detect ",nrow(barcode_star)," barcode")
barcode_star$umi1=barcode_star$chain1 %>% ifelse(. == "*","1,1,1,1,1,NA,0",.)  %>% str_split(",") %>% sapply("[",7) %>% as.numeric() 
barcode_star$umi2=barcode_star$chain2 %>% ifelse(. == "*","1,1,1,1,1,NA,0",.)  %>% str_split(",") %>% sapply("[",7) %>% as.numeric()
barcode_star$umi=barcode_star$umi1+barcode_star$umi2
barcode_star= filter(barcode_star,umi>6)
message("7.The sum of umi of chain one and chain two is greater than 6")
message("After filter step 7 the barcode number is: ",nrow(barcode_star))


# 9.If barcode X's two chains CDR3s are identical to another barcode Y, and X's chain abundance is significantly lower(100 times) than Y's, filter X
barcode_star$aa_cdr3_1=barcode_star$chain1 %>% ifelse(. == "*","1,1,1,1,NA,CGGGGGGGGG,0",.)  %>% str_split(",") %>% sapply("[",6)
barcode_star$aa_cdr3_2=barcode_star$chain2 %>% ifelse(. == "*","1,1,1,1,NA,CGGGGGGGGG,0",.)  %>% str_split(",") %>% sapply("[",6)
barcode_star <-barcode_star %>%
  group_by(aa_cdr3_1, aa_cdr3_2) %>%
  filter(!(umi / max(umi) <0.01)) %>% ungroup()
message("9.If barcode X's two chains CDR3s are identical to another barcode Y, and X's chain abundance is significantly lower(100 times) than Y's, filter X")
message("After filter step 9 the barcode number is: ",nrow(barcode_star))

# 8.No stop coden or N base in CDR3 base sequence
barcode_star$dna_cdr3_1=barcode_star$chain1 %>% ifelse(. == "*","1,1,1,1,GGGGGGGGGGGGGGGGGGGGGGGGGGGG,NA,0",.)  %>% str_split(",") %>% sapply("[",5) 
barcode_star$dna_cdr3_2=barcode_star$chain2 %>% ifelse(. == "*","1,1,1,1,1,GGGGGGGGGGGGGGGGGGGGGGGGGGGG,0",.)  %>% str_split(",") %>% sapply("[",5)
barcode_star=barcode_star[str_detect(barcode_star$dna_cdr3_1,"N",negate = T),]
barcode_star=barcode_star[str_detect(barcode_star$dna_cdr3_2,"N",negate = T),]


barcode_star=barcode_star[str_detect(barcode_star$dna_cdr3_1,"UAG | UAA | UGA",negate = T) ,]
barcode_star=barcode_star[str_detect(barcode_star$dna_cdr3_1,"UAG | UAA | UGA",negate = T) ,]
barcode_star$length_1=barcode_star$chain1 %>% ifelse(. == "*","1,1,1,1,GGGGGGGGGGGGGGGGGGGGGGGGGGGG,NA,0",.)  %>% str_split(",") %>% sapply("[",5) %>% str_length() 
barcode_star$length_2=barcode_star$chain2 %>% ifelse(. == "*","1,1,1,1,1,GGGGGGGGGGGGGGGGGGGGGGGGGGGG,0",.)  %>% str_split(",") %>% sapply("[",5) %>% str_length() 
message("8.No stop coden or N base in CDR3 base sequence")
message("After filter step 8 the barcode number is: ",nrow(barcode_star))


write_tsv(barcode_star, opt$output)

