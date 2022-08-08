#! /usr/bin/env Rscript

## update UMIperCellSorted.txt for starsolo output of multiple gene reads

## usage:
## udpate_umi_file.R CellReads.stats UMIperCellSorted.multiple.txt

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

d <- read_tsv(args[1])

d %>%
    filter(CB!="CBnotInPasslist") %>%
    mutate(n=nUMIunique+nUMImulti) %>%
    filter(n > 0) %>%
    select(n) %>%
    arrange(desc(n)) %>%
    write_tsv(args[2], col_names = FALSE)
