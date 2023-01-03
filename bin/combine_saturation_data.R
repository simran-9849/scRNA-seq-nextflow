#! /usr/bin/env Rscript

## This script is to merge saturation data and preseqR input data
## into one JSON file

## usage: combine_saturation_data.R sampleID saturation.json UMI_hist.tsv gene_hist.tsv totalGeneCountFile outputJSON

library(jsonlite)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

sampleID <- args[1]

saturation_json <- fromJSON(args[2])

reads_hist <- read_tsv(
    args[3],
    col_names = c("freq", "counts")
) %>%
    as.matrix

genes_hist <- read_tsv(
    args[4],
    col_names = c("freq", "counts")
) %>%
    as.matrix

totalGeneCount <- read_tsv(args[5], col_names = "count") %>% pull(count)

d <- list(
    sampleID = sampleID,
    saturation_data = saturation_json,
    preseqR_raw_UMI_data = reads_hist,
    preseqR_raw_gene_data = genes_hist,
    totalGeneCount = totalGeneCount
)

write_json(d, path = args[6], pretty = TRUE)
