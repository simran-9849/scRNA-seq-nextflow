#! /usr/bin/env Rscript

## This script is to merge saturation data and preseqR input data
## into one JSON file

## usage: combine_saturation_data.R saturation.json UMI_hist.tsv gene_hist.tsv totalGeneCountFile outputJSON

library(jsonlite)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

pilot_json <- fromJSON(args[1])

reads_hist <- read_tsv(
    args[2],
    col_names = c("freq", "counts")
) %>%
    as.matrix

genes_hist <- read_tsv(
    args[3],
    col_names = c("freq", "counts")
) %>%
    as.matrix

totalGeneCount <- read_tsv(args[4], col_names = "count") %>% pull(count)

d <- list(pilot_metrics = pilot_json,
          prediction_raw_readsData = reads_hist,
          prediction_raw_genesData = genes_hist,
          totalGeneCount = totalGeneCount)

write_json(d, path = args[5], pretty = TRUE)
