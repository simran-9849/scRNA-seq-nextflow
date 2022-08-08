#! /usr/bin/env Rscript

## update summary.csv for starsolo output of multiple gene reads

## usage:
## udpate_starsolo_summary.R CellReads.stats mult_barcodes.tsv mult_10X_dir mult_saturation.json Summary.unique.csv Summary.multiple.csv

library(tidyverse)
library(jsonlite)
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)

d <- read_tsv(args[1])

cells <- read_tsv(args[2], col_names = "barcode")

reads_in_gene_total <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    mutate(n=featureU+featureM) %>%
    pull(n) %>% sum

reads_in_gene_in_cells <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=featureU+featureM) %>%
    pull(n) %>% sum

fraction_of_reads_in_cells <- reads_in_gene_in_cells/reads_in_gene_total

mean_reads_per_cell <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=featureU+featureM) %>%
    pull(n) %>% mean %>% round

median_reads_per_cell <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=featureU+featureM) %>%
    pull(n) %>% median %>% round

UMI_in_cells <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=nUMIunique+nUMImulti) %>%
    pull(n) %>% sum

mean_UMI_per_cell <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=nUMIunique+nUMImulti) %>%
    pull(n) %>% mean %>% round

median_UMI_per_cell <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=nUMIunique+nUMImulti) %>%
    pull(n) %>% median %>% round

mean_gene_per_cell <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=nGenesUnique+nGenesMulti) %>%
    pull(n) %>% mean %>% round

median_gene_per_cell <- d %>%
    filter(CB!="CBnotInPasslist") %>%
    filter(CB %in% cells$barcode) %>%
    mutate(n=nGenesUnique+nGenesMulti) %>%
    pull(n) %>% median %>% round

saturation <- fromJSON(args[4]) %>%
    as_tibble %>%
    filter(percentage == 1) %>%
    pull(saturation) %>%
    as.numeric()

summary_unique <- read_csv(args[5], col_names = c("term", "value"))

## create new summary file
totalInputReads <- summary_unique %>%
    filter(term == "Number of Reads") %>%
    pull(value)

reads_with_valid_barcodes <- summary_unique %>%
    filter(term == "Reads With Valid Barcodes") %>%
    pull(value)

q30_in_R1 <- summary_unique %>%
    filter(term == "Q30 Bases in CB+UMI") %>%
    pull(value)

q30_in_R2 <- summary_unique %>%
    filter(term == "Q30 Bases in RNA read") %>%
    pull(value)

genome_u_m <- summary_unique %>%
    filter(term == "Reads Mapped to Genome: Unique+Multiple") %>%
    pull(value)

genome_u <- summary_unique %>%
    filter(term == "Reads Mapped to Genome: Unique") %>%
    pull(value)

feature_u_m <- summary_unique %>%
    filter(str_detect(term , regex("Reads Mapped to Gene\\w*: Unique\\+Multiple Gene\\w*"))) %>%
    pull(value)

feature_u <- summary_unique %>%
    filter(str_detect(term, regex("Reads Mapped to Gene\\w*: Unique Gene\\w*"))) %>%
    pull(value)

## calc total gene detected
mult_10X_dir <- args[3]
m <- Read10X(mult_10X_dir)
total_gene_detected <- sum(apply(m, 1, sum) > 0)

summary_multi <- tribble(
    ~term, ~value,
    "Number of Reads", totalInputReads,
    "Reads With Valid Barcodes", reads_with_valid_barcodes,
    "Sequencing Saturation", saturation,
    "Q30 Bases in CB+UMI", q30_in_R1,
    "Q30 Bases in RNA read", q30_in_R2,
    "Reads Mapped to Genome: Unique+Multiple", genome_u_m,
    "Reads Mapped to Genome: Unique", genome_u,
    "Reads Mapped to Gene: Unique+Multiple Gene", feature_u_m,
    "Reads Mapped to Gene: Unique Gene", feature_u,
    "Estimated Number of Cells", nrow(cells),
    "Reads in Cells Mapped to Gene", reads_in_gene_in_cells,
    "Fraction of Reads in Cells", reads_in_gene_in_cells/reads_in_gene_total,
    "Mean Reads per Cell", mean_reads_per_cell,
    "Median Reads per Cell", median_reads_per_cell,
    "UMIs in Cells", UMI_in_cells,
    "Mean UMI per Cell", mean_UMI_per_cell,
    "Median UMI per Cell", median_UMI_per_cell,
    "Mean Gene per Cell", mean_gene_per_cell,
    "Median Gene per Cell", median_gene_per_cell,
    "Total Gene Detected", total_gene_detected,
)

write_csv(summary_multi, args[6], col_names = FALSE)
