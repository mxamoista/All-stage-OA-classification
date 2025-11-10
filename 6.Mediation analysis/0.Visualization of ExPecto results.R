library(dplyr)
library(biomaRt)
library(data.table)
library(qqman)
library(tidyverse)
library(ggplot2)
library(ggrepel)

#Import ExPecto predict result

# TC1.effect <- read.csv('./data/Expecto_predict/TC1_predict.csv') 
# setDT(TC1.effect)[, index := index + 1]
# TC1.effect <- TC1.effect[, -c(1)]

TC5.effect <- read.csv('./data/Expecto_predict/ebi-a-GCST007090/TC5_predict.csv')
setDT(TC5.effect)[, index := index + 1]
TC5.effect <- TC5.effect[, -c(1)]

#Rename the colnames
# colnames(TC1.effect)[2:11] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE")
colnames(TC5.effect)[2:11] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "VALUE")

#Connect to human genes dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Get the mapping 
# results_TC1 <- getBM(
#   attributes = c("ensembl_gene_id", "external_gene_name"),
#   filters = "ensembl_gene_id",  
#   values = TC1.effect$gene,     
#   mart = ensembl
# )

results_TC5 <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",  
  values = TC5.effect$gene,     
  mart = ensembl
)

#Merge with results
# TC1.effect$external_gene_name <- results_TC1$external_gene_name[match(TC1.effect$gene, results_TC1$ensembl_gene_id)]
TC5.effect$external_gene_name <- results_TC5$external_gene_name[match(TC5.effect$gene, results_TC5$ensembl_gene_id)]

#Draw the manhattan plot
# draw_manhattan(TC1.effect, "rs376359757")
draw_manhattan(TC5.effect, "rs376359757")


