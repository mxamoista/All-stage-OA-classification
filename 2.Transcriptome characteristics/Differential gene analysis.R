library(ggplot2)
library(dplyr)
library(DESeq2)
library(biomaRt)
library(edgeR)
library(readr)
library(tidyr)
# Connect to the human genes dataset
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#==========================
#load datasets
BRBseq1 <- read.table("./Batch/BRBseq1/expressions.tsv.gz",header = TRUE , row.names =1)
BRBseq2 <- read.table("./Batch/BRBseq2/expressions.tsv.gz",header = TRUE , row.names =1)
RNAseq1 <- read.table("./Batch/RNAseq1/expressions.tsv.gz",header = TRUE , row.names =1)
RNAseq2 <- read.table("./Batch/RNAseq2/expressions.tsv.gz",header = TRUE , row.names =1)

matrix <- cbind(RNAseq1,RNAseq2,BRBseq1,BRBseq2)
QC <- read.table("./Batch/Quality_check_results.tsv",header = TRUE )
identical(QC$Sample_ID,colnames(matrix)) #[1] TRUE
matrix <- matrix[,!(QC$Final_Status) %in% c("Failed") ] 
QC <- QC[!(QC$Final_Status) %in% c("Failed") ,]

#load subtype
annotation_col <- read.table("./annotation_col/annotation_col.txt",header = TRUE , row.names =1)

#==========================
# Batch detection by umap
library(umap)
library(ggplot2)
umap_df <- umap(t(matrix))
umap_df = cbind.data.frame(umap_df$layout,QC$Batch)
colnames(umap_df) = c("UMAP1","UMAP2","Batch")
umap_df$Batch = as.factor(umap_df$Batch)
umap_df$Batch <- factor(umap_df$Batch, levels = c("BRBseq1","BRBseq2","RNAseq1","RNAseq2"))
ggplot(umap_df,aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Batch),size=4) + theme_classic(base_size = 20)+scale_color_manual(values=c("BRBseq1"="#00CED1","BRBseq2"="#0000CD","RNAseq1"="#984EA3","RNAseq2"="#FF0000")) + guides(colour = guide_legend(override.aes = list(size=4)))

#==========================
##Remove batch effect
library(sva)
batch <- c(rep('RNAseq1',87),rep('RNAseq2',14),rep('BRBseq1',20),rep('BRBseq2',18))

# Using ComBat to remove batch effect
combat_edata = ComBat(dat = matrix, batch=batch, mean.only = F) #mean.only=F work better than mean.only=T

# Using the UMAP plot to verify that the batch effect has been removed
umap_df2 <- umap(t(combat_edata))
umap_df2 = cbind.data.frame(umap_df2$layout,QC$Batch,annotation_col$Subtype)
colnames(umap_df2) = c("UMAP1","UMAP2","Batch","Subtype")
umap_df2$Batch = as.factor(umap_df2$Batch)
umap_df2$Batch <- factor(umap_df2$Batch, levels = c("BRBseq1","BRBseq2","RNAseq1","RNAseq2"))
umap_df2$Subtype = as.factor(umap_df2$Subtype)
umap_df2$Subtype <- factor(umap_df2$Subtype, levels = c("TC1","TC2","TC3","TC4","TC5","TC6"))
ggplot(umap_df2,aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Batch),size=4) + theme_classic(base_size = 20)+scale_color_manual(values=c("BRBseq1"="#00CED1","BRBseq2"="#0000CD","RNAseq1"="#984EA3","RNAseq2"="#FF0000")) + guides(colour = guide_legend(override.aes = list(size=4)))
ggplot(umap_df2,aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Subtype),size=4) + theme_classic(base_size = 20)+scale_color_manual(values=c("TC1"="#87CEEB","TC2"="#984EA3","TC3"="#FF0000","TC4"="#0000CD","TC5"="#EE82EE","TC6"="#00CED1")) + guides(colour = guide_legend(override.aes = list(size=4)))

# Using varFilter to remove low-variance genes 
library(genefilter)
var_cutoff <- 0.25
combat_edata.filter <- varFilter(as.matrix(combat_edata), var.cutoff = 0.25,filterByQuantile=TRUE)

# Convert negative values to 0
combat_edata.filter[combat_edata.filter < 0] <- 0 # fix negative values

# Select batch-corrected counts of TC1-6 samples
identical(rownames(annotation_col),colnames(combat_edata.filter)) #[1] TRUE
for (i in 1:6) {
  subtype <- paste0("TC", i)  
  pos <- which(annotation_col$Subtype == subtype)  
  annotation_col_subset <- annotation_col[pos, ]  
  combat_edata_subset <- combat_edata.filter[, pos]  
  assign(paste0("annotation_col_", subtype), annotation_col_subset)  
  assign(paste0("combat_edata_", subtype), combat_edata_subset)
  }

######################################################################
### Run DESeq2 analysis pipeline
###################################################################### 
#==========================

## Using DESEQ2 to find subtype specific differential genes 
#TC1 vs Others
TC1VSOthers<-round(cbind(combat_edata_TC2,combat_edata_TC3,combat_edata_TC4,
                  combat_edata_TC5,combat_edata_TC6,combat_edata_TC1))
conditionTC1VSOthers <- factor(c(rep("Others",ncol(cbind(combat_edata_TC2,combat_edata_TC3,combat_edata_TC4,combat_edata_TC5,combat_edata_TC6))),
                                 rep("TC1",ncol(combat_edata_TC1))))

#TC2 vs Others
TC2VSOthers<-round(cbind(combat_edata_TC1,combat_edata_TC3,combat_edata_TC4,
                        combat_edata_TC5,combat_edata_TC6,combat_edata_TC2))
conditionTC2VSOthers <- factor(c(rep("Others",ncol(cbind(combat_edata_TC1,combat_edata_TC3,combat_edata_TC4,combat_edata_TC5,combat_edata_TC6))),
                                 rep("TC2",ncol(combat_edata_TC2))))

#TC3 vs Others
TC3VSOthers<-round(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC4,
                        combat_edata_TC5,combat_edata_TC6,combat_edata_TC3))
conditionTC3VSOthers <- factor(c(rep("Others",ncol(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC4,combat_edata_TC5,combat_edata_TC6))),
                                 rep("TC3",ncol(combat_edata_TC3))))

#TC4 vs Others
TC4VSOthers<-round(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC3,
                        combat_edata_TC5,combat_edata_TC6,combat_edata_TC4))
conditionTC4VSOthers <- factor(c(rep("Others",ncol(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC3,combat_edata_TC5,combat_edata_TC6))),
                              rep("TC4",ncol(combat_edata_TC4))))

#TC5 vs Others
TC5VSOthers<-round(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC3,
                         combat_edata_TC4,combat_edata_TC6,combat_edata_TC5))
conditionTC5VSOthers <- factor(c(rep("Others",ncol(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC3,combat_edata_TC4,combat_edata_TC6))),
                                 rep("TC5",ncol(combat_edata_TC5))))

#TC6 vs Others
TC6VSOthers<-round(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC3,
                         combat_edata_TC4,combat_edata_TC5,combat_edata_TC6))
conditionTC6VSOthers <- factor(c(rep("Others",ncol(cbind(combat_edata_TC1,combat_edata_TC2,combat_edata_TC3,combat_edata_TC4,combat_edata_TC5))),
                                 rep("TC6",ncol(combat_edata_TC6))))

# Iterate from TC1 to TC6
tc_groups <- paste0("TC", 1:6)

for (tc in tc_groups) {
  
  count_data <- get(paste0(tc, "VSOthers"))
  condition_col <- get(paste0("condition", tc, "VSOthers"))
  
  # 1. Create colData
  colData <- data.frame(
    sample_id = colnames(count_data),  
    condition_col = condition_col,
    stringsAsFactors = FALSE
  )
  colData$condition_col <- factor(
    colData$condition_col, levels = c("Others", tc))
  
  # 3. Run DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = colData,
    design = ~ condition_col
  )
  dds <- DESeq(dds)
  
  # 4. Get results
  res <- results(dds)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # 5. Get gene ID mapping
  gene_ids <- getBM(
    attributes = c("hgnc_symbol", "entrezgene_id"),
    filters = "hgnc_symbol",
    values = rownames(res_df),
    mart = mart
  )
  
  # 6. Create dz_signature by merging results with gene IDs
  dz_signature <- merge(res_df, gene_ids, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
  dz_signature <- subset(dz_signature, !is.na(entrezgene_id) & entrezgene_id != '?')
  
  # Save the dz_signature as DESeq2 results file
  output_file2 <- paste0("./results/", tc, "_DESeq2result.txt")
  write.table(dz_signature, output_file2, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # 7. Create and save MAGMA input file
  genecovar <- dz_signature %>%
    dplyr::select(entrezgene_id, log2FoldChange, pvalue, padj) %>%
    dplyr::mutate(
      log2FoldChange = ifelse(is.na(padj), 0, log2FoldChange),
      log2fc = abs(log2FoldChange),
      pval = -log10(pvalue)
    ) %>%
    dplyr::select(entrezgene_id, log2fc, pval) %>%
    dplyr::distinct(entrezgene_id, .keep_all = TRUE)
  
  output_file <- paste0("./results/", tc, "_genecovar.txt")
  write.table(genecovar, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # 8. Save DESeq2 object
  save(dds, file = paste0("./results/dds", tc, "VSOthers.RData"))
  
  # 9. Clean up all iteration objects
  rm(dds, res, res_df, gene_ids, dz_signature, genecovar, condition_values, colData)
  #gc()  # Garbage collection to free up memory
}


#========================== Analysis of differentially expressed genes for each of the six OA subtypes versus the sum of all other subtypes ==========================
# Initialize lists to store DEGs
subtype_specific_DEGs_up <- list()
subtype_specific_DEGs_down <- list()
subtype_specific_DEGs_total <- list()

for (i in 1:6) {
  # Read the DESeq2 results file
  file_name <- paste0("./results/TC", i, "_DESeq2result.txt")
  DESeq2_result <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  # Remove duplicate genes
  DESeq2_result <- distinct(DESeq2_result, gene, .keep_all = TRUE)
  # Filter for significant DEGs (p < 0.05 & |log2FC| > 0.58 & finite log2FC)
  subtype_specific_DEGs <- subset(DESeq2_result, pvalue < 0.05 & abs(log2FoldChange) > 0.58 & is.finite(log2FoldChange))
  # Classify up/down-regulated genes
  subtype_specific_DEGs$up_down <- ifelse(subtype_specific_DEGs$log2FoldChange > 0, "up", "down")
  # Extract up-regulated genes
  subtype_specific_DEGs_up[[paste0("TC", i, "-Up")]] <- subset(subtype_specific_DEGs, up_down == "up")$gene
  # Extract down-regulated genes
  subtype_specific_DEGs_down[[paste0("TC", i, "-Down")]] <- subset(subtype_specific_DEGs, up_down == "down")$gene
  # Extract differential expressed genes
  subtype_specific_DEGs_total[[paste0("TC", i, "_DEGs")]] <- subset(subtype_specific_DEGs)$gene
  rm(DESeq2_result)
  rm(subtype_specific_DEGs)
}

#Draw venn diagrams showing the common or unique targets of up/down DEGs.
library(ggVennDiagram)
library(VennDiagram)
library(ggvenn)
library(patchwork) 
library(UpSetR)
library(dplyr)

# Convert upregulated DEGs to a binary matrix
palette <- c("#2c9061", "#147ab0", "#cf0e41","#cdf101ff", "#dc7320", "#6610f2")
upset_DEGs_up <- subtype_specific_DEGs_up %>%
  lapply(function(x) data.frame(gene = x)) %>%
  bind_rows(.id = "TC") %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = TC, values_from = present, values_fill = 0) %>%
  as.data.frame()

# Plot UpSet for upregulated DEGs
UpSetR::upset(
  upset_DEGs_up,
  nsets = 6,
  nintersects = 200,  
  mainbar.y.label = " ",
  sets.x.label = " ",
  sets = names(palette),
  sets.bar.color = palette,
  main.bar.color = "#E44C4C",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(2, 2, 1.8, 1.8, 2, 2),
  point.size = 4,
  line.size = 1.5
)#ylabel:Number of Shared upregulated DEGs, xlabel:Upregulated DEG counts per subtype


# Plot UpSet for downregulated DEGs 
upset_DEGs_down <- subtype_specific_DEGs_down %>%
  lapply(function(x) data.frame(gene = x)) %>%
  bind_rows(.id = "TC") %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = TC, values_from = present, values_fill = 0) %>%
  as.data.frame()

# Plot UpSet for downregulated DEGs 
UpSetR::upset(
  upset_DEGs_down,
  nsets = 6,
  nintersects = 200,  
  mainbar.y.label = " ",
  sets.x.label = " ",
  sets = names(palette),
  sets.bar.color = palette,
  main.bar.color = "#3d70bd",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(2, 2, 1.8, 1.8, 2, 2),
  point.size = 4,
  line.size = 1.5
)


