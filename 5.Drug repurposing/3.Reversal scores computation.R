###########################################################################
# Draw the heatmap showing drug candidates that "flips" the signature of each subtype -----------------------------------------
###########################################################################

# Define the subtype names
subtypes <- paste0("TC", 1:6)

Drugs_with_negative_cMap_scores<- list()

# Store the drug candidates with negative Cmap scores for each subtype
for (subtype in subtypes) {
  # Get the data for current subtype
  subtype_data <- get(paste0(subtype, "_drug_instances_all"))
  subtype_data <- subtype_data[subtype_data$cell_id == "PC3", ]
  
  # Merge with OA drugs
  label_compounds <- subtype_data[subtype_data$drug_name %in% OA_DRUG$DrugBank.ID, ]
  label_compounds <- merge(label_compounds, OA_DRUG, by.x = "drug_name", by.y = "DrugBank.ID")
  label_compounds <- subset(label_compounds, cmap_score < 0)
  
  #Select the columns of drug_name and cmap_score
  label_compounds <- label_compounds[, c("Drug", "cmap_score")]
  colnames(label_compounds)[2] <- paste0("cmap_score_", subtype)
  
  Drugs_with_negative_cMap_scores[[subtype]] <- label_compounds
  rm(label_compounds)
}

# Merge the drug candidates among six subtypes into one dataframe
all_drugs <- data.frame(drug_name = unique(unlist(lapply(Drugs_with_negative_cMap_scores, function(x) x$Drug))))
merged_negative_cMap_scores <- all_drugs

for (subtype in subtypes) {
  subtype_df <- Drugs_with_negative_cMap_scores[[subtype]]
  colnames(subtype_df)[colnames(subtype_df) == "Drug"] <- "drug_name"
  merged_negative_cMap_scores <- merge(merged_negative_cMap_scores, subtype_df, by = "drug_name", all = TRUE)
  rm(subtype_df)
}
rm(all_drugs)

# Draw heatmap for cMap scores
library("pheatmap")
library(ComplexHeatmap)
library(circlize)
rownames(merged_negative_cMap_scores) <- merged_negative_cMap_scores$drug_name
merged_negative_cMap_scores[,-1][is.na(merged_negative_cMap_scores[,-1])] <- 0

colPal <- colorRamp2(
  breaks = c(min(merged_negative_cMap_scores[,-1]), max(merged_negative_cMap_scores[,-1])),
  colors = c("#D44C3C","white")  
)

Heatmap(merged_negative_cMap_scores[,-1],cluster_columns = FALSE,
        cluster_rows = FALSE,height = unit(8, "cm"),width = unit(4, "cm"),
        col = colPal,border_gp = gpar(col = "black"))

###########################################################################
# Draw the heatmap showing reversal relationship between the drug candidates with negative cMap scores and the subtype-specific profiles -----------------------------------------
###########################################################################
library(dplyr)
all_drugs <- data.frame(Drug = unique(c(TC1_cor_results[TC1_cor_results$Correlation < 0,]$Drug, 
                                        TC2_cor_results[TC2_cor_results$Correlation < 0,]$Drug, 
                                        TC3_cor_results[TC3_cor_results$Correlation < 0,]$Drug,
                                        TC4_cor_results[TC4_cor_results$Correlation < 0,]$Drug,
                                        TC5_cor_results[TC5_cor_results$Correlation < 0,]$Drug, 
                                        TC6_cor_results[TC6_cor_results$Correlation < 0,]$Drug)))

#Combine all drug candidates that have negative correlation scores with any subtype into a single dataframe
merged_negative_correlation_scores <- all_drugs
for (subtype in subtypes) {
  subtype_df <- get(paste0(subtype, "_cor_results"))
  colnames(subtype_df)[colnames(subtype_df) == "Correlation"] <- paste0(subtype, "_Correlation")
  merged_negative_correlation_scores <- merge(merged_negative_correlation_scores, subtype_df[subtype_df[[paste0(subtype, "_Correlation")]] < 0, ][,-1], by = "Drug", all = TRUE)
  rm(subtype_df)
}
rm(all_drugs)


# Draw heatmap for reversal relationship
library("pheatmap")
library(ComplexHeatmap)
library(circlize)
rownames(merged_negative_correlation_scores) <- merged_negative_correlation_scores$Drug
merged_negative_correlation_scores[,-1][is.na(merged_negative_correlation_scores[,-1])] <- 0

colPal2 <- colorRamp2(
  breaks = c(min(merged_negative_correlation_scores[,-1]), max(merged_negative_correlation_scores[,-1])),
  colors = c("#D44C3C","white")  
)

Heatmap(merged_negative_correlation_scores[,-1],cluster_columns = FALSE,
        cluster_rows = FALSE,height = unit(8, "cm"),width = unit(4, "cm"),
        col = colPal2,border_gp = gpar(col = "black"))

###########################################################################
# Generate a heatmap illustrating the reversal effects of HA on TC4 and of Vitamin D3 on TC5 -----------------------------------------
###########################################################################

# reversal effects of HA on TC4
HA_TC4_signature <- merge(
  TC4_dz_signature[, c("Symbol", "value")],
  drug_sig[, as.character(which(as.data.frame(drug_list) == OA_DRUG[which(OA_DRUG$Drug == "Hyaluronic acid"),]$DrugBank.ID)), drop = FALSE],
  by.x = "Symbol",
  by.y = 0
)
HA_TC4_signature <- HA_TC4_signature[order(HA_TC4_signature$value, decreasing = TRUE), ]
# HA_TC4_signature <- HA_TC4_signature[, -1]
for (i in 1:ncol(HA_TC4_signature[, -1])) {
  HA_TC4_signature[, -1][, i] <- rank(HA_TC4_signature[, -1][, i])
}# Rank normalize each column, except the "Symbol" column

# Generate Heatmap to show reversal effects of HA on TC4
tiff_filename <- paste0("./picture/Reversal heatmap of HA and VitaminD3/The reversal effects of HA on TC4.tiff")

tiff(
  filename = tiff_filename,
  width = 4, 
  height = 8, 
  units = "in",
  res = 600, 
  compression = "lzw"
)

# Set plot margins
par(mar = c(8, 4, 4, 2))  # Bottom, Left, Top, Right

# Create color palette
colPal <- colorRampPalette(c("red", "white", "#000099"))(nrow(HA_TC4_signature[, -1]))

# Calculate x-axis positions
n_cols <- ncol(HA_TC4_signature[, -1])
x_positions <- seq(0, 1, length.out = max(2, n_cols))

# Plot heatmap
image(t(HA_TC4_signature[, -1]), col = colPal, axes = FALSE)

# Prepare and add labels
drug_labels <- "Hyaluronic acid"
text(
  x = x_positions,
  y = -0.1,
  labels = c("TC4", drug_labels),
  srt = 45,
  pos = 2,
  offset = 1,
  xpd = TRUE,
  cex = 1.2
)

dev.off()
message("Successfully saved: ", tiff_filename)

# reversal effects of VitaminD3 on TC5
VitaminD3_TC5_signature <- merge(
  TC5_dz_signature[, c("Symbol", "value")],
  drug_sig[, as.character(which(as.data.frame(drug_list) == OA_DRUG[which(OA_DRUG$Drug == "Vitamin D3"),]$DrugBank.ID)), drop = FALSE],
  by.x = "Symbol",
  by.y = 0
)
VitaminD3_TC5_signature <- VitaminD3_TC5_signature[order(VitaminD3_TC5_signature$value, decreasing = TRUE), ]
# VitaminD3_TC5_signature <- VitaminD3_TC5_signature[, -1]
for (i in 1:ncol(VitaminD3_TC5_signature[, -1])) {
  VitaminD3_TC5_signature[, -1][, i] <- rank(VitaminD3_TC5_signature[, -1][, i])
}# Rank normalize each column, except the "Symbol" column

# Generate Heatmap to show reversal effects of HA on TC4
tiff_filename <- paste0("./picture/Reversal heatmap of HA and VitaminD3/The reversal effects of VitaminD3 on TC5.tiff")

tiff(
  filename = tiff_filename,
  width = 4, 
  height = 8, 
  units = "in",
  res = 600, 
  compression = "lzw"
)

# Set plot margins
par(mar = c(8, 4, 4, 2))  # Bottom, Left, Top, Right

# Create color palette
colPal <- colorRampPalette(c("red", "white", "#000099"))(nrow(VitaminD3_TC5_signature[, -1]))

# Calculate x-axis positions
n_cols <- ncol(VitaminD3_TC5_signature[, -1])
x_positions <- seq(0, 1, length.out = max(2, n_cols))

# Plot heatmap
image(t(VitaminD3_TC5_signature[, -1]), col = colPal, axes = FALSE)

# Prepare and add labels
drug_labels <- "Vitamin D3"
text(
  x = x_positions,
  y = -0.1,
  labels = c("TC5", drug_labels),
  srt = 45,
  pos = 2,
  offset = 1,
  xpd = TRUE,
  cex = 1.2
)

dev.off()
message("Successfully saved: ", tiff_filename)

###########################################################################
# Perform GSEA-Preranked in HA_TC4_signature and VitaminD3_TC5_signature -----------------------------------------
###########################################################################
library(fgsea)          
library(data.table)     
library(ggplot2)        
library(dplyr)          
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)

# Save the table for Cytoscape
write.table(HA_TC4_signature,"./result/HA_TC4_signature.txt",sep = '\t')
write.table(HA_TC4_signature[,-2],"./Cytoscape/HA_signature.txt",sep = '\t',row.names = F, quote = FALSE)
write.table(HA_TC4_signature[,-3],"./Cytoscape/TC4_signature.txt",sep = '\t',row.names = F, quote = FALSE)

# Transform Symbol into Entrzid
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl") # Connect to the human genes dataset
# Get the mapping
gene_ids_HA_TC4 <- getBM(attributes=c("hgnc_symbol", "entrezgene_id"),
                      filters="hgnc_symbol",
                      values=HA_TC4_signature$Symbol,
                      mart=mart)
# Merge with mapping result
HA_TC4_signature_entrezid <- merge(HA_TC4_signature, gene_ids_HA_TC4, by.x="Symbol", by.y="hgnc_symbol", all.x=TRUE)
c2.cp_gmt <- read.gmt("./GSEA/c2.cp.v2024.1.Hs.symbols.gmt") 
c5.go.bp_gmt <- read.gmt("./GSEA/c5.go.bp.v2024.1.Hs.symbols.gmt") 

geneList_TC4 <- HA_TC4_signature[order(HA_TC4_signature$value, decreasing = TRUE), ]$value
names(geneList_TC4) <- HA_TC4_signature[order(HA_TC4_signature$log2FoldChange, decreasing = TRUE), ]$Symbol

gsea_TC4 <- GSEA(geneList_TC4, TERM2GENE = c2.cp_gmt, by = "fgsea") 


