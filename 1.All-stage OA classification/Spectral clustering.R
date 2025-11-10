# Load necessary library
library(ggplot2)
library(genefilter)
library(preprocessCore)
library(GSVA)
library(impute)
library(GSEABase)

#Load datasets
RNAseq1 <- read.table("./Batch/RNAseq1/expressions.tsv.gz",header = TRUE , row.names =1)
RNAseq2 <- read.table("./Batch/RNAseq2/expressions.tsv.gz",header = TRUE , row.names =1)
BRBseq1 <- read.table("./Batch/BRBseq1/expressions.tsv.gz",header = TRUE , row.names =1)
BRBseq2 <- read.table("./Batch/BRBseq2/expressions.tsv.gz",header = TRUE , row.names =1)
# BRBseq3 <- read.table("./Batch/BRBseq3/expressions.tsv.gz",header = TRUE , row.names =1)

var_cutoff <- 0.25
RNAseq1.filter <- varFilter(as.matrix(RNAseq1), var.cutoff = 0.25,filterByQuantile=TRUE)
RNAseq2.filter <- varFilter(as.matrix(RNAseq2), var.cutoff = 0.25,filterByQuantile=TRUE)
BRBseq1.filter <- varFilter(as.matrix(BRBseq1), var.cutoff = 0.25,filterByQuantile=TRUE)
BRBseq2.filter <- varFilter(as.matrix(BRBseq2), var.cutoff = 0.25,filterByQuantile=TRUE)
# BRBseq3.filter <- varFilter(as.matrix(BRBseq3), var.cutoff = 0.25,filterByQuantile=TRUE)

# Store the four dataframes in a named list
dataframes <- list(
  BRBseq1 = BRBseq1.filter,
  BRBseq2 = BRBseq2.filter,
  # BRBseq3 = BRBseq3.filter,
  RNAseq1 = RNAseq1.filter,
  RNAseq2 = RNAseq2.filter
)

# Store the gene set lists in a named list
CP.geneSets = getGmt("./ssGSEA/c2.cp.v2024.1.Hs.symbols.gmt")
H.geneSets = getGmt("./ssGSEA/h.all.v2024.1.Hs.symbols.gmt")
GO.geneSets = getGmt("./ssGSEA/c5.go.bp.v2024.1.Hs.symbols.gmt")
geneSetLists <- list(
  CP = CP.geneSets,
  H = H.geneSets,
  GO = GO.geneSets
)

# Initialize a list to store all enrichment scores
allEnrichmentScores <- list()

# Iterate over each dataset
for (dfName in names(dataframes)) {
  currentDF <- dataframes[[dfName]]
  enrichmentScoresList <- list()
  for (gsName in names(geneSetLists)) {
    scores <- gsva(
      expr = currentDF,                   # Current dataframe
      gset.idx.list = geneSetLists[[gsName]], # Current gene set list
      method = "ssgsea",                 # ssGSEA method
      ssgsea.norm = TRUE                 # Normalize ssGSEA scores
    )
    enrichmentScoresList[[gsName]] <- scores
  }
  allEnrichmentScores[[dfName]] <- enrichmentScoresList
}

## Combine all enrichment scores into one dataframe
RNAseq1_finalCombined <- do.call(rbind, allEnrichmentScores[["RNAseq1"]])
RNAseq2_finalCombined <- do.call(rbind, allEnrichmentScores[["RNAseq2"]])
BRBseq1_finalCombined <- do.call(rbind, allEnrichmentScores[["BRBseq1"]])
BRBseq2_finalCombined <- do.call(rbind, allEnrichmentScores[["BRBseq2"]])

overlap <- Reduce(intersect, list(rownames(RNAseq1_finalCombined),
                                  rownames(RNAseq2_finalCombined),
                                  rownames(BRBseq1_finalCombined),
                                  rownames(BRBseq2_finalCombined)))
RNAseq1_Combined_overlap <- RNAseq1_finalCombined[overlap,]
RNAseq2_Combined_overlap <- RNAseq2_finalCombined[overlap,]
BRBseq1_Combined_overlap <- BRBseq1_finalCombined[overlap,]
BRBseq2_Combined_overlap <- BRBseq2_finalCombined[overlap,]
# BRBseq3_Combined_overlap <- BRBseq3_finalCombined[overlap,]

matrix <- cbind(RNAseq1_Combined_overlap,RNAseq2_Combined_overlap,
                BRBseq1_Combined_overlap,BRBseq2_Combined_overlap)

#Filter by Quality check result
QC <- read.table("./Quality_check_results.tsv",header = TRUE )
identical(QC$Sample_ID,colnames(matrix)) #[1] TRUE
matrix <- matrix[,!(QC$Final_Status) %in% c("Failed") ] 
QC <- QC[!(QC$Final_Status) %in% c("Failed") ,]

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

##Remove batch effect
library(sva)
batch <- c(rep('RNAseq1',87),
           rep('RNAseq2',14),
           rep('BRBseq1',20),
           rep('BRBseq2',18))

# Using ComBat to remove batch effect
combat_edata = ComBat(dat = matrix, batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
write.table(combat_edata, './SpectralClustering/combat_edata.txt', col.names = NA, sep = '\t', quote = FALSE)

# Using the UMAP plot to verify that the batch effect has been removed
umap_df2 <- umap(t(combat_edata))
umap_df2 = cbind.data.frame(umap_df2$layout,QC$Batch)
colnames(umap_df2) = c("UMAP1","UMAP2","Batch")
umap_df2$Batch = as.factor(umap_df2$Batch)
umap_df2$Batch <- factor(umap_df2$Batch, levels = c("BRBseq1","BRBseq2","RNAseq1","RNAseq2"))
ggplot(umap_df2,aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Batch),size=4) + theme_classic(base_size = 20)+scale_color_manual(values=c("BRBseq1"="#00CED1","BRBseq2"="#0000CD","RNAseq1"="#984EA3","RNAseq2"="#FF0000")) + guides(colour = guide_legend(override.aes = list(size=4)))

#==========================
## Determining optimal number of clusters
library(kernlab)  # For spectral clustering
library(ggplot2)  # For plotting
library(rstatix)  # For Kruskal-Wallis test

# Function to calculate significant features using Kruskal-Wallis
calculate_kruskal_features <- function(data, k) {
  clusters <- specc(data, centers = k)@.Data
  p_values <- numeric(ncol(data))
  for (feature in 1:ncol(data)) {
    test_result <- kruskal.test(data[, feature] ~ as.factor(clusters))
    p_values[feature] <- test_result$p.value
  }
  p_adj <- p.adjust(p_values, method = "bonferroni")
  sum(p_adj < 0.05)
}

# Test k from 2 to 14
k_values <- 2:14
significant_features <- sapply(k_values, function(k) {
  calculate_kruskal_features(t(combat_edata), k)
})

# Create results table
Cluster_number_evaluation <- data.frame(
  Clusters = k_values,
  Significant_Features = significant_features)

# Visualize the Results of Cluster Number Evaluation
ggplot(Cluster_number_evaluation, 
       aes(x = Clusters, y = Significant_Features)) +
  geom_col(aes(fill = ifelse(Clusters == 6, "6", "Other")),  # Conditional fill
           width = 0.6, alpha = 0.9) +
  scale_fill_manual(values = c("6" = "#D44C3C", "Other" = "#999696"), 
                    guide = "none") +  # Apply colors and hide legend
  labs(
    x = "Number of clusters",
    y = "Number of features"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

#==========================
# Using spectral clustering 
specc <- specc(t(combat_edata), centers = 6)
ht_ann <- specc@.Data
annotation_col = data.frame(row.names = colnames(combat_edata),Subtype = ht_ann)
annotation_col2 <- annotation_col
annotation_col2$Subtype  <- sub("^1","TC2",annotation_col2$Subtype)
annotation_col2$Subtype  <- sub("^2","TC4",annotation_col2$Subtype)
annotation_col2$Subtype  <- sub("^3","TC3",annotation_col2$Subtype)
annotation_col2$Subtype  <- sub("^4","TC6",annotation_col2$Subtype)
annotation_col2$Subtype  <- sub("^5","TC5",annotation_col2$Subtype)
annotation_col2$Subtype  <- sub("^6","TC1",annotation_col2$Subtype)
write.table(annotation_col2, './SpectralClustering/annotation_col.txt', col.names = NA, sep = '\t', quote = FALSE)
umap_Subtype = cbind.data.frame(umap_df2,annotation_col2$Subtype)
colnames(umap_Subtype) = c("UMAP1","UMAP2","Batch","Subtype")
# umap_Subtype$Subtype <- factor(umap_Subtype$Subtype, levels = c("TC1","TC2","TC3","TC4","TC5","TC6"))
write.table(umap_Subtype, './SpectralClustering/umap_Subtype.txt', col.names = NA, sep = '\t', quote = FALSE)
ggplot(umap_Subtype,aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Subtype),size=4) + theme_classic(base_size = 20)+scale_color_manual(values=c("TC1"="#FF0000","TC2"="#984EA3","TC3"="#87CEad","TC4"="#00CED1","TC5"="#EE82EE","TC6"="#0000CD")) + guides(colour = guide_legend(override.aes = list(size=4)))

