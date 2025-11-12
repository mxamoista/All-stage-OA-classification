##### Protein differential analysis #####

## Differential protein analysis: TC4(baseline)/Others(HA treatment) ##

# read protein quantification files
Protein_quantification <- read_excel("./Expressed_annotation.xlsx")
# Select protein quantification values of TC1-6 samples
Protein_quantification_value <- Protein_quantification[, -c(1:3, 33:41)]
rownames(Protein_quantification_value) <- Protein_quantification$Accession
identical(as.character(Information$IndexP),colnames(Protein_quantification_value)) #[1] TRUE

## Using limma to find TC4 specific differential proteins at baseline
library(limma)
Information_TC4 <- Information %>% filter(Points == "Baseline" )
Protein_quantification_value_baseline <- Protein_quantification_value[,colnames(Protein_quantification_value) %in% Information_TC4$IndexP]
rownames(Protein_quantification_value_baseline) <- rownames(Protein_quantification_value)
Information_TC4$TC4 <- ifelse(Information_TC4$Subtype == "TC4", "yes", "no")

# Specify the sample condition (TC4(baseline) vs. TC4(HA treatment))
design_TC4 <- model.matrix(~0 + TC4 , Information_TC4)
fit <- lmFit(Protein_quantification_value_baseline, design_TC4)

## The variability of protein expression is compared between these groups
contr <- makeContrasts(TC4yes - TC4no, levels = design_TC4)
contr

## Estimate contrast for each protein
fit2 <- contrasts.fit(fit, contr)

## Calculate t-stats
## Empirical Bayes smoothing of standard errors (shrinks standard errors
## that are much larger or smaller than those from other proteins towards the average standard error)
fit2 <- eBayes(fit2)

## Extract results
fit_results_TC4 <- as.data.frame(topTable(fit2, sort.by = "P", n = Inf))
fit_results_TC4$Accession <- rownames(fit_results_TC4)

## compile results for plotting
map <- Protein_quantification %>% select(Accession,Gene) %>% distinct() 

results_limma_TC4 <- merge(fit_results_TC4, map, by = "Accession", all.x=T) 
write.table(results_limma_TC4, "./results/results_limma_TC4(baseline VS HA treatment).txt",sep="\t",quote=F,row.names=F,col.names=T)

## Using DESEQ2 to find TC4 specific differential proteins at baseline
library(DESeq2)

TC4VSOthers <- round(cbind(Protein_quantification_value_baseline[ ,which(Information_TC4$TC4 == "no") ],
            Protein_quantification_value_baseline[ ,which(Information_TC4$TC4 == "yes") ]))
rownames(TC4VSOthers) <- rownames(Protein_quantification_value_baseline)
conditionTC4VSOthers <- factor(c(rep("Others",ncol(Protein_quantification_value_baseline[ ,which(Information_TC4$TC4 == "no") ])),
                                 rep("TC4",ncol(Protein_quantification_value_baseline[ ,which(Information_TC4$TC4 == "yes") ]))))
# Create colData
colDataTC4VSOthers <- data.frame(
  sample_id = colnames(TC4VSOthers),  
  conditionTC4VSOthers = conditionTC4VSOthers,
  stringsAsFactors = FALSE
)
colDataTC4VSOthers$conditionTC4VSOthers <- factor(
  colDataTC4VSOthers$conditionTC4VSOthers, levels = c("Others", "TC4"))

# Run DESeq2 analysis
dds <- DESeqDataSetFromMatrix(
  countData = TC4VSOthers,
  colData = colDataTC4VSOthers,
  design = ~ conditionTC4VSOthers
)
dds <- DESeq(dds)

# Get results
res_dds <- results(dds)
res_dds <- as.data.frame(res_dds)
res_dds$Accession <- rownames(res_dds)
results_DEseq2_TC4 <- merge(res_dds, map, by = "Accession", all.x=T) 


# plot results
## using DESEQ2 results 
toplot_TC4 <- results_DEseq2_TC4 %>% select(Gene,log2FoldChange,pvalue)
toplot_TC4$ToLabel = factor(ifelse(toplot_TC4$pvalue < 0.1 & abs(toplot_TC4$log2FoldChange) >= 0.5, ifelse(toplot_TC4$log2FoldChange >= 0.5,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))

ggplot(toplot_TC4) + 
  geom_point(aes(x=log2FoldChange,y=-log10(pvalue), color=ToLabel),alpha=0.5, size=2.5) +
  scale_color_manual(values=c("#d01c8b","#4dac26","grey"))+
  geom_vline(xintercept=c(-0.1,0.1), linewidth=0.3,linetype='dotted', color='gray50') +
  geom_text_repel(data =toplot_TC4[toplot_TC4$Gene %in% c('APOE'),],
            aes(x=log2FoldChange,y=-log10(pvalue), label=Gene), size = 5,
              box.padding = unit(2, "lines"),
              point.padding = unit(2, "lines"), 
              segment.color = "black", show.legend =FALSE) +
  ggtitle("Static/Others") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

## Using DESEQ2 to find TC4 specific differential proteins at HA treatment
Information_TC4_HA <- Information %>% filter(Points != "Baseline" ) %>% filter(IndexP != "25" )
Protein_quantification_value_HA <- Protein_quantification_value[,colnames(Protein_quantification_value) %in% Information_TC4_HA$IndexP]
rownames(Protein_quantification_value_HA) <- rownames(Protein_quantification_value)
Information_TC4_HA$TC4 <- ifelse(Information_TC4_HA$Subtype == "TC4", "yes", "no")

TC4VSOthers_HA <- round(cbind(Protein_quantification_value_HA[ ,which(Information_TC4_HA$TC4 == "no") ],
                              Protein_quantification_value_HA[ ,which(Information_TC4_HA$TC4 == "yes") ]))
rownames(TC4VSOthers_HA) <- rownames(Protein_quantification_value_HA)
conditionTC4VSOthers_HA <- factor(c(rep("Others",ncol(Protein_quantification_value_HA[ ,which(Information_TC4_HA$TC4 == "no") ])),
                                 rep("TC4",ncol(Protein_quantification_value_HA[ ,which(Information_TC4_HA$TC4 == "yes") ]))))
# Create colData
colDataTC4VSOthers_HA <- data.frame(
  sample_id = colnames(TC4VSOthers_HA),  
  conditionTC4VSOthers = conditionTC4VSOthers_HA,
  stringsAsFactors = FALSE
)
colDataTC4VSOthers_HA$conditionTC4VSOthers <- factor(
  colDataTC4VSOthers_HA$conditionTC4VSOthers, levels = c("Others", "TC4"))

# Run DESeq2 analysis
dds_HA <- DESeqDataSetFromMatrix(
  countData = TC4VSOthers_HA,
  colData = colDataTC4VSOthers_HA,
  design = ~ conditionTC4VSOthers
)
dds_HA <- DESeq(dds_HA)

# Get results
res_dds_HA <- results(dds_HA)
res_dds_HA <- as.data.frame(res_dds_HA)
res_dds_HA$Accession <- rownames(res_dds_HA)
results_DEseq2_TC4_HA <- merge(res_dds_HA, map, by = "Accession", all.x=T) 

# plot results
## using DESEQ2 results 
toplot_TC4_HA <- results_DEseq2_TC4_HA %>% select(Gene,log2FoldChange,pvalue)
toplot_TC4_HA$ToLabel = factor(ifelse(toplot_TC4_HA$pvalue < 0.1 & abs(toplot_TC4_HA$log2FoldChange) >= 0.5, ifelse(toplot_TC4_HA$log2FoldChange >= 0.5,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
toplot_TC4_HA <- toplot_TC4_HA[!(toplot_TC4_HA$Gene) == "HYI", ] # exclude outlier
ggplot(toplot_TC4_HA) + 
  geom_point(aes(x=log2FoldChange,y=-log10(pvalue), color=ToLabel),alpha=0.5, size=2.5) +
  scale_color_manual(values=c("#DC2B18","#29419E","grey"))+
  geom_vline(xintercept=c(-0.1,0.1), linewidth=0.3,linetype='dotted', color='gray50') +
  geom_text_repel(data =toplot_TC4_HA[toplot_TC4_HA$Gene %in% c('APOE'),],
                  aes(x=log2FoldChange,y=-log10(pvalue), label=Gene), size = 5,
                  box.padding = unit(2, "lines"),
                  point.padding = unit(2, "lines"), 
                  segment.color = "black", show.legend =FALSE) +
  ggtitle("Static/Others") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

