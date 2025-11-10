# Purpose: Generate boxplots show ssGSEA scores related to OA pathology,drivers or symptoms among subtypes 
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(ggpubr)
library(tidyr)
library(tibble) 
###########################################################################
#  Step 0:Load gene_sets associated with OA pathology,drivers or symptoms -------------------------------------------------------
###########################################################################
gene_sets <- c("KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_CHONDROITIN_SULFATE",
               "REACTOME_GLYCOSAMINOGLYCAN_METABOLISM",
               "KEGG_GLYCOSAMINOGLYCAN_DEGRADATION",
               "GOBP_GLYCOSAMINOGLYCAN_CATABOLIC_PROCESS",
               "WP_GLYCOSAMINOGLYCAN_DEGRADATION",
               "WP_MATRIX_METALLOPROTEINASES",
               "REACTOME_ACTIVATION_OF_MATRIX_METALLOPROTEINASES",
               "REACTOME_COLLAGEN_DEGRADATION",
               "REACTOME_COLLAGEN_FORMATION",
               "GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS",
               "GOBP_COLLAGEN_CATABOLIC_PROCESS",
               "GOBP_COLLAGEN_METABOLIC_PROCESS",
               "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
               "GOBP_POSITIVE_REGULATION_OF_COLLAGEN_METABOLIC_PROCESS",
               "WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE",
               "WP_INFLAMMATORY_RESPONSE_PATHWAY",
               "GOBP_ORGAN_OR_TISSUE_SPECIFIC_IMMUNE_RESPONSE",
               "GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION",
               "GOBP_POSITIVE_REGULATION_OF_COMPLEMENT_ACTIVATION",
               "GOBP_REGULATION_OF_NERVOUS_SYSTEM_PROCESS",
               "GOBP_NEUROGENESIS",
               "REACTOME_NERVOUS_SYSTEM_DEVELOPMENT",
               "GOBP_POSITIVE_REGULATION_OF_NEUROGENESIS",
               "GOBP_REGULATION_OF_NEUROGENESIS",
               "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION",
               "GOBP_POSITIVE_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT",
               "GOBP_REGULATION_OF_NERVOUS_SYSTEM_DEVELOPMENT",
               "REACTOME_NERVOUS_SYSTEM_DEVELOPMENT",
               "GOBP_GENERATION_OF_NEURONS",
               "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
               "REACTOME_FATTY_ACIDS")

palette <- c("#2c9061", "#147ab0", "#cf0e41","#cdf101ff", "#dc7320", "#6610f2")

###########################################################################
#  Step 1:Draw boxplots towards gene_sets among subtypes -------------------------------------------------------
###########################################################################
for (gs in gene_sets) {
  df_plot <- data.frame(
    Score = t(combat_edata)[, gs],
    group = annotation_col2$Subtype
  )
  
  p <- ggplot(df_plot, aes(x = group, y = Score,
                           fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.8,
                 position = position_dodge(width = 0.75)) +
    geom_jitter(shape = 16, size = 1.5, alpha = 0.8,
                position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = palette) +
    scale_colour_manual(values = palette) +
    labs(title = gs, x = NULL, y = "Normalized ssGSEA score") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold"),
          axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
          axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
          plot.title = element_text(family = "Helvetica", size = 14, face = "bold", hjust = 0.5),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.line.x = element_line())
  
  ggsave(filename = file.path("./Picture/Boxplot", paste0(gs, ".tiff")),
         plot = p, width = 5, height = 4, dpi = 300)
}

###########################################################################
#  Step 2:Wilcoxon rank-sum test to calculate significance -------------------------------------------------------------
###########################################################################
Calculate_significance <- function(vec, grp, gs_name) {
  lvls <- levels(grp)
  n    <- length(lvls)
  mat <- matrix(NA_real_, n, n,dimnames = list(lvls, lvls))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- wilcox.test(vec[grp == lvls[i]],
                               vec[grp == lvls[j]])$p.value
    }
  }
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  diag(mat) <- 1
  tidy_df <- mat %>%
    as.data.frame() %>%
    rownames_to_column("Row") %>%
    pivot_longer(-Row, names_to = "Col", values_to = "p")
  ggplot(tidy_df,
         aes(Row, Col, fill = p)) +
    geom_tile(colour = "NA") +
    geom_text(aes(label = sprintf("%.2g", p)), size = 3.5) +
    scale_fill_gradientn(
      colours = c("#dc7320","#e89d62","#edb385","#f5d3b8","#fff0e5","white","#f0f6fc","#cfe2f6","#9ec3ed","#7dafe6","#147ab0"),
      name    = " "
    ) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    labs(title = gs_name) +
    theme_minimal() +
    theme(
      panel.grid   = element_blank(),
      axis.title   = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(face = "bold"),
      axis.text    = element_text(face = "bold"),
      plot.title   = element_text(hjust = .5, face = "bold")
    )
}

# Iterate gene_set -----------------------------------------------------
plots <- list()
for (gs in gene_sets) {
  plots[[gs]] <- Calculate_significance(t(combat_edata)[, gs], as.factor(annotation_col2$Subtype), gs)
}

for (gs in names(plots)){
  ggsave(filename = file.path("./Picture/Heatmap", paste0(gs, ".tiff")),plot = plots[[gs]], width = 5, height = 4, dpi = 300)
}
