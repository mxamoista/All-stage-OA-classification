library(tidyverse)
library(ggplot2)
library(scales)
library(tibble)
library(dplyr)
library(tidyr)

Protein_quantification_filtered <- Protein_quantification[, -c(2:3, 33:41)]
Protein_quantification_filtered2 <- as.data.frame(t(Protein_quantification_filtered))
colnames(Protein_quantification_filtered2) <- Protein_quantification_filtered2["Accession",]
Protein_quantification_filtered2 <- Protein_quantification_filtered2[-which(rownames(Protein_quantification_filtered2) == "Accession"), ]

#Delete the proteomics data of patients 12 months after HA treatment, indexed as P_74
Protein_quantification_filtered2 <- Protein_quantification_filtered2[!rownames(Protein_quantification_filtered2) == "25", ]
Information2 <- Information[match(rownames(Protein_quantification_filtered2) , Information$IndexP), ]
Information2$Group = factor(c(Information2$Group), levels = c("no", "HA"))

###########################################################################
# Draw the change of potential molecular mediators in knee osteoarthritis pathology between "baseline" and "HA treatment"  -----------------------------------------
###########################################################################
#Select potential molecular mediators of knee osteoarthritis pathology in synovial fluid
Inflammation = c("TD02", "MIF", "CD74", "IL-17A", "IL-17","IL-6","TNF","IL-40")
Cartilage_degeneration = c("MMP2", "RANKL", "SPARC", "APOE")
Pain = c("IL-6","IL-17A","COMP","PIICP")
OA_symptoms = c("MMP3","sVCAM1","sICAM1","TIMP1","VEGF","MCP1")
HA_receptor = c("CD44")
HA_binder = c("PRG4")
upregulated_proteins_in_OASF = c("CCL26","ANGPTL7","PROK1","CD27","VEGFD","CXCL12",
                                 "EGF","CSF3R","IL1RAPL2","TNFAIP6","CX3CL1","CXCL13",
                                 "CCL24","BMP15","EDAR","MMP25","KIT","NRP2",
                                 "PLAUR","CTNNB1")

molecular_mediators <- c(Inflammation, Cartilage_degeneration, Pain, OA_symptoms,HA_receptor,HA_binder,upregulated_proteins_in_OASF)

subtypes <- unique(Information2$Subtype)

#Draw the change between "baseline" and "HA treatment"
for (subtype in subtypes) {
  dir_path <- file.path("./picture/change", subtype)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message(paste("Created directory:", dir_path))
  }
}

for (feature in molecular_mediators) {
  if (!feature %in% map$Gene) {
    message(paste("Feature", feature, "not found, skipping..."))
    next
  }
  accessions <- map[which(map$Gene == feature), ]$Accession
  for (accession in accessions) {
    if (!accession %in% colnames(Protein_quantification_filtered2)) {
      message(paste("Accession", accession, "not found in data, skipping..."))
      next
    }
    plot_data <- data.frame(
      Quantification = Protein_quantification_filtered2[[accession]],
      Group = Information2$Group,
      Patient = Information2$IndexT,  
      Subtype = Information2$Subtype  
    )
    plot_data$Quantification <- as.numeric(plot_data$Quantification)
    
    for (highlight_subtype in subtypes) {
      subtype_colors <- setNames(rep("gray", length(subtypes)), subtypes)
      subtype_colors[highlight_subtype] <- "#cf0e41"  
      
      p <- ggplot(plot_data, aes(x = Group, y = Quantification, group = Patient)) +
        geom_line(aes(color = Subtype), linewidth = 1, alpha = 0.7) +
        geom_point(aes(color = Subtype, shape = Subtype), size = 3) +
        scale_y_continuous(labels = function(x) format(x, digits = 3, scientific = FALSE)) +
        scale_x_discrete(expand = expansion(mult = 0.1))+
        labs(
          title = paste("Expression of", feature, "(", accession, ")"),
          subtitle = paste("Highlighted:", highlight_subtype),  
          x = "Group",
          y = " ",
          color = "Subtype",
          shape = "Subtype"
        ) +
        scale_color_manual(values = subtype_colors) +
        theme_minimal() +
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          legend.position = "right",
          axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
          axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),
          axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
          plot.title = element_text(family = "Helvetica", size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(family = "Helvetica", size = 12, face = "bold", hjust = 0.5),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.line.x = element_line() 
        )
      
      # pdf_file <- file.path("./picture/change", highlight_subtype, paste0(feature, "_", accession, ".pdf"))
      # pdf(file = pdf_file, width = 5, height = 6)
      # print(p)
      # dev.off()
      png_file <- file.path("./picture/change", highlight_subtype, paste0(feature, "_", accession, ".png"))
      png(file = png_file, width = 1000, height = 1200,res = 300)
      print(p)
      dev.off()
      
      message(paste("Plot for", feature, "(", accession, ") in subtype", highlight_subtype, "saved successfully"))
    }
  }
}

#Draw the change value between "baseline" and "HA treatment"
#Select the protein ,the accession with "H0YDX6","P02649","P10721","J3KP74"
#Select group "TC1-TC3","TC4","TC5","TC6"

selected_accessions <- c("H0YDX6", "P02649","Q92954", "J3KP74")

###########################################################################
#  Function to extract and merge protein quantification data with patient information -----------------------------------------
###########################################################################
extract_protein_data <- function(protein_data, info_data, accessions) {
  # Convert row names (IndexP) to a column and select specified accessions
  protein_data_with_id <- protein_data %>%
    rownames_to_column("IndexP") %>%
    mutate(IndexP = as.character(IndexP)) %>%
    mutate(across(all_of(accessions), as.numeric)) %>%
    select(IndexP, all_of(accessions))
  
  # Merge with patient information and calculate mean values per patient per group
  merged_data <- info_data %>%
    mutate(IndexP = as.character(IndexP)) %>%
    select(IndexP, IndexT, Subtype, Group) %>%
    left_join(protein_data_with_id, by = "IndexP") %>%
    # Group by patient (IndexT) and experimental group to handle multiple measurements
    group_by(IndexT, Group) %>%
    # Calculate mean protein quantification for each accession per patient per group
    summarise(across(all_of(accessions), mean, na.rm = TRUE), .groups = "drop") %>%
    # Add subtype information (taking the first subtype record for each patient)
    left_join(
      info_data %>% 
        select(IndexT, Subtype) %>% 
        distinct(IndexT, .keep_all = TRUE),
      by = "IndexT"
    )
  
  return(merged_data)
}
protein_wide <- extract_protein_data(Protein_quantification_filtered2, Information2, selected_accessions)

# Function to calculate differences between "baseline" and "HA treatment" for each accession
calculate_differences <- function(wide_data, accessions) {
  diff_list <- list()
  
  # Calculate differences for each accession separately
  for (accession in accessions) {
    # Reshape data from wide to long format for current accession
    accession_data <- wide_data %>%
      select(IndexT, Group, Subtype, all_of(accession)) %>%
      pivot_wider(
        names_from = Group,
        values_from = all_of(accession)
      )
    
    # Calculate difference (HA - no) and filter valid cases
    accession_data <- accession_data %>%
      mutate(
        Difference = HA - no,  # Calculate difference between conditions
        Accession = accession   # Add accession identifier
      ) %>%
      # Remove cases with missing values in either condition
      filter(!is.na(Difference) & !is.na(no) & !is.na(HA))
    
    diff_list[[accession]] <- accession_data
  }
  
  # Combine results from all accessions
  all_differences <- bind_rows(diff_list)
  return(all_differences)
}

# Calculate differences for all selected accessions
diff_results <- calculate_differences(protein_wide, selected_accessions)

#  Step 1:Draw boxplots towards selcted mediators among subtypes -------------------------------------------------------
palette <- c("#cdf101ff", "#dc7320", "#6610f2")
for (accession in selected_accessions) {
  df_plot <- data.frame(
    Change = (diff_results %>% filter(Accession == accession) %>% filter(!Subtype %in% c("TC1", "TC2", "TC3")))$Difference,
    group = (diff_results %>% filter(Accession == accession) %>% filter(!Subtype %in% c("TC1", "TC2", "TC3")))$Subtype
  )
  
  p <- ggplot(df_plot, aes(x = group, y = Change,
                           fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, fill = NA, linewidth = 0.8,
                 position = position_dodge(width = 0.75)) +
    geom_jitter(shape = 16, size = 1.5, alpha = 0.8,
                position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = palette) +
    scale_colour_manual(values = palette) +
    labs(title = accession, x = NULL, y = "Change") +
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
  
  ggsave(filename = file.path("./picture/boxplot", paste0(accession, ".pdf")),
         plot = p, width = 5, height = 4, dpi = 300)
  rm(p)
}

#  Step 2:Wilcoxon rank-sum test to calculate significance -------------------------------------------------------------
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

# Iterate selected_accessions -----------------------------------------------------
plots <- list()
for (accession in selected_accessions) {
  plots[[accession]] <- Calculate_significance((diff_results %>% filter(Accession == accession) %>% filter(!Subtype %in% c("TC1", "TC2", "TC3")))$Difference,
                                               as.factor((diff_results %>% filter(Accession == accession) %>% filter(!Subtype %in% c("TC1", "TC2", "TC3")))$Subtype), 
                                               accession)
}

for (accession in names(plots)){
  ggsave(filename = file.path("./picture/heatmap", paste0(accession, ".pdf")),plot = plots[[accession]], width = 5, height = 4, dpi = 300)
}

# Compared the change in selected accessions of TC4 with other subtypes from baseline to HA treatment. -----------------------------------------------------
diff_results$is_TC4 <- ifelse(diff_results$Subtype == "TC4", "TC4", "Other")
unique_accessions <- unique(diff_results$Accession)
Wilcoxon_results_TC4andOthers <- data.frame(
  Accession = character(),
  p_value = numeric(),
  statistic = numeric(),
  stringsAsFactors = FALSE
)
for (accession in unique_accessions) {
  subset_data <- diff_results[diff_results$Accession == accession, ]
  test_result <- wilcox.test(no ~ is_TC4, data = subset_data)
  Wilcoxon_results_TC4andOthers <- rbind(Wilcoxon_results_TC4andOthers, 
                                         data.frame(
                                           Accession = accession,
                                           p_value = test_result$p.value,
                                           statistic = test_result$statistic
                                         ))
}
rm(test_result)


###########################################################################
#  Function to compare individuals in certain subtype before and after HA treatment using paired two-tailed Student t test   -----------------------------------------
###########################################################################
library(dplyr)
library(tidyr)

protein_long <- protein_wide %>% pivot_longer(cols = c("H0YDX6", "P02649", "J3KP74","Q92954"),names_to = "Protein",values_to = "Expression")

paired_Student_t_results <- protein_long %>%
  group_by(Subtype, Protein) %>%
  filter(n_distinct(IndexT) >= 3) %>%  
  do({
    no_data <- filter(., Group == "no") %>% arrange(IndexT)
    ha_data <- filter(., Group == "HA") %>% arrange(IndexT)
    
    common_patients <- intersect(no_data$IndexT, ha_data$IndexT)
    no_data <- filter(no_data, IndexT %in% common_patients)
    ha_data <- filter(ha_data, IndexT %in% common_patients)
    
    if(length(common_patients) >= 3) {
      t_test <- t.test(no_data$Expression, ha_data$Expression, 
                       paired = TRUE, alternative = "two.sided")
      
      data.frame(
        Patient_count = length(common_patients),
        t_statistic = t_test$statistic,
        p_value = t_test$p.value,
        mean_no = mean(no_data$Expression),
        mean_ha = mean(ha_data$Expression),
        mean_diff = mean(ha_data$Expression - no_data$Expression)
      )
    } else {
      data.frame(
        Patient_count = length(common_patients),
        t_statistic = NA,
        p_value = NA,
        mean_no = NA,
        mean_ha = NA,
        mean_diff = NA
      )
    }
  })

#  Change of Quantification value of selected protein before and after HA treatment in each subtype respectively   -----------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

no_color <- "#FF69B4"     
ha_color <- "#40e0d0"   

subtypes <- unique(protein_long$Subtype)

for(subtype in subtypes) {
  
  subtype_dir <- paste0("./picture/barplot/", subtype)
  if(!dir.exists(subtype_dir)) {
    dir.create(subtype_dir)
  }
  
  subtype_data <- protein_long %>% 
    filter(Subtype == subtype)
  
  proteins <- unique(subtype_data$Protein)
  
  for(protein in proteins) {
    
    protein_data <- subtype_data %>% 
      filter(Protein == protein)
    
    p <- ggplot(protein_data, aes(x = Group, y = Expression, group = IndexT)) +
      geom_line(color = "black", alpha = 0.6, linewidth = 0.8) +
      geom_point(aes(color = Group), size = 3, alpha = 0.8) +
      scale_fill_manual(values = c("no" = no_color, "HA" = ha_color)) +
      scale_color_manual(values = c("no" = no_color, "HA" = ha_color)) +
      scale_x_discrete(expand = expansion(mult = 0.2))+
      theme_minimal() +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(family = "Helvetica", size = 10, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title.x = element_text(family = "Helvetica", size = 12, face = "bold"),
        axis.title.y = element_text(family = "Helvetica", size = 12, face = "bold"),
        plot.title = element_text(family = "Helvetica", size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(family = "Helvetica", size = 12, face = "bold", hjust = 0.5),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line() 
      ) +
      scale_y_continuous(limits = c(0, NA))
    
    # filename <- paste0(subtype_dir, "/", subtype, "_", protein, "_expression.png")
    # ggsave(filename, p, width = 5.5, height = 8, dpi = 300)
    filename <- paste0(subtype_dir, "/", subtype, "_", protein, "_expression.pdf")
    ggsave(filename, p, width = 5.5, height = 8, dpi = 300)
    
    cat("Save to:", filename, "\n")
  }
}

