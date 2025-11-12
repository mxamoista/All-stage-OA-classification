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
Inflammation = c("TDO2","CD74", "IL-17A", "IL-17","IL-6","TNF","IL-40",
                 "MMP3","sVCAM1","sICAM1","TIMP1","VEGF","MCP1")
Cartilage_degeneration = c("MMP2", "RANKL", "SPARC", "APOE")
Pain = c("IL-6","IL-17A","COMP","PIICP")
# OA_symptoms = c("MMP3","sVCAM1","sICAM1","TIMP1","VEGF","MCP1")
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
#Select the proteins CD44, APOE, KIT, and PRG4, with corresponding UniProt IDs H0YDX6, P02649, P10721, and J3KP74, respectively.
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
      geom_line(color = "black", alpha = 0.6, linewidth = 0.5) +
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

