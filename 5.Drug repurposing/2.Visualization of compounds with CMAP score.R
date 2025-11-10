library(ggplot2)
library(ggrepel)

# Define the subtype names
subtypes <- paste0("TC", 1:6)

for (subtype in subtypes) {
  # Get the data for current subtype
  subtype_data <- get(paste0(subtype, "_drug_instances_all"))
  subtype_data <- subtype_data[subtype_data$cell_id == "PC3", ]
  subtype_data <- subtype_data[order(subtype_data$cmap_score), ]
  subtype_data$rank <- 1:nrow(subtype_data)
  subtype_data$cmap_score <- as.numeric(as.character(subtype_data$cmap_score))
  subtype_data$rank <- as.numeric(subtype_data$rank)
  
  # Merge with OA drugs
  label_compounds <- subtype_data[subtype_data$drug_name %in% OA_DRUG$DrugBank.ID, ]
  label_compounds <- merge(label_compounds, OA_DRUG, by.x = "drug_name", by.y = "DrugBank.ID")
  label_compounds <- label_compounds[order(label_compounds$rank), ]
  
  assign(paste0(subtype,"_label_compounds"),label_compounds)
  
  # Create the base plot
  p <- ggplot(subtype_data, aes(x = rank, y = cmap_score)) +
    geom_point(size = 2, color = "black") +
    geom_hline(yintercept = 0, color = "gray50") +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.background  = element_rect(fill = "white"),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black"),
      axis.ticks       = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      plot.margin      = margin(10, 10, 10, 10)
    ) +
    labs(x = "Compound rank", y = "CMap score",
         title = paste("CMap score for", subtype)) +
    scale_x_continuous(limits = c(0, max(subtype_data$rank)),
                       expand = c(0.02, 0)) +
    scale_y_continuous(limits = c(min(subtype_data$cmap_score), max(subtype_data$cmap_score)),
                       expand = c(0.02, 0))
  
  # Add highlights for negative scores
  if (nrow(subset(label_compounds, cmap_score < 0)) > 0) {
    p <- p +
      geom_point(data = subset(label_compounds, cmap_score < 0),
                 aes(x = rank, y = cmap_score), 
                 color = "red", size = 3) +
      geom_text_repel(data = subset(label_compounds, cmap_score < 0),
                      aes(x = rank, y = cmap_score, 
                          label = paste0(Drug, "\n", round(cmap_score, 4))),
                      color = "red", 
                      box.padding = unit(1, "lines"),  
                      point.padding = unit(2, "lines"),
                      segment.size = 0.4, 
                      min.segment.length = 0.2,
                      max.overlaps = Inf,
                      force = 5.0)
  }
  
  # Save as TIFF
  tiff_filename <- paste0("./picture/Graphs of compounds/", subtype, "_cmap_score.tiff")
  ggsave(tiff_filename, plot = p, width = 6, height = 6, dpi = 300, compression = "lzw")
  assign(paste0("p_", subtype), p)
  rm(p,subtype,label_compounds)
  message("Saved: ", tiff_filename)
}