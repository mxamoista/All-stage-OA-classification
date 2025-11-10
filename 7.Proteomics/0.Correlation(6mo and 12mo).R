# Pearson correlation coefficients (R) were calculated to assess the linear association 
# between the proteomics of synovial fluid at 6mo and 1year after HA treatment 
library(ggplot2)
correlation <- cor(x = Protein_quantification_filtered$`15`,y = Protein_quantification_filtered$`25`, method = "pearson")

#remove the largest and smallest value of protein quantification value
library(genefilter)
Protein_quantification_correlation <- data.frame(row.names = Protein_quantification_filtered$Accession,
                                                 sixmonths = Protein_quantification_filtered$`15`,
                                                 oneyear = Protein_quantification_filtered$`25`)
Protein_quantification_correlation_filtered <- Protein_quantification_correlation %>%
  mutate(across(everything(), ~ {
    lower_cutoff <- quantile(., 0.01, na.rm = TRUE)
    upper_cutoff <- quantile(., 0.99, na.rm = TRUE)
    ifelse(. <= lower_cutoff | . >= upper_cutoff, NA, .)
  })) %>%
  na.omit()

# Perform correlation analysis
correlation <- cor.test(x = Protein_quantification_correlation_filtered$sixmonths,
                        y = Protein_quantification_correlation_filtered$oneyear, 
                        method = "pearson")

# Extract correlation coefficient and p-value
correlation_value <- round(correlation$estimate, 4)
p_value <- format.pval(correlation$p.value, digits = 4)

ggplot(Protein_quantification_correlation_filtered, aes(x = sixmonths, y = oneyear)) +
  geom_point(shape = 16, color = "#87a8fa", size = 3) +  # Increase point size
  geom_smooth(method = "glm", color = "gray", se = TRUE) +  # Solid black regression line
  labs(x = "6 mo", y = "12 mo") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title = element_text(face = "bold"),  
    axis.line = element_line(color = "black", linewidth = 1.2),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    text = element_text(size = 18, face = "bold")
  ) +
  annotate("text", x = 10, y = max(Protein_quantification_correlation_filtered), 
           label = paste("R =", correlation_value, "\nP-value ", p_value), 
           hjust = 0, size = 3, fontface = "bold", vjust = 1)
ggsave("./picture/correlation/6mo_1year.pdf", limitsize = F, width = 8, height = 5)


