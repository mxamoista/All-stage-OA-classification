library(readxl)
library(dplyr)
library(ggplot2)
library(DescTools)

#Import information 
OA_drivers <- read_excel("./annotation/Information_update.xlsx")

# names(OA_drivers)[names(OA_drivers) == "LYM \r\npercentage"] <- "LYM percentage"
names(OA_drivers)[names(OA_drivers) == "NEU \r\ncount"] <- "NEU count"

#Add a line to record the neutrophil-to-lymphocyte ratio (NLR) to reflect the systemic inflammation status.
OA_drivers$NLR <- as.numeric(as.character(OA_drivers$`NEU count`)) / as.numeric(as.character(OA_drivers$`LYM count`))

#Add a line to record leukocytes count.
OA_drivers$`WBC count` <- as.numeric(as.character(OA_drivers$`NEU count`))+as.numeric(as.character(OA_drivers$`LYM count`))+
                                        as.numeric(as.character(OA_drivers$`MON count`))+as.numeric(as.character(OA_drivers$`EOS count`))+
                                        as.numeric(as.character(OA_drivers$`BAS percentage`))
#Add a line to record total bilirubin.
OA_drivers$TBIL <- as.numeric(as.character(OA_drivers$DBIL))+as.numeric(as.character(OA_drivers$IBIL))

#Add a line to record OA stage.
# OA_drivers <- OA_drivers %>% mutate(Stage = case_when(
#   KL %in% c(0, 1) ~ "mild",
#   KL %in% c(2, 3) ~ "middle",
#   KL == 4 ~ "severe",
#   TRUE ~ NA_character_  
# ))

#------------------------------------------------------------------------------------------------------#
#Draw a boxplot to compare BMI distributions across subtypes
BMI <- OA_drivers %>%
  filter(Subtype != "N" & BMI != "NA") %>%  
  mutate(BMI = as.numeric(as.character(BMI))) %>% 
  ggplot(aes(x = Subtype, y = BMI, fill = Subtype)) +
  geom_boxplot() +
  labs(
    title = " ",
    x = "Subtype",
    y = "BMI"
  ) +
  scale_fill_manual(values = c("TC1" = "#66c2a5", 
                               "TC2" = "#fc8d62",
                               "TC3" = "#8da0cb",
                               "TC4" = "#e78ac3",
                               "TC5" = "#e7c2a5",
                               "TC6" = "#a78ac4"
  )) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(size = 10,face = "bold"),
    axis.text.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_text(size = 10,face = "bold"),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.line = element_line(size = 0.5), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  
  )
pairwise.wilcox.test(
  x = (OA_drivers %>%
         filter(Subtype != "N" & BMI != "NA" & BMI != "0" ) %>%
         mutate(value = as.numeric(as.character(BMI))))$value,
  g = (OA_drivers %>%
         filter(Subtype != "N" & BMI != "NA" & BMI != "0" ) %>%
         mutate(value = as.numeric(as.character(BMI))))$Subtype,
  p.adjust.method = "BH"
)

#Classify obesity by BMI and report proportions and counts
OA_drivers <- OA_drivers %>% mutate(Obesity = ifelse(BMI >= 27.5, "Y", "N"))
Obesity_plot <- OA_drivers %>%
  filter(Subtype != "N" & Obesity != "NA") %>%
  group_by(Subtype, Obesity) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(label = sprintf("%.1f%%\n(n=%d)", percentage, count)) %>%
    mutate(gender = factor(Obesity, levels = c("Y", "N"))) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Obesity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Y" = "#f0768b",  
                               "N" = "#2dca86"), 
                    labels = c("No", "Yes"))  +
  theme_minimal() +
  theme(legend.position = "right",  
        axis.text = element_text(size = 10),  
        axis.title = element_text(size = 10),  
        legend.text = element_text(size = 8),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) 

# Create a contingency table
Obesity_data <- xtabs(~ Subtype + Obesity, data = OA_drivers %>%filter(Subtype != "N" & Obesity != "NA"))

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Obesity_data[c(1,2),])) 
CochranArmitageTest(t(Obesity_data[c(1,3),])) 
CochranArmitageTest(t(Obesity_data[c(1,4),])) 
CochranArmitageTest(t(Obesity_data[c(1,5),])) 
CochranArmitageTest(t(Obesity_data[c(1,6),])) 
CochranArmitageTest(t(Obesity_data[c(2,3),])) 
CochranArmitageTest(t(Obesity_data[c(2,4),]))#p-value = 0.0916
CochranArmitageTest(t(Obesity_data[c(2,5),])) 
CochranArmitageTest(t(Obesity_data[c(2,6),])) 
CochranArmitageTest(t(Obesity_data[c(3,4),])) 
CochranArmitageTest(t(Obesity_data[c(3,5),]))
CochranArmitageTest(t(Obesity_data[c(3,6),]))
CochranArmitageTest(t(Obesity_data[c(4,5),]))
CochranArmitageTest(t(Obesity_data[c(4,6),]))#p-value = 0.08779
CochranArmitageTest(t(Obesity_data[c(5,6),]))

#------------------------------------------------------------------------------------------------------#
#Draw a boxplot to compare GLU concentration distributions across subtypes
GLU <- OA_drivers %>%
  filter(Subtype != "N" & GLU != "NA") %>%  
  mutate(GLU = as.numeric(as.character(GLU))) %>% 
  ggplot(aes(x = Subtype, y = GLU, fill = Subtype)) +
  geom_boxplot() +
  labs(
    title = " ",
    x = "Subtype",
    y = "GLU"
  ) +
  scale_fill_manual(values = c("TC1" = "#66c2a5", 
                               "TC2" = "#fc8d62",
                               "TC3" = "#8da0cb",
                               "TC4" = "#e78ac3",
                               "TC5" = "#e7c2a5",
                               "TC6" = "#a78ac4"
  )) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(size = 10,face = "bold"),
    axis.text.y = element_text(size = 10,face = "bold"),
    axis.title.x = element_text(size = 10,face = "bold"),
    axis.title.y = element_text(size = 10,face = "bold"),
    axis.line = element_line(size = 0.5), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  
  )
pairwise.wilcox.test(
  x = (OA_drivers %>%
         filter(Subtype != "N" & GLU != "NA" & GLU != "0" ) %>%
         mutate(value = as.numeric(as.character(GLU))))$value,
  g = (OA_drivers %>%
         filter(Subtype != "N" & GLU != "NA" & GLU != "0" ) %>%
         mutate(value = as.numeric(as.character(GLU))))$Subtype,
  p.adjust.method = "BH"
)

#Classify diabetes by GLU and report proportions and counts
OA_drivers <- OA_drivers %>% mutate(Diabetes = ifelse(GLU >= 6.11, "Y", "N"))
Diabetes_plot <- OA_drivers %>%
  filter(Subtype != "N" & Diabetes != "NA") %>%
  group_by(Subtype, Diabetes) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(label = sprintf("%.1f%%\n(n=%d)", percentage, count)) %>%
  mutate(gender = factor(Diabetes, levels = c("Y", "N"))) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Diabetes)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Y" = "#f0768b",  
                               "N" = "#2dca86"), 
                    labels = c("No", "Yes"))  +
  theme_minimal() +
  theme(legend.position = "right",  
        axis.text = element_text(size = 10),  
        axis.title = element_text(size = 10),  
        legend.text = element_text(size = 8),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) 

Diabetes_data <- xtabs(~ Subtype + Diabetes, data = OA_drivers %>%filter(Subtype != "N" & Diabetes != "NA"))

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Diabetes_data[c(1,2),])) 
CochranArmitageTest(t(Diabetes_data[c(1,3),])) 
CochranArmitageTest(t(Diabetes_data[c(1,4),])) 
CochranArmitageTest(t(Diabetes_data[c(1,5),])) 
CochranArmitageTest(t(Diabetes_data[c(1,6),])) 
CochranArmitageTest(t(Diabetes_data[c(2,3),])) 
CochranArmitageTest(t(Diabetes_data[c(2,4),])) 
CochranArmitageTest(t(Diabetes_data[c(2,5),])) 
CochranArmitageTest(t(Diabetes_data[c(2,6),])) 
CochranArmitageTest(t(Diabetes_data[c(3,4),])) 
CochranArmitageTest(t(Diabetes_data[c(3,5),]))
CochranArmitageTest(t(Diabetes_data[c(3,6),]))
CochranArmitageTest(t(Diabetes_data[c(4,5),]))
CochranArmitageTest(t(Diabetes_data[c(4,6),]))
CochranArmitageTest(t(Diabetes_data[c(5,6),]))

#------------------------------------------------------------------------------------------------------#
#Create a Trauma Plot with Proportion and Count
Trauma_plot <- OA_drivers %>%
  filter(Subtype != "N" & Trauma != "NA") %>%
  group_by(Subtype, Trauma) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(label = sprintf("%.1f%%\n(n=%d)", percentage, count)) %>%
  mutate(gender = factor(Trauma, levels = c("Y", "N"))) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Trauma)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Y" = "#f0768b",  
                               "N" = "#2dca86"), 
                    labels = c("No", "Yes"))  +
  theme_minimal() +
  theme(legend.position = "right",  
        axis.text = element_text(size = 10),  
        axis.title = element_text(size = 10),  
        legend.text = element_text(size = 8),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) 

# Create a contingency table
Trauma_data <- xtabs(~ Subtype + Trauma, data = OA_drivers %>%filter(Subtype != "N" & Trauma != "NA"))

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Trauma_data[c(1,2),])) 
CochranArmitageTest(t(Trauma_data[c(1,3),])) 
CochranArmitageTest(t(Trauma_data[c(1,4),])) 
CochranArmitageTest(t(Trauma_data[c(1,5),])) #p-value = 0.06525
CochranArmitageTest(t(Trauma_data[c(1,6),])) 
CochranArmitageTest(t(Trauma_data[c(2,3),]))
CochranArmitageTest(t(Trauma_data[c(2,4),]))
CochranArmitageTest(t(Trauma_data[c(2,5),])) 
CochranArmitageTest(t(Trauma_data[c(2,6),])) 
CochranArmitageTest(t(Trauma_data[c(3,4),])) 
CochranArmitageTest(t(Trauma_data[c(3,5),]))
CochranArmitageTest(t(Trauma_data[c(3,6),]))
CochranArmitageTest(t(Trauma_data[c(4,5),]))
CochranArmitageTest(t(Trauma_data[c(4,6),]))
CochranArmitageTest(t(Trauma_data[c(5,6),])) #p-value = 0.09426

#------------------------------------------------------------------------------------------------------#
#Create a Menopause Plot with Proportion and Count
Menopause_plot <- OA_drivers %>%
  filter(Sex == "F" & Subtype != "N" & Menopause != "NA") %>%
  group_by(Subtype, Menopause) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(label = sprintf("%.1f%%\n(n=%d)", percentage, count)) %>%
  mutate(gender = factor(Menopause, levels = c("Y", "N"))) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Menopause)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Y" = "#f0768b",  
                               "N" = "#2dca86"), 
                    labels = c("No", "Yes"))  +
  theme_minimal() +
  theme(legend.position = "right",  
        axis.text = element_text(size = 10),  
        axis.title = element_text(size = 10),  
        legend.text = element_text(size = 8),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) 

# Create a contingency table
Menopause_data <- xtabs(~ Subtype + Menopause, data = OA_drivers %>%filter(Subtype != "N" & Menopause != "NA"))

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Menopause_data[c(1,2),])) 
CochranArmitageTest(t(Menopause_data[c(1,3),])) 
CochranArmitageTest(t(Menopause_data[c(1,4),])) 
CochranArmitageTest(t(Menopause_data[c(1,5),])) 
CochranArmitageTest(t(Menopause_data[c(1,6),])) 
CochranArmitageTest(t(Menopause_data[c(2,3),]))
CochranArmitageTest(t(Menopause_data[c(2,4),]))
CochranArmitageTest(t(Menopause_data[c(2,5),])) #p-value = 0.09673
CochranArmitageTest(t(Menopause_data[c(2,6),])) 
CochranArmitageTest(t(Menopause_data[c(3,4),])) 
CochranArmitageTest(t(Menopause_data[c(3,5),]))
CochranArmitageTest(t(Menopause_data[c(3,6),]))
CochranArmitageTest(t(Menopause_data[c(4,5),]))
CochranArmitageTest(t(Menopause_data[c(4,6),]))
CochranArmitageTest(t(Menopause_data[c(5,6),])) #p-value = 0.01723

#------------------------------------------------------------------------------------------------------#
# Visualize 14 immunoinflammatory markers (CRP, platelets, leukocytes, basophils, eosinophils, lymphocytes, monocytes, neutrophils; absolute counts and leukocyte percentages)  
# Compare distributions across subtypes with boxplots
immunoinflammatory_markers <- c(
  "platelet count", "CRP", "WBC count", "NLR", "BAS count", "BAS percentage",
  "EOS count", "EOS percentage", "LYM count", "LYM percentage",
  "MON count", "MON percentage", "NEU count", "NEU percentage"
)

subtype_colors <- c("TC1" = "#2c9061", "TC2" = "#147ab0","TC3" = "#cf0e41","TC4" = "#cdf101ff","TC5" = "#dc7320","TC6" = "#6610f2")
for (col_name in immunoinflammatory_markers) {
  p <- OA_drivers %>%
    filter(Subtype != "N" & !is.na(!!sym(col_name))) %>%  
    mutate(value = as.numeric(as.character(!!sym(col_name)))) %>% 
    ggplot(aes(x = Subtype, y = value, fill = Subtype)) +
    geom_boxplot() +
    labs(
      title = " ",
      x = "Subtype",
      y = col_name  
    ) +
    scale_fill_manual(values = subtype_colors) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      axis.text.x = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.line = element_line(size = 0.5), 
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  
    )
  file_name <- paste0('./picture/Non_heritable_factors/immunoinflammatory_markers/', gsub(" ", "_", col_name), ".pdf")  
  ggsave(
    filename = file_name,
    plot = p,
    device = "pdf",
    width = 336 / 72,  
    height = 231 / 72,
    units = "in",
    dpi = 300
  )
  message("Save to: ", file_name)
}

#------------------------------------------------------------------------------------------------------#
#Significance caculate
#inflammatory markers & liver function & glucose level & clinical scores of joint narrowing and osteophyte 
all_results <- list()

for (col in c(colnames(OA_drivers)[6:49], colnames(OA_drivers)[53:55])) {
  filtered_data <- OA_drivers %>%
    filter(Subtype != "N" & !is.na(.data[[col]])) %>%
    mutate(value = as.numeric(as.character(.data[[col]])))
  test_result <- pairwise.wilcox.test(
    x = filtered_data$value,
    g = filtered_data$Subtype,
    p.adjust.method = "BH"
  )
  result_df <- as.data.frame(test_result$p.value) %>%
    tibble::rownames_to_column("Group1") %>%
    tidyr::pivot_longer(
      cols = -Group1,
      names_to = "Group2",
      values_to = "p.value"
    ) %>%
    filter(!is.na(p.value)) %>%
    mutate(
      Variable = col,
      Method = "Mann–Whitney U test",
      p.adjust.method = "BH"
    ) %>%
    select(Variable, Group1, Group2, p.value, Method, p.adjust.method)
  all_results[[col]] <- result_df
}

final_results <- bind_rows(all_results)


