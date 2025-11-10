library(dplyr)
library(ggplot2)
library(DescTools)
library(readxl)
library(tidyr)
library(tibble) 

# 1. Statistics analysis of clinical information including age,gender,side and KL score  -------------------------------------------------
#Import clinical scores 
Information <- read_excel("./annotation/Information_update.xlsx")
Basic_Information =  data.frame(Index = Information$Index, Age = as.factor(Information$Age), Sex = Information$Sex , Side = Information$Side,KL = Information$KL)

#Overlap with the result of Classification result
annotation_col <- read.table("./annotation/annotation_col.txt",header = TRUE , row.names =1,sep = "\t",comment.char = "")
Basic_Information <- Basic_Information[Basic_Information$Index %in% rownames(annotation_col), ]
Basic_Information <- Basic_Information[match(rownames(annotation_col),Basic_Information$Index),]
identical(Basic_Information$Index,rownames(annotation_col)) #[1] TRUE
Basic_Information$Subtype <- annotation_col$Subtype

#Age at diagnosis
Age <- Basic_Information %>%
  filter(Subtype != "N" & Age != "NA" & Age != "0") %>%  
  mutate(Age = as.numeric(as.character(Age))) %>% 
  ggplot(aes(x = Subtype, y = Age, fill = Subtype)) +
  geom_boxplot() +
  labs(
    title = " ",
    x = "Subtype",
    y = "Age at diagnosis"
  ) +
  scale_fill_manual(values = c("TC1" = "#2c9061", "TC2" = "#147ab0","TC3" = "#cf0e41","TC4" = "#cdf101ff","TC5" = "#dc7320","TC6" = "#6610f2")) +
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

# Wilcoxon rank-sum test to calculate significance for age among subtypes
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

Significance_for_age <- Calculate_significance((Basic_Information %>%
                                                  filter(Subtype != "N" & Age != "NA" & Age != "0" ) %>%
                                                  mutate(value = as.numeric(as.character(Age))))$value,
                                               as.factor((Basic_Information %>%
                                                  filter(Subtype != "N" & Age != "NA" & Age != "0" ) %>%
                                                  mutate(value = as.numeric(as.character(Age))))$Subtype),
                                               "Age at diagnosis")

#KL score
KL_plot <- Basic_Information %>%
  filter(Subtype != "N" & KL != "N" & KL != "Null") %>%
  mutate(KL = factor(KL, levels = c("0", "1", "2", "3", "4"))) %>%
  group_by(Subtype, KL) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = KL)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0"= "#f6f8fb",
                               "1" = "#dbe4f0",
                               "2" = "#c9d7e9", 
                               "3" = "#b7cae2", 
                               "4" = "#a6bddb")) +
  labs(title = "KL", x = "", y = "Percentage %", fill = "KL score") +
  theme_minimal()


# Create a contingency table
KL_data <- xtabs(~ Subtype + KL, data = Basic_Information %>% filter(Subtype != "N" & KL != "N" & KL != "Null"))

# Perform the Cochran-Armitage test
CochranArmitageTest(t(KL_data[c(1,2),])) 
CochranArmitageTest(t(KL_data[c(1,3),])) 
CochranArmitageTest(t(KL_data[c(1,4),])) 
CochranArmitageTest(t(KL_data[c(1,5),])) 
CochranArmitageTest(t(KL_data[c(1,6),])) 
CochranArmitageTest(t(KL_data[c(2,3),])) 
CochranArmitageTest(t(KL_data[c(2,4),])) 
CochranArmitageTest(t(KL_data[c(2,5),])) 
CochranArmitageTest(t(KL_data[c(2,6),])) 
CochranArmitageTest(t(KL_data[c(3,4),])) 
CochranArmitageTest(t(KL_data[c(3,5),]))
CochranArmitageTest(t(KL_data[c(3,6),]))
CochranArmitageTest(t(KL_data[c(4,5),]))
CochranArmitageTest(t(KL_data[c(4,6),]))
CochranArmitageTest(t(KL_data[c(5,6),]))


#OA stage,the criteria for OA stage classification based on KL grade originated from paper below.
# "Metabolite asymmetric dimethylarginine (ADMA) functions as a destabilization enhancer of SOX9 mediated by DDAH1 in osteoarthritis"
OA_stage <- Basic_Information %>% mutate(
  Stage = case_when(
    KL %in% c(0, 1) ~ "mild",
    KL %in% c(2, 3) ~ "middle",
    KL == 4 ~ "severe",
    TRUE ~ NA_character_  
  )
)
write.table(OA_stage,"./annotation/OA_stage.txt",sep = '\t')

OA_stage_plot <- OA_stage %>%
  filter(Stage != "N" & Stage != "Null") %>%
  group_by(Subtype, Stage) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Stage)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("mild" = "#099396", "middle" = "#4642B1", "severe" = "#EBD7A5"),
    breaks = c("mild", "middle", "severe"),  
    labels = c("Mild", "Moderate", "Severe")  
  ) +
  labs(title = "OA stage", x = "", y = "Percentage %", fill = "OA stage") +
  theme_minimal()

# Create a contingency table
OA_stage_data <- xtabs(~ Subtype + Stage, data = OA_stage %>% filter(Stage != "N" & KL != "N" & KL != "Null"))
CochranArmitageTest(t(OA_stage_data[c(1,2),])) 
CochranArmitageTest(t(OA_stage_data[c(1,3),])) 
CochranArmitageTest(t(OA_stage_data[c(1,4),])) 
CochranArmitageTest(t(OA_stage_data[c(1,5),])) 
CochranArmitageTest(t(OA_stage_data[c(1,6),])) 
CochranArmitageTest(t(OA_stage_data[c(2,3),])) 
CochranArmitageTest(t(OA_stage_data[c(2,4),])) #p-value = 0.06427
CochranArmitageTest(t(OA_stage_data[c(2,5),])) #p-value = 0.06699
CochranArmitageTest(t(OA_stage_data[c(2,6),])) 
CochranArmitageTest(t(OA_stage_data[c(3,4),])) 
CochranArmitageTest(t(OA_stage_data[c(3,5),]))
CochranArmitageTest(t(OA_stage_data[c(3,6),]))
CochranArmitageTest(t(OA_stage_data[c(4,5),]))
CochranArmitageTest(t(OA_stage_data[c(4,6),]))
CochranArmitageTest(t(OA_stage_data[c(5,6),]))

#gender
gender_plot <- Basic_Information %>%
  filter(Sex != "N") %>%
  group_by(Subtype, Sex) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Sex)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("M"= "#2dca86",
                               "F" = "#f0768b")) +
  labs(title = "Sex", x = "", y = "Percentage %", fill = "Sex") +
  theme_minimal()

print(gender_plot)
# Create a contingency table
Sex_data <- xtabs(~ Subtype + Sex, data = Basic_Information %>% filter(Subtype != "N" & Sex != "N" ))

# Perform the Chi-square test to check association between the proportion of gender among subtypes
chisq.test(t(Sex_data[c(1,2),])) 
chisq.test(t(Sex_data[c(1,3),])) 
chisq.test(t(Sex_data[c(1,4),])) 
chisq.test(t(Sex_data[c(1,5),])) 
chisq.test(t(Sex_data[c(1,6),])) 
chisq.test(t(Sex_data[c(2,3),])) 
chisq.test(t(Sex_data[c(2,4),])) 
chisq.test(t(Sex_data[c(2,5),])) 
chisq.test(t(Sex_data[c(2,6),])) 
chisq.test(t(Sex_data[c(3,4),])) 
chisq.test(t(Sex_data[c(3,5),]))
chisq.test(t(Sex_data[c(3,6),]))
chisq.test(t(Sex_data[c(4,5),]))
chisq.test(t(Sex_data[c(4,6),]))
chisq.test(t(Sex_data[c(5,6),]))


#Side
Side_plot <- Basic_Information %>%
  filter(Side != "N") %>%
  group_by(Subtype, Side) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Side)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("R"= "#2dca86",
                               "L" = "#f0768b")) +
  labs(title = "Side", x = "", y = "Percentage %", fill = "Side") +
  theme_minimal()
print(Side_plot)

# Create a contingency table
Side_data <- xtabs(~ Subtype + Side, data = Basic_Information %>% filter(Subtype != "N" & Side != "N" ))

# Perform the Cochran-Armitage test
chisq.test(t(Side_data[c(1,2),])) 
chisq.test(t(Side_data[c(1,3),])) 
chisq.test(t(Side_data[c(1,4),])) 
chisq.test(t(Side_data[c(1,5),])) 
chisq.test(t(Side_data[c(1,6),])) 
chisq.test(t(Side_data[c(2,3),])) 
chisq.test(t(Side_data[c(2,4),])) 
chisq.test(t(Side_data[c(2,5),])) 
chisq.test(t(Side_data[c(2,6),])) 
chisq.test(t(Side_data[c(3,4),]))
chisq.test(t(Side_data[c(3,5),]))
chisq.test(t(Side_data[c(3,6),]))
chisq.test(t(Side_data[c(4,5),]))
chisq.test(t(Side_data[c(4,6),]))
chisq.test(t(Side_data[c(5,6),]))

# 2. Statistical analysis was conducted on imaging phenotypes, 
# such as Cartilage loss, bone marrow lesions,osteophytes and effusion-synovitis
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(DescTools)
MOAKS <- read_excel("./MOAKS/Radiology-MRI.xlsx")
MOAKS <- MOAKS  %>%
  filter(Index2 != "NA") %>%
  filter(Subtype != "#N/A")

##Cartilage loss
Cartilage_loss1_plot <- MOAKS %>%
  mutate(`Cartilage loss` = factor(`Cartilage loss`,levels = c("0","1","2","3","4","5","6","12"))) %>%
  group_by(Subtype, `Cartilage loss`) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = `Cartilage loss`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0" = "#cff0ec",
                               "1" = "#d9e8e6",
                               "2" = "#e2e1e1",
                               "3" = "#fdf1f3",
                               "4" = "#fceaed",
                               "5" = "#f1d1d6",
                               "6" = "#f8c8d0",
                               "12" = "#f296a6")) + 
  labs(title = "Cartilage loss", x = "", y = "Percentage %", fill = "Cartilage loss") +
  theme_minimal()

print(Cartilage_loss1_plot)
# Create a contingency table
Cartilage_loss1_data <- xtabs(~ Subtype + `Cartilage loss`, data = MOAKS)

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Cartilage_loss1_data [c(1,2),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(1,3),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(1,4),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(1,5),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(1,6),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(2,3),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(2,4),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(2,5),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(2,6),])) 
CochranArmitageTest(t(Cartilage_loss1_data [c(3,4),])) # p-value = 0.09331
CochranArmitageTest(t(Cartilage_loss1_data [c(3,5),]))
CochranArmitageTest(t(Cartilage_loss1_data [c(3,6),]))
CochranArmitageTest(t(Cartilage_loss1_data [c(4,5),]))
CochranArmitageTest(t(Cartilage_loss1_data [c(4,6),]))
CochranArmitageTest(t(Cartilage_loss1_data [c(5,6),]))

##Cartilage loss(full thickness)
Cartilage_loss2_plot <- MOAKS %>%
  mutate(`Cartilage loss (full thickness)` = as.factor(`Cartilage loss (full thickness)`)) %>%
  group_by(Subtype, `Cartilage loss (full thickness)`) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = `Cartilage loss (full thickness)`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0" = "#cff0ec",
                               "1" = "#d9e8e6",
                               "2" = "#e2e1e1",
                               "3" = "#ead9db",
                               "4" = "#f1d1d6",
                               "5" = "#f8c8d0",
                               "12" = "#eabac2")) + 
  labs(title = "Cartilage loss(full thickness)", x = "", y = "Percentage %", fill = "Cartilage loss(full thickness)") +
  theme_minimal()

print(Cartilage_loss2_plot)
# Create a contingency table
Cartilage_loss2_data <- xtabs(~ Subtype + `Cartilage loss (full thickness)`, data = MOAKS)

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Cartilage_loss2_data [c(1,2),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(1,3),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(1,4),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(1,5),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(1,6),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(2,3),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(2,4),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(2,5),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(2,6),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(3,4),])) 
CochranArmitageTest(t(Cartilage_loss2_data [c(3,5),]))
CochranArmitageTest(t(Cartilage_loss2_data [c(3,6),]))
CochranArmitageTest(t(Cartilage_loss2_data [c(4,5),]))
CochranArmitageTest(t(Cartilage_loss2_data [c(4,6),]))
CochranArmitageTest(t(Cartilage_loss2_data [c(5,6),]))

##BML
BML_plot <- MOAKS %>%
  mutate(`BML` = as.factor(`BML`)) %>%
  group_by(Subtype, `BML`) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = `BML`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0" = "#cff0ec",
                               "1" = "#d9e8e6",
                               "2" = "#e2e1e1",
                               "3" = "#ead9db",
                               "6" = "#f1d1d6",
                               "11" = "#f8c8d0")) + 
  labs(title = "BML", x = "", y = "Percentage %", fill = "BML") +
  theme_minimal()

print(BML_plot)
# Create a contingency table
BML_data <- xtabs(~ Subtype + `BML`, data = MOAKS)

# Perform the Cochran-Armitage test
CochranArmitageTest(t(BML_data [c(1,2),])) 
CochranArmitageTest(t(BML_data [c(1,3),])) 
CochranArmitageTest(t(BML_data [c(1,4),])) 
CochranArmitageTest(t(BML_data [c(1,5),])) 
CochranArmitageTest(t(BML_data [c(1,6),])) #p-value = 0.0542
CochranArmitageTest(t(BML_data [c(2,3),])) 
CochranArmitageTest(t(BML_data [c(2,4),])) 
CochranArmitageTest(t(BML_data [c(2,5),])) 
CochranArmitageTest(t(BML_data [c(2,6),])) 
CochranArmitageTest(t(BML_data [c(3,4),])) 
CochranArmitageTest(t(BML_data [c(3,5),]))
CochranArmitageTest(t(BML_data [c(3,6),]))
CochranArmitageTest(t(BML_data [c(4,5),])) 
CochranArmitageTest(t(BML_data [c(4,6),]))
CochranArmitageTest(t(BML_data [c(5,6),]))

##Osteophytes in X-ray
Information$Osteophyte <- 
  as.numeric(as.character(Information$`Medial femoral osteophyte`))+as.numeric(as.character(Information$`Medial tibial osteophyte`))+
  as.numeric(as.character(Information$`Lateral femoral osteophyte`))+as.numeric(as.character(Information$`Lateral tibial osteophyte`))

Osteophyte_plot <- Information %>%
  filter(Subtype != "N" & Osteophyte != "N" & Osteophyte != "Null") %>%
  mutate(Osteophyte = factor(Osteophyte, levels = c("2","3","4","5","6","7","8","9","10","11","12"))) %>%
  group_by(Subtype, Osteophyte) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = Osteophyte)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("2" = "#cff0ec",
                               "3" = "#d9e8e6",
                               "4" = "#e2e1e1",
                               "5" = "#ead9db",
                               "6" = "#f1d1d6",
                               "7" = "#f8c8d0",
                               "8" = "#ebbdc5",
                               "9" = "#e5aab3",
                               "10" = "#df96a2",
                               "11" = "#CD5C6F",
                               "12" = "#AA3448")) +
  labs(title = "Osteophyte", x = "", y = "Percentage %", fill = "Osteophyte score") +
  theme_minimal()


# Create a contingency table
Osteophyte_data <- xtabs(~ Subtype + Osteophyte, data = Information %>% filter(Subtype != "N" & Osteophyte != "N" & Osteophyte != "Null"))

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Osteophyte_data[c(1,2),])) 
CochranArmitageTest(t(Osteophyte_data[c(1,3),])) #p-value = 0.03633
CochranArmitageTest(t(Osteophyte_data[c(1,4),])) 
CochranArmitageTest(t(Osteophyte_data[c(1,5),])) 
CochranArmitageTest(t(Osteophyte_data[c(1,6),])) 
CochranArmitageTest(t(Osteophyte_data[c(2,3),])) 
CochranArmitageTest(t(Osteophyte_data[c(2,4),])) #p-value = 0.003179
CochranArmitageTest(t(Osteophyte_data[c(2,5),])) 
CochranArmitageTest(t(Osteophyte_data[c(2,6),])) 
CochranArmitageTest(t(Osteophyte_data[c(3,4),])) #p-value = 0.001056
CochranArmitageTest(t(Osteophyte_data[c(3,5),]))
CochranArmitageTest(t(Osteophyte_data[c(3,6),]))
CochranArmitageTest(t(Osteophyte_data[c(4,5),]))
CochranArmitageTest(t(Osteophyte_data[c(4,6),])) #p-value = 0.002704
CochranArmitageTest(t(Osteophyte_data[c(5,6),]))

##Osteophytes in MRI
Osteophytes_plot2 <- MOAKS %>%
  mutate(Osteophytes = factor(Osteophytes, levels = c("0","1","2","3","4","5","6","7","10"))) %>%
  group_by(Subtype, `Osteophytes`) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = `Osteophytes`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0" = "#cff0ec",
                               "1" = "#d9e8e6",
                               "2" = "#e2e1e1",
                               "3" = "#ead9db",
                               "4" = "#f1d1d6",
                               "5" = "#f8c8d0",
                               "6" = "#ebbdc5",
                               "7" = "#e5aab3",
                               "10" = "#df96a2")) + 
  labs(title = "Osteophytes", x = "", y = "Percentage %", fill = "Osteophytes") +
  theme_minimal()

print(Osteophytes_plot2)
# Create a contingency table
Osteophytes_data <- xtabs(~ Subtype + `Osteophytes`, data = MOAKS)

# Perform the Cochran-Armitage test
CochranArmitageTest(t(Osteophytes_data [c(1,2),])) #p-value = 0.09999
CochranArmitageTest(t(Osteophytes_data [c(1,3),])) 
CochranArmitageTest(t(Osteophytes_data [c(1,4),])) 
CochranArmitageTest(t(Osteophytes_data [c(1,5),])) 
CochranArmitageTest(t(Osteophytes_data [c(1,6),])) #p-value = 0.07481
CochranArmitageTest(t(Osteophytes_data [c(2,3),])) 
CochranArmitageTest(t(Osteophytes_data [c(2,4),])) 
CochranArmitageTest(t(Osteophytes_data [c(2,5),])) #p-value = 0.09287
CochranArmitageTest(t(Osteophytes_data [c(2,6),])) 
CochranArmitageTest(t(Osteophytes_data [c(3,4),])) 
CochranArmitageTest(t(Osteophytes_data [c(3,5),]))
CochranArmitageTest(t(Osteophytes_data [c(3,6),]))
CochranArmitageTest(t(Osteophytes_data [c(4,5),])) 
CochranArmitageTest(t(Osteophytes_data [c(4,6),]))
CochranArmitageTest(t(Osteophytes_data [c(5,6),])) #p-value = 0.03684


##effusion-synovitis
effusion_synovitis_plot <- MOAKS %>%
  mutate(`effusion-synovitis` = as.factor(`effusion-synovitis`)) %>%
  group_by(Subtype, `effusion-synovitis`) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Subtype) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ggplot(aes(x = Subtype, y = percentage, fill = `effusion-synovitis`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("0" = "#cff0ec",
                               "1" = "#d9e8e6",
                               "2" = "#e2e1e1",
                               "3" = "#ead9db")) + 
  labs(title = "effusion-synovitis", x = "", y = "Percentage %", fill = "effusion-synovitis") +
  theme_minimal()

print(effusion_synovitis_plot)
# Create a contingency table
effusion_synovitis_data <- xtabs(~ Subtype + `effusion-synovitis`, data = MOAKS)

# Perform the Cochran-Armitage test
CochranArmitageTest(t(effusion_synovitis_data [c(1,2),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(1,3),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(1,4),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(1,5),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(1,6),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(2,3),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(2,4),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(2,5),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(2,6),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(3,4),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(3,5),])) 
CochranArmitageTest(t(effusion_synovitis_data [c(3,6),]))
CochranArmitageTest(t(effusion_synovitis_data [c(4,5),]))
CochranArmitageTest(t(effusion_synovitis_data [c(4,6),]))
CochranArmitageTest(t(effusion_synovitis_data [c(5,6),]))

