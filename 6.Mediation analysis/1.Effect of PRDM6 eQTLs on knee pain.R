library(TwoSampleMR)
library(dplyr)
library("readxl")
#library(coloc)
#library(VariantAnnotation)
library(gwasglue)
library(plyr)
library(ggplot2)
library(MRPRESSO)

#add token to renviron_path
renviron_path <- file.path(Sys.getenv("HOME"), ".Renviron")
write(paste0("openGWAS_jwt=", "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxMjMxODAxN0B6anUuZWR1LmNuIiwiaWF0IjoxNzQ3MjI0NjQ2LCJleHAiOjE3NDg0MzQyNDZ9.GXmArVA8A0xz0uYz1Sqx50F3RBl7kf50mpxBbbBJqBtblAGwwTCepS5r0aJ8PsKwnjT77WEFNunGJvLrrDHqx_G3wnmcKpkosFJhvyVJ7wol3QFXlwRYeqSGSzRFagvlikTgFpz7viRMrNCAM_kt5THc5ZRnV3jXyu8Mt6lTn1ifb4L3_KM2W04E8of3X7pvUneKfAVaPgSmaxYg9XmQoHZYraUX9wVJzJIR0bpwx4XZVkSwUBgv-kHc2UiKfpNVplfPrm-IeGeqaGOMv3W_A9DgZzkSxLMal7i4axP3zMjwI5RzslKgvb4YTgFBK9oHZmAUPeEi_oogCBF9YoebZg"), 
      file = renviron_path, 
      append = TRUE)
readRenviron(renviron_path)
Sys.getenv("openGWAS_jwt")

eQTLs <- read.table("./data/PRDM6.Muscle_Skeletal.eqtl.txt", header = TRUE)
eQTLs <- eQTLs[eQTLs$p < 5e-8 & eQTLs$Freq > 0.01, ]

Ins <- format_data(eQTLs, type = "exposure", header = TRUE,
                 phenotype_col = "Probe", snp_col = "SNP", beta_col = "b",
                 se_col = "SE", eaf_col = "Freq", effect_allele_col = "A1",
                 other_allele_col = "A2", pval_col = "p")

#The extracted SNPs were used as instrument variables and clumped for independence at r2 > 0.8 with a window of 10,000 kb.
Ins_clumped <- clump_data(Ins, clump_kb=5000, clump_r2=0.8)

# GWAS data were used as the outcome
batch_extract <- function(snps, outcomes, batch_size = 500) {
  batches <- split(snps, ceiling(seq_along(snps)/batch_size))
  result <- lapply(batches, function(batch_snps) {
    extract_outcome_data(snps = batch_snps, outcomes = outcomes)
  })
  do.call(rbind, result)
}

outcome_dat <- batch_extract(
  snps = Ins_clumped$SNP,
  outcomes = c("ukb-b-6725", "ukb-b-16254", "ukb-b-476", "ukb-b-8906", "ukb-b-353")
)

#"ukb-b-16254", "ukb-b-8906"

##harmonise the exposure and outcome data
dat <- NULL
dat <- harmonise_data(
  exposure_dat = Ins_clumped, 
  outcome_dat = outcome_dat
)

##run the MR and sensitivity analyses 
mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL
try(mr_results <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw","mr_ivw_mre","mr_ivw_fe",
                                        "mr_simple_median","mr_weighted_median")))  # main MR analysis
mr_hetero <- mr_heterogeneity(dat) # heterogeneity test across instruments
mr_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test  
try(mr_single <- mr_singlesnp(dat)) #single SNP MR using Wald ratio
mr_scatter_plot1(mr_results, dat) #Scatterplot shows SnP effects on outcomes

mr_loo <- mr_leaveoneout(dat) #Leave-one-SNP-out analysis was performed to assess if the overall effect was driven by any single SNP 
mr_leaveoneout_plot(mr_loo)
result_or <- generate_odds_ratios(mr_results)
mr_presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
          SignifThreshold = 0.05)

################################## Plot MR
library(forestplot)
library(dplyr)
library(ggthemes)

#result_or1=result_or
#result_or1$or=log(result_or1$or)
#result_or1$or_lci95=log(result_or1$or_lci95)
#result_or1$or_uci95=log(result_or1$or_uci95)
#result_or1$pval=p.adjust(result_or1$pval,method='bonferroni',n=nrow(result_or1))

#For chronic knee pain 
result_or_chronic_pain <- result_or[result_or$id.outcome == 'ukb-b-8906' & result_or$method != "Inverse variance weighted",]

result_OR_chronic_pain <- tibble(
  OR  = log(result_or_chronic_pain$or),
  method = result_or_chronic_pain$method,
  p_value = as.character(result_or_chronic_pain$pval),
  lower = log(result_or_chronic_pain$or_lci95),
  upper = log(result_or_chronic_pain$or_uci95)
)

header <- tibble(
  method = "Method",
  p_value = "P value",
)

combined_data <- bind_rows(header, result_OR_chronic_pain)

colors <- c("red", "blue", "green", "purple", "orange")

result_OR_chronic_pain |>
  forestplot(labeltext = method,
             mean = OR,
             lower = lower,
             upper = upper,
             zero = 0,
             boxsize=0.25,
             col = fpColors(box = c("red"),line = "black",summary = "grey"),
             clip = c(0.0, 0.08), 
             xticks = c(0.0,0.04,0.08),
             lwd.ci = 2.0,
               txt_gp = fpTxtGp(
                 ticks = gpar(cex = 0.8, fontface = "bold"),
                 xlab = gpar(cex = 1, fontface = "bold")  
               ),
             xlab = "log OR") |>
  fp_add_header("Method") |>
  fp_set_style(lines = gpar(col = "darkblue"))

#For acute knee pain 
result_or_acute_pain <- result_or[result_or$id.outcome == 'ukb-b-16254' & result_or$method != "Inverse variance weighted",]

result_OR_acute_pain <- tibble(
  OR  = log(result_or_acute_pain$or),
  method = result_or_acute_pain$method,
  p_value = as.character(result_or_acute_pain$pval),
  lower = log(result_or_acute_pain$or_lci95),
  upper = log(result_or_acute_pain$or_uci95)
)

header <- tibble(
  method = "Method",
  p_value = "P value",
)

combined_data <- bind_rows(header, result_OR_acute_pain)

colors <- c("red", "blue", "green", "purple", "orange")

result_OR_acute_pain |>
  forestplot(labeltext = method,
             mean = OR,
             lower = lower,
             upper = upper,
             zero = 0,
             boxsize=0.25,
             col = fpColors(box = c("red"),line = "black",summary = "grey"),
             clip = c(0.0, 0.08), 
             xticks = c(0.0,0.04,0.08),
             lwd.ci = 2.0,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 0.8, fontface = "bold"),
               xlab = gpar(cex = 1, fontface = "bold")  
             ),
             xlab = "log OR") |>
  fp_add_header("Method") |>
  fp_set_style(lines = gpar(col = "darkblue"))
