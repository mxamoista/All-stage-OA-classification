########################
library(data.table)
library(readr)
library(readxl)
library(TwoSampleMR)
#Perform two-step MR
################################### Process1,the first step was to estimate the causal effect of the PRDM6 on the mediator using UVMR.
mediator <- extract_outcome_data(
  snps=dat$SNP,
  outcomes=c("ukb-b-8147","ukb-a-194","finn-b-RX_PARACETAMOL_NSAID",
             "ukb-b-16288","ukb-b-12506","ukb-a-463","ukb-b-6888","ukb-b-14518",
             "ukb-b-6549","ukb-b-4991","ukb-b-12648",
             "ukb-b-9386",
             "ukb-b-6911","ukb-a-156",
             "ukb-b-14037","ukb-b-11535","ukb-a-188","ukb-a-494",
             "ukb-b-16452"
             ),
  proxies = TRUE,
  maf_threshold = 0.01
)

exposure_to_mediator <- harmonise_data(
  exposure_dat = Ins_clumped,
  outcome_dat = mediator,
  action= 2
)

res_exposure_to_mediator <- mr(exposure_to_mediator)
res_exposure_to_mediator_or <- generate_odds_ratios(res_exposure_to_mediator)


################################### Process2,the second step was to estimate the causal effect of each mediator on the knee pain. 
mediator_as_exposure<- format_data(mediator, type = "exposure", header = TRUE, snp_col = "SNP", beta_col = "beta.outcome",
                   phenotype_col = "outcome",se_col = "se.outcome", eaf_col = "eaf.outcome", effect_allele_col = "effect_allele.outcome",
                   other_allele_col = "other_allele.outcome", pval_col = "pval.outcome")

mediator_to_outcome <- harmonise_data(
  exposure_dat = mediator_as_exposure,
  outcome_dat = outcome_dat,
  action= 2
)
res_mediator_to_outcome <- mr(mediator_to_outcome)
res_mediator_to_outcome_or <- generate_odds_ratios(res_mediator_to_outcome)


################################## Plot MR
library(forestplot)
library(dplyr)
# mapping_df <- mediator %>%
#   distinct(originalname.outcome, outcome) %>%
#   rename(exposure = originalname.outcome, 
#          mapped_outcome = outcome)
# 
# res_mediator_to_outcome_or2 <- res_mediator_to_outcome_or %>%
#   mutate(exposure_name = mapping_df$mapped_outcome[match(exposure, mapping_df$exposure)]) %>%
#   relocate(exposure, .before = 1)  

#Select the common mediator from process1 and process2
intersect(res_exposure_to_mediator_or[
  res_exposure_to_mediator_or$pval < 0.05 & 
    res_exposure_to_mediator_or$method == "Inverse variance weighted", 
]$outcome,
res_mediator_to_outcome_or[
  res_mediator_to_outcome_or$pval < 0.05 & 
    res_mediator_to_outcome_or$method == "Inverse variance weighted", 
]$exposure)

#[1] "Paracetamol of NSAID medication || id:finn-b-RX_PARACETAMOL_NSAID" [2] "Vitamin and mineral supplements: Vitamin D || id:ukb-b-12648"

# Draw the forest plot showing the effect of PRDM6 expression on each mediator.
res_exposure_to_mediator_or_filter <- res_exposure_to_mediator_or[
  res_exposure_to_mediator_or$pval < 0.05 & 
    res_exposure_to_mediator_or$method == "Inverse variance weighted" & 
    (res_exposure_to_mediator_or$outcome == "Paracetamol of NSAID medication || id:finn-b-RX_PARACETAMOL_NSAID" | 
       res_exposure_to_mediator_or$outcome == "Vitamin and mineral supplements: Vitamin D || id:ukb-b-12648"), 
]

result_OR_exposure_to_mediator <- tibble(
  outcome = as.character(res_exposure_to_mediator_or_filter$id.outcome),
  OR  = log(res_exposure_to_mediator_or_filter$or),
  p_value = as.character(res_exposure_to_mediator_or_filter$pval),
  lower = log(res_exposure_to_mediator_or_filter$or_lci95),
  upper = log(res_exposure_to_mediator_or_filter$or_uci95)
)

result_OR_exposure_to_mediator |>
  forestplot(labeltext = outcome,
             mean = OR,
             lower = lower,
             upper = upper,
             zero = 0,
             boxsize=0.05,
             col = fpColors(box = c("red"),line = "black",summary = "grey"),
             clip = c(-0.002, 0.04), 
             xticks = c(-0.002,0.0,0.02,0.04),
             lwd.ci = 2.0,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 0.8, fontface = "bold"),
               xlab = gpar(cex = 1, fontface = "bold")  
             ),
             xlab = "log OR") |>
  fp_add_header("outcome") |>
  fp_set_style(lines = gpar(col = "darkblue"))
    
# Draw the forest plot showing the effect of mediator on chronic knee pain.
res_mediator_to_outcome_or_filter <- res_mediator_to_outcome_or[
  res_mediator_to_outcome_or$pval < 0.05 & 
    res_mediator_to_outcome_or$method == "Inverse variance weighted" & 
    (res_mediator_to_outcome_or$exposure == "Paracetamol of NSAID medication || id:finn-b-RX_PARACETAMOL_NSAID" | 
       res_mediator_to_outcome_or$exposure == "Vitamin and mineral supplements: Vitamin D || id:ukb-b-12648"), 
]

result_OR_mediator_to_outcome <- tibble(
  mediator = as.character(res_mediator_to_outcome_or_filter$exposure),
  OR  = log(res_mediator_to_outcome_or_filter$or),
  p_value = as.character(res_mediator_to_outcome_or_filter$pval),
  lower = log(res_mediator_to_outcome_or_filter$or_lci95),
  upper = log(res_mediator_to_outcome_or_filter$or_uci95)
)

result_OR_mediator_to_outcome <- result_OR_mediator_to_outcome %>% 
  slice(c(2, 1))

result_OR_mediator_to_outcome |>
  forestplot(labeltext = mediator,
             mean = OR,
             lower = lower,
             upper = upper,
             zero = 0,
             boxsize=0.07,
             col = fpColors(box = c("red"),line = "black",summary = "grey"),
             clip = c(-12.5, 0.5), 
             xticks = c(-12.5,-6.25,0.0,0.25,0.5),
             lwd.ci = 2.0,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 0.8, fontface = "bold"),
               xlab = gpar(cex = 1, fontface = "bold")  
             ),
             xlab = "log OR") |>
  fp_add_header("outcome") |>
  fp_set_style(lines = gpar(col = "darkblue"))


#Caculate mediation effect
all_beta = mr_results[8,7]  
all_beta # [1] 0.01176171

##Vitamin D (ukb-b-12648)
#Calculate Mediation proportion, calculated as the product of β1 and β2 divided by the total effect of the PRDM6 eQTLs on the chronic knee pain
XM_VitD_beta = res_exposure_to_mediator[43,7]  
XM_VitD_beta # [1] -0.001012102

MY_VitD_beta = res_mediator_to_outcome[38,7]  
MY_VitD_beta # [1] -7.852006

#mediation effect
aa_VitD = XM_VitD_beta * MY_VitD_beta  
aa_VitD # [1] 0.007947033

#Caculate 95%CI of the mediation effect,the 95% CIs of the mediation effect were calculated using the delta method
se_VitD <- sqrt((res_exposure_to_mediator[43,8]*MY_VitD_beta)^2 + (res_mediator_to_outcome[38,8]*XM_VitD_beta)^2)# [1] 0.004280345

#mediation proportions
bb_VitD = aa_VitD/all_beta 
bb_VitD# [1] 0.6756701
lower_CI_VitD <- bb_VitD - 0.97*(se_VitD/all_beta)
upper_CI_VitD <- bb_VitD + 0.97*(se_VitD/all_beta)

##NSAID(finn-b-RX_PARACETAMOL_NSAID)
XM_NSAID_beta = res_exposure_to_mediator[3,7]  
XM_NSAID_beta# [1] 0.02823064

MY_NSAID_beta = res_mediator_to_outcome[28,7]  
MY_NSAID_beta# [1] 0.3044713

#mediation effect
aa_NSAID = XM_NSAID_beta * MY_NSAID_beta  
aa_NSAID# [1] 0.00859542

#Caculate 95%CI of the mediation effect,the 95% CIs of the mediation effect were calculated using the delta method
se_NSAID <- sqrt((res_exposure_to_mediator[3,7]*MY_NSAID_beta)^2 + (res_mediator_to_outcome[28,7]*XM_NSAID_beta)^2)# [1] 0.004280345

#Caculate 95%CI of the mediation effect
bb_NSAID = aa_NSAID/all_beta 
bb_NSAID# [1] 0.730797
lower_CI_NSAID <- bb_NSAID - 0.97*(se_NSAID/all_beta) #0.97 is the critical value
upper_CI_NSAID <- aa_NSAID + 0.97*(se_NSAID/all_beta)

