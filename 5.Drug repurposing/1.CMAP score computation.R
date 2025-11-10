library(qvalue)

#Cmap_score function
cmap_score <- function(sig_up, sig_down, drug_signature) {
  #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this, not fully validated
  #Note. I think that creating the anonymous functions in each iteration of the sapply's below is slowing things down. Predefine them eventually.
  num_genes <- nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  
  # I think we are re-ranking because the GeneID mapping changed the original rank range
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
  
  # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
  
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  
  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  
  # 
  if(num_tags_up > 1) {
    a_up <- 0
    b_up <- 0
    
    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
  }else{
    ks_up <- 0
  }
  
  if (num_tags_down > 1){
    
    a_down <- 0
    b_down <- 0
    
    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
  }else{
    ks_down <- 0
  }
  
  if (ks_up == 0 & ks_down != 0){ #only down gene inputed
    connectivity_score <- -ks_down 
  }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
    connectivity_score <- ks_up
  }else if (sum(sign(c(ks_down,ks_up))) == 0) {
    connectivity_score <- ks_up - ks_down # different signs
  }
  
  return(connectivity_score)
}

# Define subtypes (subtype1 to subtype6)
subtypes <- paste0("TC", 1:6)

# Parameters
N_PERMUTATIONS <- 100000  # Default: 100,000 permutations
landmark <- 1             # Assuming landmark genes are used

for (subtype in subtypes) {
  cat("\n===== Processing", subtype, "=====\n")
  #===========================================
  # 1. Compute CMap scores for real signatures
  #===========================================
  cat("\nComputing CMap scores for", subtype, "...\n")
  dz_cmap_scores <- NULL
  count <- 0
  
  # Load subtype-specific DEGs
  dz_genes_up <- get(paste0(subtype, "_dz_genes_up"))  
  dz_genes_down <- get(paste0(subtype, "_dz_genes_down"))  
  
  for (exp_id in rownames(drug_info)) {
    count <- count + 1
    if (landmark == 1) {
      cmap_exp_signature <- data.frame(
        gene_name, 
        rank(-1 * drug_sig[, as.character(exp_id)], ties.method = "random")
      )
    } else {
      id <- which(colnames(drug_sig) == exp_id)
      cmap_exp_signature <- data.frame(gene_name, drug_sig[, id])
    }
    
    colnames(cmap_exp_signature) <- c("ids", "rank")
    dz_cmap_scores <- c(dz_cmap_scores, cmap_score(dz_genes_up,dz_genes_down, cmap_exp_signature))
  }
  
  assign(paste0(subtype, "_dz_cmap_scores"), dz_cmap_scores)
  save(list = paste0(subtype, "_dz_cmap_scores"), file = paste0("./result/", subtype, "_dz_cmap_scores.RData"))
  
  #===========================================
  # 2. Compute the significance
  #===========================================
  cat("\nComputing the significance", subtype, "...\n")
  random_sig_ids <- sample(colnames(drug_sig), N_PERMUTATIONS, replace = TRUE)
  random_cmap_scores <- NULL
  count <- 0
  
  for (exp_id in random_sig_ids) {
    count <- count + 1
    id <- which(colnames(drug_sig) == exp_id)
    if (landmark == 1) {
      cmap_exp_signature_randoms <- data.frame(
        gene_name, 
        rank(-1 * drug_sig[, id], ties.method = "random")
      )
    } else {
      cmap_exp_signature_randoms <- data.frame(gene_name, drug_sig[, id])
    }
    colnames(cmap_exp_signature_randoms) <- c("ids", "rank")
    random_input_genes <- sample(gene_name, nrow(dz_genes_up) + nrow(dz_genes_down))
    rand_dz_gene_up <- data.frame(GeneID = random_input_genes[1:nrow(dz_genes_up)])
    rand_dz_gene_down <- data.frame(GeneID = random_input_genes[(nrow(dz_genes_up) + 1):length(random_input_genes)])
    random_cmap_scores <- c(random_cmap_scores,cmap_score(rand_dz_gene_up, rand_dz_gene_down, cmap_exp_signature_randoms))
    rm(random_input_genes,rand_dz_gene_up,rand_dz_gene_down)
    }
  
  # Save random scores
  assign(paste0(subtype, "_random_cmap_scores"), random_cmap_scores)
  save(list = paste0(subtype, "_random_cmap_scores"), 
       file = paste0("./result/", subtype, "_random_cmap_scores.RData"))
  
  #===========================================
  # 3. Compute p-values and q-values
  #===========================================
  cat("\nCalculating p/q-values for", subtype, "...\n")
  random_scores <- unlist(random_cmap_scores)
  
  # Two-tailed p-value
  p_values <- sapply(dz_cmap_scores, function(score) {
    sum(abs(random_scores) >= abs(score)) / length(random_scores)
  })
  
  # FDR correction
  q_values <- qvalue(p_values)$qvalues
  
  # Save results
  assign(paste0("p_values_", subtype), p_values)
  assign(paste0("q_values_", subtype), q_values)
  
  #===========================================
  # 4. Merge with drug metadata
  #===========================================
  drugs <- data.frame(
    exp_id = rownames(drug_info),
    cmap_score = dz_cmap_scores,
    p = p_values,
    q = q_values
  )
  
  drug_instances_all <- merge(drugs, drug_info, by.x = "exp_id", by.y = 0)
  assign(paste0(subtype, "_drug_instances_all"), drug_instances_all)
  rm(dz_genes_up,dz_genes_down,dz_cmap_scores,random_cmap_scores,p_values,q_values,drugs,drug_instances_all)
  cat("\nCompleted", subtype, "!\n")
}