# 3. Feature Construction

library(psych)
library(survival)
library(dplyr)

# --- 3.1 Factor Analysis ---
fa_result <- fa(protein_data_final, nfactors = 18, rotate = "oblimin", fm = "minres", maxit = 10000)
loadings_df <- as.data.frame(unclass(fa_result$loadings))

# --- 3.2 Top 5 proteins per factor by SS loading ---
top_per_factor <- list()
for (i in 1:ncol(loadings_df)) {
  factor_name <- paste0("Factor", i)
  factor_loadings <- loadings_df[[i]]
  names(factor_loadings) <- rownames(loadings_df)
  
  ss_loading <- factor_loadings^2 # squared loadings (variance explained per protein)
  top_idx <- order(ss_loading, decreasing = TRUE)[1:5]
  top_per_factor[[factor_name]] <- names(factor_loadings)[top_idx]
}

# --- 3.3 Cox regression for each of those proteins ---
cox_results <- data.frame(Factor = character(),
                          Protein = character(),
                          Loading = numeric(),
                          HR = numeric(),
                          lower95 = numeric(),
                          upper95 = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

for (fact in names(top_per_factor)) {
  fact_num <- as.numeric(gsub("Factor", "", fact))
  
  for (prot in top_per_factor[[fact]]) {
    loading_val <- loadings_df[prot, fact_num]
    
    fml <- as.formula(paste0("Surv(OS.time, OS) ~ `", prot, "`"))
    fit <- coxph(fml, data = protein_data_clean)
    summ <- summary(fit)
    
    cox_results <- rbind(cox_results, data.frame(
      Factor = fact,
      Protein = prot,
      Loading = loading_val,
      HR = summ$coef[,"exp(coef)"],
      lower95 = summ$conf.int[,"lower .95"],
      upper95 = summ$conf.int[,"upper .95"],
      p_value = summ$coef[,"Pr(>|z|)"],
      stringsAsFactors = FALSE
    ))
  }
}

# --- 3.4 Pick best per factor, resolve duplicates ---
best_per_factor <- cox_results %>%
  group_by(Factor) %>%
  arrange(p_value, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

while(any(duplicated(best_per_factor$Protein))) { # prevent duplicates
  dup_prots <- best_per_factor$Protein[duplicated(best_per_factor$Protein)]
  
  for (prot in dup_prots) {
    dup_rows <- best_per_factor %>% filter(Protein == prot)
    keep_row <- dup_rows %>% slice_min(p_value, n = 1, with_ties = FALSE)
    lose_factors <- setdiff(dup_rows$Factor, keep_row$Factor) # other factors lose duplicate
    
    for (lf in lose_factors) { # pick next best protein
      replacement <- cox_results %>%
        filter(Factor == lf, Protein != prot) %>%
        arrange(p_value) %>%
        slice(1)
      
      best_per_factor <- best_per_factor %>%
        filter(!(Factor == lf)) %>%
        bind_rows(replacement)
    }
  }
}

# --- 3.5 Final table: ordered by factor ---
final_top_proteins <- best_per_factor %>%
  mutate(HR_CI = paste0(sprintf("%.2f", HR),
                        " (", sprintf("%.2f", lower95),
                        "–", sprintf("%.2f", upper95), ")")) %>%
  select(Factor, Protein, Loading, HR_CI, p_value) %>%
  arrange(as.numeric(gsub("Factor", "", Factor)))

print(final_top_proteins)
