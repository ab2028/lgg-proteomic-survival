# 6. Univariate and Multivariable Cox comparison table

library(survival)
library(dplyr)
library(broom)

proteins <- c(final_top_proteins$Protein, "age", "gender")

age_midpoints <- c(
  "10-20" = 15,
  "21-30" = 25.5,
  "31-40" = 35.5,
  "41-50" = 45.5,
  "51-60" = 55.5,
  "61-70" = 65.5,
  "71-80" = 75.5,
  "81-90" = 85.5
)

protein_data_c2 <- protein_data_clean
protein_data_c2$age <- age_midpoints[as.character(protein_data_c2$age)]
protein_data_c2$age <- as.numeric(protein_data_c2$age)


#1
univ_results <- lapply(proteins, function(prot) {
  
  fml <- as.formula(paste0("Surv(OS.time, OS) ~ `", prot, "`"))
  fit <- coxph(fml, data = protein_data_c2)
  summ <- summary(fit)
  
  data.frame(
    Protein = prot,
    HR = summ$coef[,"exp(coef)"],
    lower95 = summ$conf.int[,"lower .95"],
    upper95 = summ$conf.int[,"upper .95"],
    p_value = summ$coef[,"Pr(>|z|)"],
    Concordance = summ$concordance[1]
  )
})

univ_cox_table <- bind_rows(univ_results) %>%
  mutate(
    HR_CI = paste0(
      sprintf("%.2f", HR), 
      " (", sprintf("%.2f", lower95), 
      "–", sprintf("%.2f", upper95), ")"
    )
  ) %>%
  select(Protein, HR_CI, p_value, Concordance)

print(univ_cox_table)


# 2

multi_formula <- as.formula(
  paste(
    "Surv(OS.time, OS) ~",
    paste(predictor_vars, collapse = " + ")
  )
)

multi_fit <- coxph(multi_formula, data = df_model, ties = "efron")
multi_summary <- summary(multi_fit)


multi_cox_table <- data.frame(
  Variable = rownames(multi_summary$coef),
  HR = multi_summary$coef[,"exp(coef)"],
  lower95 = multi_summary$conf.int[,"lower .95"],
  upper95 = multi_summary$conf.int[,"upper .95"],
  p_value = multi_summary$coef[,"Pr(>|z|)"]
) %>%
  mutate(
    HR_CI = paste0(
      sprintf("%.2f", HR),
      " (", sprintf("%.2f", lower95),
      "–", sprintf("%.2f", upper95), ")"
    )
  ) %>%
  select(Variable, HR_CI, p_value)

print(multi_cox_table)


c_index <- multi_summary$concordance
print(paste0(
  "C-index: ", round(c_index[1], 3),
  " (SE = ", round(c_index[2], 3), ")"
))

global_tests <- data.frame(
  Test = c("Likelihood ratio", "Wald", "Score (log-rank)"),
  Statistic = c(
    multi_summary$logtest["test"],
    multi_summary$waldtest["test"],
    multi_summary$sctest["test"]
  ),
  df = c(
    multi_summary$logtest["df"],
    multi_summary$waldtest["df"],
    multi_summary$sctest["df"]
  ),
  p_value = c(
    multi_summary$logtest["pvalue"],
    multi_summary$waldtest["pvalue"],
    multi_summary$sctest["pvalue"]
  )
)

print(global_tests)


ph_test <- cox.zph(multi_fit)
print(ph_test)


