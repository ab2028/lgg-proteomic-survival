# 6. Univariate and Multivariable Cox comparison table

library(survival)
library(dplyr)
library(broom)


# --- 6.1 Prepare Modeling Dataset ---

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


# --- 6.2 Univariate Table ---
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


# --- 6.3 Multivariate Table ---

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


# --- 6.4 PH Test (optional) ---

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

# print(global_tests)


ph_test <- cox.zph(multi_fit)
# print(ph_test)

# --- 6.5 Save combined univariate + multivariable Cox table in LaTeX ---

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

fmt_hr_ci_p <- function(hr, low, high, p) {
  paste0(
    sprintf("%.2f", hr),
    " (", sprintf("%.2f", low), "--", sprintf("%.2f", high), "), ",
    fmt_p(p)
  )
}

escape_tex <- function(x) {
  gsub("_", "\\\\_", x)
}

# --- rebuild UNIVARIATE table from raw univ_results ---

univ_raw <- bind_rows(univ_results)

univ_fmt <- univ_raw %>%
  mutate(
    Variable = escape_tex(Protein),
    Univ = fmt_hr_ci_p(HR, lower95, upper95, p_value)
  ) %>%
  select(Variable, Univ)

# --- rebuild multivariable table from multi_summary ---

multi_raw <- data.frame(
  Variable = rownames(multi_summary$coef),
  HR = multi_summary$coef[,"exp(coef)"],
  lower95 = multi_summary$conf.int[,"lower .95"],
  upper95 = multi_summary$conf.int[,"upper .95"],
  p_value = multi_summary$coef[,"Pr(>|z|)"],
  stringsAsFactors = FALSE
)

multi_fmt <- multi_raw %>%
  mutate(
    Variable = escape_tex(Variable),
    Multi = fmt_hr_ci_p(HR, lower95, upper95, p_value)
  ) %>%
  select(Variable, Multi)

# --- combine, preserve univariate order ---

combined_table <- univ_fmt %>%
  left_join(multi_fmt, by = "Variable") %>%
  mutate(
    Multi = ifelse(is.na(Multi), "--", Multi)
  )

# --- global test stats ---

c_index_val <- round(multi_summary$concordance[1], 3)
c_index_se  <- round(multi_summary$concordance[2], 3)

lrt   <- multi_summary$logtest
wald  <- multi_summary$waldtest
score <- multi_summary$sctest

# --- write and save LaTeX table ---

rows <- apply(combined_table, 1, function(x) {
  paste(x["Variable"], "&", x["Univ"], "&", x["Multi"], "\\\\")
})

writeLines(c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Univariate and full multivariable Cox proportional hazards regression results.}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Variable & Univariate HR (95\\% CI), $p$ & Full multivariable HR (95\\% CI), $p$ \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{0.5em}",
  "\\footnotesize{",
  paste0(
    "Full multivariable model fit statistics: Concordance = ",
    sprintf("%.3f", c_index_val),
    " (SE = ",
    sprintf("%.3f", c_index_se),
    ")."
  ),
  paste0(
    " Likelihood ratio test = ",
    round(lrt["test"], 1),
    " on ",
    lrt["df"],
    " df ($p < 0.001$);"
  ),
  paste0(
    " Wald test = ",
    round(wald["test"], 1),
    " on ",
    wald["df"],
    " df ($p < 0.001$);"
  ),
  paste0(
    " Score (log-rank) test = ",
    round(score["test"], 1),
    " on ",
    score["df"],
    " df ($p < 0.001$)."
  ),
  "}",
  "\\footnotesize{",
  "Note: YAP, BAK, EIF4G, and Annexin A7 were excluded from multivariable modeling due to violations of the proportional hazards assumption (see Section 3.1 of paper)",
  "}",
  "\\label{tab:cox_all_combined}",
  "\\end{table}"
), "tables/table_cox_all_combined.tex")



