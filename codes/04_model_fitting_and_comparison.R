# 4. Model Fitting and Comparison
# run this line by line to define the functions

library(survival)
library(flexsurv)
library(randomForestSRC)
library(dplyr)

# --- 4.0 Helper functions ---

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

escape_tex <- function(x) {
  gsub("_", "\\\\_", x)
}

if (!dir.exists("tables")) {
  dir.create("tables")
}

aicc <- function(aic, k, n) {
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

clean_var_names <- function(x) {
  gsub("^genderMALE$", "gender", x)
}


# rounding

fmt_num <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), x)
}

fmt_z <- function(x) {
  sprintf("%.2f", x)
}

fmt_se <- function(x) {
  sprintf("%.3f", x)
}

fmt_ci <- function(est, low, high, digits = 2) {
  paste0(
    sprintf(paste0("%.", digits, "f"), est),
    " (",
    sprintf(paste0("%.", digits, "f"), low),
    ",",
    sprintf(paste0("%.", digits, "f"), high),
    ")"
  )
}

# LaTeX table writer

write_coef_table <- function(df, caption, label, effect_label) {
  
  rows <- apply(df, 1, function(x) {
    paste(
      x["Variable"], "&",
      x["Effect_CI"], "&",
      x["SE"], "&",
      x["z"], "&",
      x["p"], "\\\\"
    )
  })
  
  writeLines(c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    "\\begin{tabular}{",
    "    l",
    "    l",
    "    S[table-format=1.3]",
    "    S[table-format=-1.2]",
    "    l",
    "}",
    "\\toprule",
    paste0(
      "\\textbf{Variable} & \\textbf{", effect_label,
      "} & {\\textbf{SE}} & {\\textbf{z}} & {\\textbf{$p$}} \\\\"
    ),
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\multicolumn{5}{l}{\\textsl{Significance codes: p $<$ 0.001 (***), p $<$ 0.01 (**), p $<$ 0.05 (*)}} \\\\",
    "\\end{tabular}",
    paste0("\\label{", label, "}"),
    "\\end{table}"
  ), paste0("tables/", label, ".tex"))
}

# --- 4.1 Prepare modeling dataset ---

excluded_proteins <- c("ANNEXINVII", "YAP", "BAK", "EIF4G")
filtered_proteins <- setdiff(final_top_proteins$Protein, excluded_proteins)

predictor_vars <- c(filtered_proteins, "age", "gender")

df_model <- protein_data_clean[, c("OS.time", "OS", predictor_vars)]

df_model$gender <- as.factor(df_model$gender)
df_model <- df_model %>% mutate(across(where(is.character), as.factor))

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

df_model$age <- as.numeric(age_midpoints[as.character(df_model$age)])
df_model <- df_model %>% filter(!is.na(OS.time), OS.time > 0)

n <- nrow(df_model)

rhs_formula <- as.formula(
  paste("~", paste(predictor_vars, collapse = " + "))
)

cox_formula <- as.formula(
  paste("Surv(OS.time, OS) ~", paste(predictor_vars, collapse = " + "))
)

# --- 4.2 Cox proportional hazards model ---

cox_fit <- coxph(cox_formula, data = df_model, ties = "efron")
cox_sum <- summary(cox_fit)
cox_ci  <- confint(cox_fit)

cox_tbl <- data.frame(
  Variable = escape_tex(clean_var_names(rownames(cox_sum$coef))),
  Effect_CI = fmt_ci(
    exp(cox_sum$coef[, "coef"]),
    exp(cox_ci[, 1]),
    exp(cox_ci[, 2]),
    digits = 2
  ),
  SE = fmt_se(cox_sum$coef[, "se(coef)"]),
  z  = fmt_z(cox_sum$coef[, "z"]),
  p  = fmt_p(cox_sum$coef[, "Pr(>|z|)"]),
  stringsAsFactors = FALSE
)


write_coef_table(
  cox_tbl,
  caption = "Cox proportional hazards model results",
  label = "cox_full_results",
  effect_label = "HR (95\\% CI)"
)

# --- 4.3 Parametric survival models ---

weib_fit   <- flexsurvreg(cox_formula, data = df_model, dist = "weibull")
llog_fit   <- flexsurvreg(cox_formula, data = df_model, dist = "llogis")
ggamma_fit <- flexsurvreg(cox_formula, data = df_model, dist = "gengamma")

# AFT table extractor (ETR)

extract_aft_table <- function(fit) {
  
  res <- fit$res
  X <- model.matrix(rhs_formula, df_model)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  res <- res[rownames(res) %in% colnames(X), , drop = FALSE]
  
  est <- res[, "est"]
  se  <- res[, "se"]
  z   <- est / se
  p   <- fmt_p(2 * (1 - pnorm(abs(z))))
  
  etr  <- exp(-est)
  low  <- exp(-(est + 1.96 * se))
  high <- exp(-(est - 1.96 * se))
  
  data.frame(
    Variable = escape_tex(clean_var_names(rownames(res))),
    Effect_CI = fmt_ci(etr, low, high, digits = 2),
    SE = fmt_se(se),
    z  = fmt_z(z),
    p  = p,
    stringsAsFactors = FALSE
  )
}


write_coef_table(
  extract_aft_table(weib_fit),
  caption = "Weibull accelerated failure-time model results",
  label = "weibull_full_results",
  effect_label = "ETR (95\\% CI)"
)

write_coef_table(
  extract_aft_table(llog_fit),
  caption = "Log-logistic accelerated failure-time model results",
  label = "loglogistic_full_results",
  effect_label = "ETR (95\\% CI)"
)

# Generalized gamma table extractor (coef)

extract_gengamma_table <- function(fit) {
  
  res <- fit$res
  X <- model.matrix(rhs_formula, df_model)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  res <- res[rownames(res) %in% colnames(X), , drop = FALSE]
  
  est <- res[, "est"]
  se  <- res[, "se"]
  z   <- est / se
  p   <- fmt_p(2 * (1 - pnorm(abs(z))))
  
  data.frame(
    Variable = escape_tex(clean_var_names(rownames(res))),
    Effect_CI = fmt_ci(
      est,
      est - 1.96 * se,
      est + 1.96 * se,
      digits = 3   # <-- usually 3 for raw coefficients
    ),
    SE = fmt_se(se),
    z  = fmt_z(z),
    p  = p,
    stringsAsFactors = FALSE
  )
}


write_coef_table(
  extract_gengamma_table(ggamma_fit),
  caption = "Generalized gamma survival model results",
  label = "gengamma_full_results",
  effect_label = "Coef (95\\% CI)"
)

# --- 4.4 Random survival forest ---

set.seed(123)

rsf_fit <- rfsrc(
  formula = cox_formula,
  data = df_model,
  ntree = 1000,
  mtry = floor(sqrt(length(predictor_vars))),
  nodesize = 15,
  splitrule = "logrank"
)

rsf_pred <- predict(rsf_fit)
rsf_risk <- rsf_pred$chf[, ncol(rsf_pred$chf)]

rsf_cindex <- 1 - concordance(
  Surv(OS.time, OS) ~ rsf_risk,
  data = df_model
)$concordance

# --- 4.5 Model comparison table ---

# c index helper

get_cindex_flexsurv <- function(fit, df, rhs_formula) {
  
  X <- model.matrix(rhs_formula, df)
  
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  
  beta_names <- intersect(colnames(X), rownames(fit$res))
  beta <- fit$res[beta_names, "est"]
  
  X <- X[, beta_names, drop = FALSE]
  
  lp <- as.numeric(X %*% beta)
  
  concordance(Surv(OS.time, OS) ~ lp, data = df)$concordance
}


get_k <- function(fit) nrow(fit$res)

model_results <- data.frame(
  Model = c(
    "Cox PH",
    "Weibull",
    "Log-logistic",
    "Generalized gamma",
    "Random Survival Forest"
  ),
  AIC = c( # cox AIC, AICc, BIC = NA as not comparable
    NA,
    AIC(weib_fit),
    AIC(llog_fit),
    AIC(ggamma_fit),
    NA
  ),
  AICc = c(
    NA,
    aicc(AIC(weib_fit), get_k(weib_fit), n),
    aicc(AIC(llog_fit), get_k(llog_fit), n),
    aicc(AIC(ggamma_fit), get_k(ggamma_fit), n),
    NA
  ),
  BIC = c(
    NA,
    BIC(weib_fit),
    BIC(llog_fit),
    BIC(ggamma_fit),
    NA
  ),
  C_index = c(
    cox_sum$concordance[1],
    get_cindex_flexsurv(weib_fit, df_model, rhs_formula),
    get_cindex_flexsurv(llog_fit, df_model, rhs_formula),
    get_cindex_flexsurv(ggamma_fit, df_model, rhs_formula),
    rsf_cindex
  ),
  stringsAsFactors = FALSE
)

model_results_fmt <- model_results %>%
  mutate(
    AIC  = ifelse(is.na(AIC), "---", sprintf("%.2f", AIC)),
    AICc = ifelse(is.na(AICc), "---", sprintf("%.2f", AICc)),
    BIC  = ifelse(is.na(BIC), "---", sprintf("%.2f", BIC)),
    C_index = ifelse(is.na(C_index), "---", sprintf("%.3f", C_index))
  )

rows <- apply(model_results_fmt, 1, function(x) {
  paste(x["Model"], "&", x["AIC"], "&", x["AICc"], "&", x["BIC"], "&", x["C_index"], "\\\\")
})

# save LaTeX pdf 

writeLines(c(
  "\\begin{table}[H]",
  "\\centering",
  "\\begin{tabular}{lrrrr}",
  "\\hline",
  "Model & AIC & AICc & BIC & C-index \\\\",
  "\\hline",
  rows,
  "\\hline",
  "\\end{tabular}",
  "\\caption{Comparison of survival models by information criteria and concordance index.}",
  "\\label{tab:fullc}",
  "\\end{table}"
), "tables/full_model_comparison.tex")

