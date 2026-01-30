# 7. Backwards Elimination

library(survival)
library(flexsurv)
library(dplyr)

# --- 7.1 ---

# force age and gender
backward_eliminate <- function(model_type, df, predictor_vars, alpha = 0.05,
                               forced_vars = c("age", "gender")) {
  
  # store original names
  original_vars <- predictor_vars
  
  # create safe names for modeling
  safe_names <- make.names(predictor_vars)
  names(safe_names) <- predictor_vars
  
  # forced variables must always stay
  forced_safe <- safe_names[names(safe_names) %in% forced_vars]
  
  # current working set
  current_safe <- safe_names
  
  repeat {
    
    # build formula with safe names
    surv_formula <- as.formula(
      paste("Surv(OS.time, OS) ~", paste(current_safe, collapse = " + "))
    )
    
    # fit model depending on family
    if (model_type == "cox") {
      fit <- coxph(surv_formula, data = df, ties = "efron")
      sm <- summary(fit)
      coef_table <- sm$coefficients
      
      p_vals <- coef_table[, "Pr(>|z|)"]
      names(p_vals) <- rownames(coef_table)
      
    } else {
      fit <- flexsurvreg(surv_formula, data = df, dist = model_type)
      
      # determine non covariate parameters to drop
      if (model_type == "weibull") drop_params <- c("shape", "scale")
      else if (model_type == "llogis") drop_params <- c("shape", "scale")
      else if (model_type == "gengamma") drop_params <- c("shape", "sigma", "Q", "mu") 
      
      res <- fit$res[!(rownames(fit$res) %in% drop_params), , drop = FALSE]
      
      # compute p values if missing
      if (!"p" %in% colnames(res)) {
        z_val <- res[, "est"] / res[, "se"]
        p_vals <- 2 * (1 - pnorm(abs(z_val)))
        names(p_vals) <- rownames(res)
      } else {
        p_vals <- res[, "p"]
      }
    }
    
    # handle NA p values which indicate aliasing or dropped terms
    if (any(is.na(p_vals))) {
      na_name <- names(p_vals)[which(is.na(p_vals))[1]]
      message("Dropping term with NA p value ", na_name)
      
      # map back to original
      keep_original <- names(current_safe)[current_safe == na_name]
      
      # do not allow dropping forced variables
      if (keep_original %in% forced_vars) {
        stop("Model attempted to drop a forced variable. Check predictors.")
      }
      
      current_safe <- current_safe[current_safe != na_name]
      if (length(current_safe) == length(forced_safe)) break
      next
    }
    
    # sort variables by p value from worst to best
    ordered_names <- names(sort(p_vals, decreasing = TRUE))
    
    # find first non forced variable in the list
    worst_name <- NULL
    worst_original <- NULL
    
    for (nm in ordered_names) {
      original_candidate <- names(safe_names)[safe_names == nm]
      
      # skip NA or unmapable names
      if (length(original_candidate) != 1) next
      
      # skip forced variables
      if (original_candidate %in% forced_vars) next
      
      worst_name <- nm
      worst_original <- original_candidate
      break
    }
    
    # if no non forced vars remain above threshold stop
    if (is.null(worst_name)) {
      message("Only forced variables remain. Stopping elimination.")
      break
    }
    
    worst_p <- p_vals[worst_name]
    
    # stopping rule
    if (worst_p <= alpha) {
      message("All remaining non forced variables significant at threshold ", alpha)
      break
    }
    
    message("Dropping ", worst_original, " (p ", round(worst_p, 5), ")")
    
    # drop the worst variable
    current_safe <- current_safe[current_safe != worst_name]
    
    # stop if nothing except forced vars remains
    if (length(current_safe) == length(forced_safe)) break
    
  }
  
  # final mapping back
  final_original <- names(safe_names)[safe_names %in% current_safe]
  
  # refit final model
  final_formula <- as.formula(
    paste("Surv(OS.time, OS) ~", paste(current_safe, collapse = " + "))
  )
  
  final_fit <- if (model_type == "cox") {
    coxph(final_formula, data = df, ties = "efron")
  } else {
    flexsurvreg(final_formula, data = df, dist = model_type)
  }
  
  # produce clean summary table
  if (model_type == "cox") {
    sm <- summary(final_fit)
    tbl <- data.frame(
      coef = sm$coefficients[, "coef"],
      HR = sm$coefficients[, "exp(coef)"],
      se = sm$coefficients[, "se(coef)"],
      z = sm$coefficients[, "z"],
      p = sm$coefficients[, "Pr(>|z|)"]
    )
    tbl$Variable <- rownames(tbl)
  } else {
    res <- final_fit$res[!(rownames(final_fit$res) %in% drop_params), , drop = FALSE]
    z_val <- res[, "est"] / res[, "se"]
    p_val <- 2 * (1 - pnorm(abs(z_val)))
    
    tbl <- data.frame(
      coef = res[, "est"],
      HR = exp(res[, "est"]),
      se = res[, "se"],
      z = z_val,
      p = p_val,
      Variable = rownames(res)
    )
  }
  
  # compute model metrics
  if (model_type == "cox") {
    c_index <- summary(final_fit)$concordance[1]
  } else {
    X <- model.matrix(final_formula, data = df)[, -1, drop = FALSE]
    coefs <- final_fit$res[colnames(X), "est"]
    lp <- as.numeric(X %*% coefs)
    c_index <- concordance(Surv(OS.time, OS) ~ lp, data = df)$concordance
  }
  
  aic_val <- AIC(final_fit)
  k_val <- length(current_safe)
  n_val <- nrow(df)
  aicc_val <- aic_val + (2 * k_val * (k_val + 1)) / (n_val - k_val - 1)
  bic_val <- BIC(final_fit)
  
  list(
    final_vars = final_original,
    final_fit = final_fit,
    summary_table = tbl,
    model_info = data.frame(
      Model = model_type,
      AIC = aic_val,
      AICc = aicc_val,
      BIC = bic_val,
      C_index = c_index
    )
  )
}

#newcode
df_model <- df_model %>%
  mutate(
    gender = as.integer(gender == "MALE")
  )


# predictor_vars <- predictor_vars[predictor_vars != "genderMALE"]
# predictor_vars <- c(predictor_vars, "gender")



result_cox <- backward_eliminate(
  model_type = "cox",
  df = df_model,
  predictor_vars = predictor_vars,
  alpha = 0.05,
  forced_vars = c("age", "gender")
)

result_weibull <- backward_eliminate(
  model_type = "weibull",
  df = df_model,
  predictor_vars = predictor_vars,
  alpha = 0.05,
  forced_vars = c("age", "gender")
)

result_llogis <- backward_eliminate(
  model_type = "llogis",
  df = df_model,
  predictor_vars = predictor_vars,
  alpha = 0.05,
  forced_vars = c("age", "gender")
)

result_gengamma <- backward_eliminate(
  model_type = "gengamma",
  df = df_model,
  predictor_vars = predictor_vars,
  alpha = 0.05,
  forced_vars = c("age", "gender")
)

# Final model comparison table for reduced models
# Uses results from backward_eliminate for Cox, Weibull, Loglogistic, Generalized Gamma

library(dplyr)

# function for AICc
aicc <- function(aic_value, k_value, n_value) {
  aic_value + (2 * k_value * (k_value + 1)) / (n_value - k_value - 1)
}

# extract model and sample size
n_value <- nrow(df_model)

# construct table manually from your result objects
final_model_results <- data.frame(
  Model = c("Cox PH", "Weibull", "Log-logistic", "Generalized gamma"),
  
  AIC = c(
    AIC(result_cox$final_fit),
    AIC(result_weibull$final_fit),
    AIC(result_llogis$final_fit),
    AIC(result_gengamma$final_fit)
  ),
  
  AICc = c(
    aicc(AIC(result_cox$final_fit), length(result_cox$final_vars), n_value),
    aicc(AIC(result_weibull$final_fit), length(result_weibull$final_vars), n_value),
    aicc(AIC(result_llogis$final_fit), length(result_llogis$final_vars), n_value),
    aicc(AIC(result_gengamma$final_fit), length(result_gengamma$final_vars), n_value)
  ),
  
  BIC = c(
    BIC(result_cox$final_fit),
    BIC(result_weibull$final_fit),
    BIC(result_llogis$final_fit),
    BIC(result_gengamma$final_fit)
  ),
  
  C_index = c(
    result_cox$model_info$C_index,
    result_weibull$model_info$C_index,
    result_llogis$model_info$C_index,
    result_gengamma$model_info$C_index
  ),
  
  stringsAsFactors = FALSE
)

print(final_model_results)

# MODEL OUTPUTS

if (TRUE) {

cox_fit <- result_cox$final_fit
cox_sum <- summary(cox_fit)
cox_ci  <- confint(cox_fit)

cox_table <- data.frame(
  Variable = rownames(cox_sum$coef),
  `HR (95% CI)` = sprintf(
    "%.3f (%.3f–%.3f)",
    exp(cox_sum$coef[, "coef"]),
    exp(cox_ci[, 1]),
    exp(cox_ci[, 2])
  ),
  SE = cox_sum$coef[, "se(coef)"],
  z  = cox_sum$coef[, "z"],
  p  = cox_sum$coef[, "Pr(>|z|)"],
  row.names = NULL,
  check.names = FALSE
)

print(cox_table)


weib_fit <- result_weibull$final_fit
weib_res <- weib_fit$res[!rownames(weib_fit$res) %in% c("shape", "scale"), , drop = FALSE]

z_val <- weib_res[, "est"] / weib_res[, "se"]
p_val <- 2 * (1 - pnorm(abs(z_val)))

ci_low  <- weib_res[, "est"] - 1.96 * weib_res[, "se"]
ci_high <- weib_res[, "est"] + 1.96 * weib_res[, "se"]

weibull_table <- data.frame(
  Variable = rownames(weib_res),
  `HR (95% CI)` = sprintf(
    "%.3f (%.3f–%.3f)",
    exp(weib_res[, "est"]),
    exp(ci_low),
    exp(ci_high)
  ),
  SE = weib_res[, "se"],
  z  = z_val,
  p  = p_val,
  row.names = NULL,
  check.names = FALSE
  
)

print(weibull_table)






llog_fit <- result_llogis$final_fit
llog_res <- llog_fit$res[!rownames(llog_fit$res) %in% c("shape", "scale"), , drop = FALSE]

z_val <- llog_res[, "est"] / llog_res[, "se"]
p_val <- 2 * (1 - pnorm(abs(z_val)))

ci_low  <- llog_res[, "est"] - 1.96 * llog_res[, "se"]
ci_high <- llog_res[, "est"] + 1.96 * llog_res[, "se"]

loglogistic_table <- data.frame(
  Variable = rownames(llog_res),
  `HR (95% CI)` = sprintf(
    "%.3f (%.3f–%.3f)",
    exp(llog_res[, "est"]),
    exp(ci_low),
    exp(ci_high)
  ),
  SE = llog_res[, "se"],
  z  = z_val,
  p  = p_val,
  row.names = NULL,
  check.names = FALSE
  
)

print(loglogistic_table)


ggamma_fit <- result_gengamma$final_fit
ggamma_res <- ggamma_fit$res[
  !rownames(ggamma_fit$res) %in% c("shape", "sigma", "Q", "mu"),
  , drop = FALSE
]

z_val <- ggamma_res[, "est"] / ggamma_res[, "se"]
p_val <- 2 * (1 - pnorm(abs(z_val)))

ci_low  <- ggamma_res[, "est"] - 1.96 * ggamma_res[, "se"]
ci_high <- ggamma_res[, "est"] + 1.96 * ggamma_res[, "se"]

gengamma_table <- data.frame(
  Variable = rownames(ggamma_res),
  `Coef (95% CI)` = sprintf(
    "%.3f (%.3f–%.3f)",
    ggamma_res[, "est"],
    ci_low,
    ci_high
  ),
  SE = ggamma_res[, "se"],
  z  = z_val,
  p  = p_val,
  row.names = NULL,
  check.names = FALSE
  
)

print(gengamma_table)

}



# LATEX TABLE SAVING

if (!dir.exists("tables")) dir.create("tables")


sig_stars <- function(p) {
  ifelse(
    p < 0.001, "{***}",
    ifelse(p < 0.01, "{**}",
           ifelse(p < 0.05, "{*}", ""))
  )
}


write_reduced_cox_table <- function(result, caption, label) {
  
  fit <- result$final_fit
  sm  <- summary(fit)
  ci  <- confint(fit)
  
  df <- data.frame(
    Variable = escape_tex(rownames(sm$coef)),
    Effect_CI = sprintf(
      "%.3f (%.3f, %.3f)",
      exp(sm$coef[, "coef"]),
      exp(ci[, 1]),
      exp(ci[, 2])
    ),
    SE = sprintf("%.3f", sm$coef[, "se(coef)"]),
    z  = sprintf("%.2f", sm$coef[, "z"]),
    p  = ifelse(
      sm$coef[, "Pr(>|z|)"] < 0.001,
      "$<0.001$",
      sprintf("%.3f", sm$coef[, "Pr(>|z|)"])
    ),
    star = sig_stars(sm$coef[, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
  
  rows <- apply(df, 1, function(x) {
    paste(
      x["Variable"], "&",
      x["Effect_CI"], "&",
      x["SE"], "&",
      x["z"], "&",
      paste0(x["p"], " ", x["star"]),
      "\\\\"
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
    "\\textbf{Variable} & \\textbf{HR (95\\% CI)} & {\\textbf{SE}} & {\\textbf{z}} & {\\textbf{Pr$(>|z|)$}} \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\multicolumn{5}{l}{\\textsl{Significance codes: p $<$ 0.001 (***), p $<$ 0.01 (**), p $<$ 0.05 (*)}} \\\\",
    "\\end{tabular}",
    paste0("\\label{", label, "}"),
    "\\end{table}"
  ), paste0("tables/", label, ".tex"))
}

write_reduced_aft_table <- function(result, caption, label) {
  
  fit <- result$final_fit
  res <- fit$res
  
  # use the retained variables from backward elimination
  covars <- result$final_vars
  
  # keep only covariate rows
  res <- res[rownames(res) %in% covars, , drop = FALSE]
  
  if (nrow(res) == 0) {
    stop("No covariates found for AFT table: ", label)
  }
  
  z <- res[, "est"] / res[, "se"]
  p <- 2 * (1 - pnorm(abs(z)))
  
  etr  <- exp(res[, "est"])
  low  <- exp(-(res[, "est"] + 1.96 * res[, "se"]))
  high <- exp(-(res[, "est"] - 1.96 * res[, "se"]))
  
  df_out <- data.frame(
    Variable = escape_tex(rownames(res)),
    Effect_CI = sprintf("%.3f (%.3f--%.3f)", etr, low, high),
    SE = sprintf("%.3f", res[, "se"]),
    z  = sprintf("%.2f", z),
    p  = ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p)),
    star = sig_stars(p),
    stringsAsFactors = FALSE
  )
  
  rows <- apply(df_out, 1, function(x) {
    paste(
      x["Variable"], "&",
      x["Effect_CI"], "&",
      x["SE"], "&",
      x["z"], "&",
      paste0(x["p"], " ", x["star"]),
      "\\\\"
    )
  })
  
  writeLines(c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    "\\begin{tabular}{l l S[table-format=1.3] S[table-format=-1.2] l}",
    "\\toprule",
    "\\textbf{Variable} & \\textbf{ETR (95\\% CI)} & {\\textbf{SE}} & {\\textbf{z}} & {\\textbf{Pr$(>|z|)$}} \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\multicolumn{5}{l}{\\parbox[t]{\\linewidth}{",
    "\\textsl{Interpretation: ETR $>$ 1 indicates longer survival}\\\\",
    "\\textsl{Significance codes: *** $p<0.001$, ** $p<0.01$, * $p<0.05$}",
    "}} \\\\",
    "\\end{tabular}",
    paste0("\\label{", label, "}"),
    "\\end{table}"
  ), paste0("tables/", label, ".tex"))
}



write_reduced_gengamma_table <- function(result, caption, label) {
  
  fit <- result$final_fit
  res <- fit$res
  
  drop_params <- intersect(
    rownames(res),
    c("shape", "sigma", "Q", "mu")
  )
  res <- res[!rownames(res) %in% drop_params, , drop = FALSE]
  
  z <- res[, "est"] / res[, "se"]
  p <- 2 * (1 - pnorm(abs(z)))
  
  df <- data.frame(
    Variable = escape_tex(rownames(res)),
    Effect_CI = sprintf(
      "%.3f (%.3f, %.3f)",
      res[, "est"],
      res[, "est"] - 1.96 * res[, "se"],
      res[, "est"] + 1.96 * res[, "se"]
    ),
    SE = sprintf("%.3f", res[, "se"]),
    z  = sprintf("%.2f", z),
    p  = ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p)),
    star = sig_stars(p),
    stringsAsFactors = FALSE
  )
  
  rows <- apply(df, 1, function(x) {
    paste(
      x["Variable"], "&",
      x["Effect_CI"], "&",
      x["SE"], "&",
      x["z"], "&",
      paste0(x["p"], " ", x["star"]),
      "\\\\"
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
    "\\textbf{Variable} & \\textbf{Coef (95\\% CI)} & {\\textbf{SE}} & {\\textbf{z}} & {\\textbf{Pr$(>|z|)$}} \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\multicolumn{5}{l}{\\textsl{Significance codes: p $<$ 0.001 (***), p $<$ 0.01 (**), p $<$ 0.05 (*)}} \\\\",
    "\\end{tabular}",
    paste0("\\label{", label, "}"),
    "\\end{table}"
  ), paste0("tables/", label, ".tex"))
}



write_reduced_cox_table(
  result_cox,
  "Reduced Cox proportional hazards model results",
  "reduced_cox"
)

write_reduced_aft_table(
  result_llogis,
  "Reduced log-logistic accelerated failure-time model results",
  "reduced_loglogistic"
)

write_reduced_aft_table(
  result_weibull,
  "Reduced Weibull accelerated failure-time model results",
  "reduced_weibull"
)

write_reduced_gengamma_table(
  result_gengamma,
  "Reduced generalized gamma survival model results",
  "reduced_gengamma"
)






