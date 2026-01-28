# 4. Model Fitting and Comparison

library(survival)
library(flexsurv)
library(rms)
library(randomForestSRC)
library(dplyr)

# --- 4.1 Prepare modeling dataset ---

# proteins excluded due to PH violations
excluded_proteins <- c("ANNEXINVII", "YAP", "BAK", "EIF4G")
filtered_proteins <- setdiff(final_top_proteins$Protein, excluded_proteins)

# add demographic predictors
predictor_vars <- c(filtered_proteins, "age", "gender")

# subset modeling dataset
df_model <- protein_data_clean[, c("OS.time", "OS", predictor_vars)]



# ensure proper variable types
df_model$gender <- as.factor(df_model$gender)
df_model <- df_model %>% mutate(across(where(is.character), as.factor))

# map age range to midpoint of range

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


df_model$age <- age_midpoints[as.character(df_model$age)]
df_model$age <- as.numeric(df_model$age)


# remove invalid survival times
n_before <- nrow(df_model)
df_model <- df_model[!is.na(df_model$OS.time) & df_model$OS.time > 0, ]
n_after <- nrow(df_model)
cat("Removed", n_before - n_after, "rows (",
    round(100 * (n_before - n_after) / n_before, 2), "% of data)\n")

write.csv(df_model, file = "survival_dataset_cleaned.csv", row.names = FALSE)


# function to compute AICc
aicc <- function(aic, k, n) aic + (2 * k * (k + 1)) / (n - k - 1)

# initialize results table
model_results <- data.frame(
  Model = character(),
  AIC = numeric(),
  AICc = numeric(),
  BIC = numeric(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)

# build formula for modeling
cox_formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(predictor_vars, collapse = " + ")))
n <- nrow(df_model)

# --- 4.2 Cox proportional hazards model ---

cox_fit <- coxph(cox_formula, data = df_model, ties = "efron")
cox_summary <- summary(cox_fit)
cox_cindex <- cox_summary$concordance[1]
k <- length(coef(cox_fit))

model_results <- rbind(
  model_results,
  data.frame(
    Model = "Cox PH",
    AIC = AIC(cox_fit),
    AICc = aicc(AIC(cox_fit), k, n),
    BIC = BIC(cox_fit),
    C_index = cox_cindex
  )
)

# helper to compute c-index for flexsurv models
get_cindex <- function(fit, df) {
  if (inherits(fit, "flexsurvreg")) {
    f <- eval(fit$call$formula)
    X <- model.matrix(f, data = df)[, -1, drop = FALSE]
    betas <- fit$res[colnames(X), "est"]
    lp <- as.numeric(X %*% betas)
    return(concordance(Surv(OS.time, OS) ~ lp, data = df)$concordance)
  }
  if (inherits(fit, "coxph")) {
    return(summary(fit)$concordance[1])
  }
  return(NA)
}

# --- 4.3 Weibull model ---

weib_fit <- flexsurvreg(cox_formula, data = df_model, dist = "weibull")
model_results <- rbind(
  model_results,
  data.frame(
    Model = "Weibull",
    AIC = AIC(weib_fit),
    AICc = aicc(AIC(weib_fit), length(coef(weib_fit)), n),
    BIC = BIC(weib_fit),
    C_index = get_cindex(weib_fit, df_model)
  )
)

# --- 4.4 Log-logistic model ---

llog_fit <- flexsurvreg(cox_formula, data = df_model, dist = "llogis")
model_results <- rbind(
  model_results,
  data.frame(
    Model = "Log-logistic",
    AIC = AIC(llog_fit),
    AICc = aicc(AIC(llog_fit), length(coef(llog_fit)), n),
    BIC = BIC(llog_fit),
    C_index = get_cindex(llog_fit, df_model)
  )
)

# --- 4.5 Generalized gamma model ---

ggamma_fit <- flexsurvreg(cox_formula, data = df_model, dist = "gengamma")
model_results <- rbind(
  model_results,
  data.frame(
    Model = "Generalized gamma",
    AIC = AIC(ggamma_fit),
    AICc = aicc(AIC(ggamma_fit), length(coef(ggamma_fit)), n),
    BIC = BIC(ggamma_fit),
    C_index = get_cindex(ggamma_fit, df_model)
  )
)


# --- 4.6 Random Survival Forest (Apparent Harrell C-index) ---

set.seed(123)

rsf_formula <- as.formula(
  paste("Surv(OS.time, OS) ~", paste(predictor_vars, collapse = " + "))
)

rsf_fit <- rfsrc(
  formula = rsf_formula,
  data = df_model,
  ntree = 1000,
  mtry = floor(sqrt(length(predictor_vars))),
  nodesize = 15,
  splitrule = "logrank",
  importance = "permute",
  forest = TRUE
)

# ---- C-index----

rsf_pred <- predict(rsf_fit, newdata = df_model)

# use cumulative hazard at final time point as risk score
rsf_risk <- rsf_pred$chf[, ncol(rsf_pred$chf)]

# invert RSF risk
rsf_cindex_harrell_apparent <- 1 - concordance(
  Surv(OS.time, OS) ~ rsf_risk,
  data = df_model
)$concordance

# append RSF (Harrell, apparent)
model_results <- rbind(
  model_results,
  data.frame(
    Model = "Random Survival Forest (Harrell)",
    AIC = NA,
    AICc = NA,
    BIC = NA,
    C_index = rsf_cindex_harrell_apparent
  )
)



# print final results
print(model_results)


if (TRUE) {
  # ============================================================
  # Coefficient tables for fitted models
  # ============================================================
  
  model_coef_tables <- list()
  
  # ==============================
  # 1. Cox PH
  # Variable | HR (95% CI) | SE | z | p
  # ==============================
  
  cox_ci <- confint(cox_fit)
  
  cox_coef_tables <- data.frame(
    Variable = rownames(cox_summary$coef),
    `HR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      exp(cox_summary$coef[, "coef"]),
      exp(cox_ci[, 1]),
      exp(cox_ci[, 2])
    ),
    SE = cox_summary$coef[, "se(coef)"],
    z  = cox_summary$coef[, "z"],
    p  = cox_summary$coef[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
    
  )
  
  model_coef_tables[["Cox_PH"]] <- cox_coef_tables
  
  
  # ==============================
  # Helper: Weibull / Log-logistic (AFT)
  # Variable | ETR (95% CI) | SE | z | p
  # ==============================
  
  extract_aft_table <- function(fit, formula, data) {
    
    res <- fit$res
    
    # design matrix to identify covariates
    X <- model.matrix(formula, data = data)[, -1, drop = FALSE]
    covars <- colnames(X)
    
    # keep only covariate coefficients
    res <- res[rownames(res) %in% covars, , drop = FALSE]
    
    z_val <- res[, "est"] / res[, "se"]
    p_val <- 2 * (1 - pnorm(abs(z_val)))
    
    ci_low  <- res[, "est"] - 1.96 * res[, "se"]
    ci_high <- res[, "est"] + 1.96 * res[, "se"]
    
    data.frame(
      Variable = rownames(res),
      `ETR (95% CI)` = sprintf(
        "%.2f (%.2f–%.2f)",
        exp(res[, "est"]),
        exp(ci_low),
        exp(ci_high)
      ),
      SE = res[, "se"],
      z  = z_val,
      p  = p_val,
      row.names = NULL,
      check.names = FALSE
      
    )
  }
  
  model_coef_tables[["Weibull"]] <-
    extract_aft_table(weib_fit, cox_formula, df_model)
  
  model_coef_tables[["Log_logistic"]] <-
    extract_aft_table(llog_fit, cox_formula, df_model)
  
  
  # ==============================
  # Helper: Generalized gamma
  # Variable | Coef (95% CI) | SE | z | p
  # ==============================
  
  extract_gengamma_table <- function(fit, formula, data) {
    
    res <- fit$res
    
    X <- model.matrix(formula, data = data)[, -1, drop = FALSE]
    covars <- colnames(X)
    
    res <- res[rownames(res) %in% covars, , drop = FALSE]
    
    z_val <- res[, "est"] / res[, "se"]
    p_val <- 2 * (1 - pnorm(abs(z_val)))
    
    ci_low  <- res[, "est"] - 1.96 * res[, "se"]
    ci_high <- res[, "est"] + 1.96 * res[, "se"]
    
    data.frame(
      Variable = rownames(res),
      `Coef (95% CI)` = sprintf(
        "%.3f (%.3f–%.3f)",
        res[, "est"],
        ci_low,
        ci_high
      ),
      SE = res[, "se"],
      z  = z_val,
      p  = p_val,
      row.names = NULL,
      check.names = FALSE
      
    )
  }
  
  model_coef_tables[["Generalized_gamma"]] <-
    extract_gengamma_table(ggamma_fit, cox_formula, df_model)
  
  

  # Print all tables
  
  for (nm in names(model_coef_tables)) {
    cat("\n==============================\n")
    cat("Model:", nm, "\n")
    cat("==============================\n")
    print(model_coef_tables[[nm]])
  }
  
}
