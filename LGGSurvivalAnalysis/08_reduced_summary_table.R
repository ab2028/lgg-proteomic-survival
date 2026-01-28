# 8. Reduced

# Signs only - for reference
if (FALSE) {
make_reference_row <- function(df, vars) {
  ref <- list()
  
  for (v in vars) {
    if (is.factor(df[[v]])) {
      ref[[v]] <- levels(df[[v]])[1]
    } else {
      ref[[v]] <- median(df[[v]], na.rm = TRUE)
    }
  }
  
  as.data.frame(ref)
}
shift_variable <- function(df_ref, df, var) {
  df_new <- df_ref
  
  if (is.factor(df[[var]])) {
    lv <- levels(df[[var]])
    if (length(lv) < 2) stop("Factor ", var, " has <2 levels")
    df_new[[var]] <- lv[2]
  } else {
    sd_val <- sd(df[[var]], na.rm = TRUE)
    df_new[[var]] <- df_ref[[var]] + sd_val
  }
  
  df_new
}


median_survival_cox <- function(fit, newdata) {
  sf <- survfit(fit, newdata = newdata)
  summary(sf)$table["median"]
}

median_survival_flex <- function(fit, newdata) {
  summary(
    fit,
    newdata = newdata,
    type = "quantile",
    quantiles = 0.5
  )[[1]]$est
}

check_direction <- function(fit, model_type, df, var, vars_in_model) {
  
  ref <- make_reference_row(df, vars_in_model)
  shifted <- shift_variable(ref, df, var)
  
  m0 <- if (model_type == "cox") {
    median_survival_cox(fit, ref)
  } else {
    median_survival_flex(fit, ref)
  }
  
  m1 <- if (model_type == "cox") {
    median_survival_cox(fit, shifted)
  } else {
    median_survival_flex(fit, shifted)
  }
  
  if (is.na(m0) || is.na(m1)) return(NA)
  
  if (m1 > m0) {
    "increase → longer survival"
  } else {
    "increase → shorter survival"
  }
}


cox_directions <- sapply(
  result_cox$final_vars,
  check_direction,
  fit = result_cox$final_fit,
  model_type = "cox",
  df = df_model,
  vars_in_model = result_cox$final_vars
)

cox_directions

weibull_directions <- sapply(
  result_weibull$final_vars,
  check_direction,
  fit = result_weibull$final_fit,
  model_type = "weibull",
  df = df_model,
  vars_in_model = result_weibull$final_vars
)

weibull_directions

llogis_directions <- sapply(
  result_llogis$final_vars,
  check_direction,
  fit = result_llogis$final_fit,
  model_type = "llogis",
  df = df_model,
  vars_in_model = result_llogis$final_vars
)

llogis_directions

gengamma_directions <- sapply(
  result_gengamma$final_vars,
  check_direction,
  fit = result_gengamma$final_fit,
  model_type = "gengamma",
  df = df_model,
  vars_in_model = result_gengamma$final_vars
)

gengamma_directions

dir_to_symbol <- function(x) {
  if (is.na(x)) return("X")
  if (x == "increase → longer survival") return("+")
  if (x == "increase → shorter survival") return("-")
  "X"
}

cox_df <- data.frame(
  Variable = names(cox_directions),
  Cox = sapply(cox_directions, dir_to_symbol),
  stringsAsFactors = FALSE
)

weibull_df <- data.frame(
  Variable = names(weibull_directions),
  Weibull = sapply(weibull_directions, dir_to_symbol),
  stringsAsFactors = FALSE
)

llogis_df <- data.frame(
  Variable = names(llogis_directions),
  LogLogistic = sapply(llogis_directions, dir_to_symbol),
  stringsAsFactors = FALSE
)

gengamma_df <- data.frame(
  Variable = names(gengamma_directions),
  GenGamma = sapply(gengamma_directions, dir_to_symbol),
  stringsAsFactors = FALSE
)

library(dplyr)
library(purrr)

direction_table <- list(
  cox_df,
  weibull_df,
  llogis_df,
  gengamma_df
) %>%
  reduce(full_join, by = "Variable") %>%
  mutate(across(-Variable, ~ ifelse(is.na(.x), "X", .x))) %>%
  arrange(Variable)

# FOR DIAGNOSTICS ONLY
print(direction_table)
}



# Summary Table

if (TRUE) 
{
  p_to_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else ""
  }
  
  fmt_est <- function(est, lo, hi, p) {
    stars <- p_to_stars(p)
    sprintf("%.2f (%.2f–%.2f)%s", est, lo, hi, stars)
  }
  extract_cox_col <- function(result) {
    fit <- result$final_fit
    sm  <- summary(fit)
    ci  <- confint(fit)
    
    df_out <- data.frame(
      Variable = rownames(sm$coef),
      Cox = mapply(
        fmt_est,
        exp(sm$coef[, "coef"]),
        exp(ci[, 1]),
        exp(ci[, 2]),
        sm$coef[, "Pr(>|z|)"]
      ),
      stringsAsFactors = FALSE
    )
    
    colnames(df_out)[2] <- "Cox (HR)"
      df_out
    
  }
  extract_aft_col <- function(result, drop_params, colname) {
    res <- result$final_fit$res
    res <- res[!(rownames(res) %in% drop_params), , drop = FALSE]
    
    z  <- res[, "est"] / res[, "se"]
    p  <- 2 * (1 - pnorm(abs(z)))
    lo <- res[, "est"] - 1.96 * res[, "se"]
    hi <- res[, "est"] + 1.96 * res[, "se"]
    
    data.frame(
      Variable = rownames(res),
      tmp = mapply(
        fmt_est,
        exp(res[, "est"]),
        exp(lo),
        exp(hi),
        p
      ),
      stringsAsFactors = FALSE
    ) |>
      setNames(c("Variable", colname))
  }
  
  extract_gengamma_col <- function(result) {
    res <- result$final_fit$res
    res <- res[!(rownames(res) %in% c("shape", "sigma", "Q", "mu")), , drop = FALSE]
    
    z  <- res[, "est"] / res[, "se"]
    p  <- 2 * (1 - pnorm(abs(z)))
    lo <- res[, "est"] - 1.96 * res[, "se"]
    hi <- res[, "est"] + 1.96 * res[, "se"]
    
    data.frame(
      Variable = rownames(res),
      "GenGamma (coef)" = mapply(
        fmt_est,
        res[, "est"],
        lo,
        hi,
        p
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  library(dplyr)
  library(purrr)
  
  cox_col <- extract_cox_col(result_cox)
  
  weibull_col <- extract_aft_col(
    result_weibull,
    drop_params = c("shape", "scale"),
    colname = "Weibull (ETR)"
  )
  
  llogis_col <- extract_aft_col(
    result_llogis,
    drop_params = c("shape", "scale"),
    colname = "LogLogistic (ETR)"
  )
  
  gengamma_col <- extract_gengamma_col(result_gengamma)
  
  final_table <- list(
    cox_col,
    weibull_col,
    llogis_col,
    gengamma_col
  ) |>
    reduce(full_join, by = "Variable") |>
    mutate(across(-Variable, ~ ifelse(is.na(.x), "---", .x))) |>
    arrange(Variable)
  
  print(final_table)
  
}


