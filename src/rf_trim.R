library(tidyverse)
library(boot)
library(ranger)

# import data
df <- read_csv("../data/causalmech.csv")

# split by gender (to replicate Table VI in report)
females <- df %>% filter(female == 1)
males <- df %>% filter(female == 0)

# helper function to extract Y (outcome), D (treatment), M (mediator), X(pre-treatment covariates), W(post-treatment, pre-mediator covariates)
extract_components <- function(df) {
  y <- df$exhealth30 # outcome Y, 1 if health is excellent 2.5 years after assignment, 0 otherwise
  d <- df$treat # treatment indicator D, 1 = program, 0 = control
  m <- df$work2year2q # mediator M, 1 if employed in first half of second year after assignment, 0 otherwise
  
  # pre-treatment covariates X
  x <- df %>% select(emplq4, emplq4full, pemplq4, pemplq4mis, vocq4, vocq4mis, 
                     health1212, health123,  pe_prb12, pe_prb12mis,  narry1, 
                     numkidhhf1zero, numkidhhf1onetwo, pubhse12, h_ins12a, h_ins12amis) %>% as.matrix()
  
  # post-treatment, pre-mediator covariates W
  w <- df %>% select(schobef, trainyrbef, jobeverbef, jobyrbef, health012, 
                     health0mis, pe_prb0, pe_prb0mis, everalc, alc12, everilldrugs, 
                     age_cat, edumis, eduhigh, rwhite, everarr, hhsize, hhsizemis, 
                     hhinc12, hhinc8, fdstamp, welf1, welf2, publicass) %>% as.matrix()
  
  list(y=y, d=d, m=m, x=x, w=w)
}

female_data <- extract_components(females)
male_data <- extract_components(males)

# Helper function to compute Absolute Standardized Mean Difference
compute_max_asmd <- function(d, ps, covariates) {
  # which weights to use
  # Over here, i decide to use different weights to compute ASMD depending on which propensity score model is estimated
  w <- ifelse(d == 1, 1/ps, 1/(1 - ps)) 
  
  # formula
  asmd_vec <- apply(covariates, 2, function(x) {
    w_mean_t  <- sum(w[d==1] * x[d==1]) / sum(w[d==1])
    w_mean_c  <- sum(w[d==0] * x[d==0]) / sum(w[d==0])
    pooled_sd <- sqrt((var(x[d==1]) + var(x[d==0])) / 2)
    if (pooled_sd == 0) return(0)
    abs(w_mean_t - w_mean_c) / pooled_sd
  })
  
  # take the maximum value, since we are minimizing the maximum ASMD
  max(asmd_vec, na.rm = TRUE)
}

# Helper function to fit RF
fit_rf <- function(d, features, params) {
  data_rf <- data.frame(d = as.factor(d), features)
  
  model <- ranger(
    d ~ .,
    data = data_rf,
    probability = TRUE,
    mtry = params$mtry,
    min.node.size = params$node_size,
    sample.fraction = params$samp_frac,
    num.trees = 500,
    seed = 4231
  )
  
  predict(model, data_rf)$predictions[, "1"]
}

# tuning the rf model
tune_ps <- function(df, type, trim = 0.05) {
  d <- df$d
  x <- df$x
  m <- df$m
  w <- df$w
  
  # select features
  features <- switch(type,
                     ps_x = x,
                     ps_mx = cbind(m,x),
                     ps_wx = cbind(w,x),
                     ps_mwx = cbind(m,w,x)
  )
  
  p <- ncol(features)
  sqrt_p <- round(sqrt(p))
  
  # grids
  mtry_grid <- unique(pmax(1, c(sqrt_p - 1, sqrt_p, sqrt_p + 1,
                                round(p / 3), round(p / 2))))
  mtry_grid <- sort(mtry_grid[mtry_grid <= p])
  
  node_size_grid <- c(5, 10, 20)
  samp_frac_grid <- c(0.5, 0.7, 1)
  
  train_df <- as.data.frame(features)
  train_df$.d <- factor(d)
  
  results <- expand.grid(
    mtry = mtry_grid,
    node_size = node_size_grid,
    samp_frac = samp_frac_grid
  ) %>%
    pmap_dfr(function(mtry, node_size, samp_frac) {
      
      params <- list(mtry = mtry, node_size = node_size, samp_frac = samp_frac)
      ps_hat <- fit_rf(d, train_df %>% select(-.d), params)
      
      # trim = 0.05
      keep <- (ps_hat > trim) & (ps_hat < (1-trim))
      d_trim <- d[keep]
      features_trim <- features[keep, , drop = FALSE]
      ps_trim <- ps_hat[keep]
      
      # compute ASMD
      asmd <- compute_max_asmd(d_trim, ps_trim, features_trim)
      
      tibble(
        mtry = mtry,
        node_size = node_size,
        samp_frac = samp_frac,
        max_asmd = asmd,
        n_kept = sum(keep)
      )
    })
  
  print(results)
  
  # pick parameters that minimizes the maximum asmd
  best_row <- results[which.min(results$max_asmd), ]
  cat("Best params:\n")
  print(best_row)
  cat("\n")
  
  best_row
}

best_ps_x_female <- tune_ps(female_data, "ps_x", trim = 0.05)
best_ps_mx_female <- tune_ps(female_data, "ps_mx", trim = 0.05)
best_ps_wx_female <- tune_ps(female_data, "ps_wx", trim = 0.05)
best_ps_mwx_female <- tune_ps(female_data, "ps_mwx", trim = 0.05)

params_list_female <- list(
  ps_x = as.list(best_ps_x_female),
  ps_mx = as.list(best_ps_mx_female),
  ps_wx = as.list(best_ps_wx_female),
  ps_mwx = as.list(best_ps_mwx_female)
)

best_ps_x_male <- tune_ps(male_data, "ps_x", trim = 0.05)
best_ps_mx_male <- tune_ps(male_data, "ps_mx", trim = 0.05)
best_ps_wx_male <- tune_ps(male_data, "ps_wx", trim = 0.05)
best_ps_mwx_male <- tune_ps(male_data, "ps_mwx", trim = 0.05)

params_list_male <- list(
  ps_x = as.list(best_ps_x_male),
  ps_mx = as.list(best_ps_mx_male),
  ps_wx = as.list(best_ps_wx_male),
  ps_mwx = as.list(best_ps_mwx_male)
)

### Same as before, just replacing with RF for propensity score model ###
# Helper function to run propensity score models (RF)
propensity_models <- function(df, params_list, trim = 0.05) {
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  # PS models with custom mtry
  ps_x   <- fit_rf(d, as.data.frame(x), params_list$ps_x)
  ps_mx  <- fit_rf(d, data.frame(m = m, x), params_list$ps_mx)
  ps_wx  <- fit_rf(d, data.frame(w, x), params_list$ps_wx)
  ps_mwx <- fit_rf(d, data.frame(m = m, w, x), params_list$ps_mwx)
  
  # following original trimming in the paper, we only trim on ps_mwx
  keep <- (ps_mwx > trim) & (ps_mwx < (1 - trim))
  
  list(ps_x = ps_x[keep], ps_mx = ps_mx[keep], ps_wx = ps_wx[keep], ps_mwx = ps_mwx[keep], keep=keep)
}

# Helper function to replicate Table VI (computation of each effects)
compute_effects <- function(df, params_list) {
  y <- df$y
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  # propensity scores
  ps <- propensity_models(df, params_list)
  ps_x <- ps$ps_x
  ps_mx <- ps$ps_mx
  ps_wx <- ps$ps_wx
  ps_mwx <- ps$ps_mwx
  keep <- ps$keep
  
  # trim
  y <- y[keep]
  d <- d[keep]
  m <- m[keep]
  x <- x[keep, , drop = FALSE]
  w <- w[keep, , drop = FALSE]
  
  # helper function to normalize weights (see Section 3 of paper/Lecture 3)
  normalize <- function(numerator, denominator) {
    sum(numerator)/sum(denominator)
  }
  
  # ate (difference in means of outcome between treated and untreated)
  ate <- mean(y[d==1]) - mean(y[d==0])
  
  ### proposition 3: average direct effect (equation 13) ###
  de_treat <- normalize(y * d/ps_x, d/ps_x) -
    normalize(y*(1-d)*ps_mwx/(ps_x*(1-ps_mwx)),
              (1-d)*ps_mwx/(ps_x*(1-ps_mwx)))
  
  de_control <- normalize(y * d * (1-ps_mwx)/((1-ps_x)*ps_mwx),
                          d * (1-ps_mwx)/((1-ps_x)*ps_mwx)) -
    normalize(y*(1-d)/(1-ps_x), (1-d)/(1-ps_x))
  
  # proposition 2: average indirect effect (equation 7)
  ie_treat <- normalize(y*d/ps_x, d/ps_x) - normalize(y*d*(1-ps_mx)/ps_mx, d*(1-ps_mx)/ps_mx)
  ie_control <- normalize(y*(1-d)*ps_mx/(ps_x*(1-ps_mx)), (1-d)*ps_mx/(ps_x*(1-ps_mx))) - normalize(y*(1-d)/(1-ps_x), (1-d)/(1-ps_x))
  
  # preposition 4: average partial indirect effects
  ie_partial_treat <- normalize(y * d/ps_x, d/ps_x) - normalize(y * d / ps_mwx * (1-ps_mwx)/(1-ps_wx) * ps_wx/ps_x, 
                                                                d / ps_mwx * (1-ps_mwx)/(1-ps_wx) * ps_wx/ps_x)
  
  ie_partial_control <- normalize(y * (1-d) / (1-ps_mwx) * ps_mwx/ps_wx * (1-ps_wx)/(1-ps_x),
                                  (1-d) / (1-ps_mwx) * ps_mwx/ps_wx * (1-ps_wx)/(1-ps_x)) - normalize(y * (1-d)/(1-ps_x), (1-d)/(1-ps_x))
  
  # preposition 5: average total indirect effect
  # fit (linear) outcome model on treated observations
  lm_treat <- lm(y[d==1] ~ cbind(m, w, x)[d==1,])
  # predict
  pred_treat <- cbind(1, mean(m*(1-d)/(1-ps_x)), w, x) %*% lm_treat$coef
  # plug into formula
  ie_total_treat <- normalize((y - pred_treat) * d/ps_x, d/ps_x)
  
  
  # fit (linear) outcome model on untreated observations
  lm_control <- lm(y[d==0] ~ cbind(m, w, x)[d==0,])
  # predict
  pred_control <- cbind(1, mean(m*d/ps_x), w, x) %*% lm_control$coef
  # plug into formula
  ie_total_control <- normalize((pred_control - y) * (1-d)/(1-ps_x), (1-d)/(1-ps_x))
  
  
  list(
    ate = ate,
    de_treat = de_treat,
    de_control = de_control,
    ie_total_treat = ie_total_treat,
    ie_total_control = ie_total_control,
    ie_partial_treat = ie_partial_treat,
    ie_partial_control = ie_partial_control,
    ie_treat = ie_treat,
    ie_control = ie_control
  )
}

female_rf_trim_est <- compute_effects(female_data, params_list = params_list_female)
male_rf_trim_est <- compute_effects(male_data, params_list = params_list_male)
female_rf_trim_est
male_rf_trim_est

### --- Bootstrap for standard errors --- ###
# Same procedure as in original, just using parallel computation to speed up
# With 5 cores, it still takes an hour for each gender
# Results are saved in rf_female_trim_results.rds and rf_male_trim_results.rds

# Use parallelization
library(parallel)
cl <- makeCluster(5)
clusterExport(cl, varlist = c("compute_effects", "propensity_models", "fit_rf"))
clusterEvalQ(cl, library(ranger))

run_bootstrap <- function(gender_df, params_list, R = 1999, seed = 4231) {
  
  compute_effects_boot <- function(data, indices, gender_data, params_list) {
    boot_df <- list(
      y = gender_data$y[indices],
      d = gender_data$d[indices],
      m = gender_data$m[indices],
      x = gender_data$x[indices, , drop = FALSE],
      w = gender_data$w[indices, , drop = FALSE]
    )
    est <- compute_effects(boot_df, params_list)
    return(unlist(est))
  }
  set.seed(seed)
  
  boot_result <- boot(
    data = 1:length(gender_df$y),
    statistic = compute_effects_boot,
    R = R,
    gender_data = gender_df,
    params_list = params_list,
    parallel = "snow",
    ncpus = 5,
    cl = cl
  )
  
  # Results table
  results_table <- data.frame(
    Effect = names(compute_effects(gender_df, params_list)),
    Estimate = boot_result$t0,
    SE = apply(boot_result$t, 2, sd, na.rm = TRUE),
    p_value = 2 * pnorm(-abs(boot_result$t0 / apply(boot_result$t, 2, sd, na.rm = TRUE)))
  )
  
  return(list(results = results_table, boot_object = boot_result))
}

# Run for females
female_rf_trim_results <- run_bootstrap(female_data, params_list_female, R = 1999, seed = 4231)
female_rf_trim_results <- readRDS("results/rf_female_trim_results.rds")
print(female_rf_trim_results$results)
#saveRDS(female_rf_trim_results, file = "results/rf_female_trim_results.rds")

# Run for males
male_rf_trim_results <- run_bootstrap(male_data, params_list_male, R = 1999, seed = 4231)
male_rf_trim_results <- readRDS("results/male_rf_trim_results.rds")
print(male_rf_trim_results$results)
#saveRDS(male_rf_trim_results, file = "results/male_rf_trim_results.rds")

stopCluster(cl)

########################################################
