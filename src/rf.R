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
compute_max_asmd <- function(d, ps_x, covariates) {
  # which weights to use
  # Over here, i decide to use different weights to compute ASMD depending on which propensity score model is estimated
  w <- ifelse(d == 1, 1/ps_x, 1/(1 - ps_x)) 
  
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

# tuning the rf model
# only tuning on mtry, because of reference paper that mentions mtry as being most significant determinant of covariate balance
tune_ps <- function(df, type, label = "") {
  d <- df$d
  x <- df$x
  m <- df$m
  w <- df$w
  
  # see which propensity score model is being tuned
  features <- switch(type,
                     ps_x = x,
                     ps_mx = cbind(m,x),
                     ps_wx = cbind(w,x),
                     ps_mwx = cbind(m,w,x)
  )
  
  p <- ncol(features) # number of features
  sqrt_p <- round(sqrt(p)) # sqrt(number of features)
  
  # parameter grid: try out [sqrt(p)+/-1, p/3 and p/2]
  mtry_grid <- unique(pmax(1, c(sqrt_p - 1, sqrt_p, sqrt_p + 1,
                                round(p / 3), round(p / 2))))
  mtry_grid <- sort(mtry_grid[mtry_grid <= p])
  
  # Using in-sample, since we are ensuring covariate balance for our propensity model
  # Not sure whether using in-sample is ideal, but again, according to reference paper, most papers use in-sample
  train_df    <- as.data.frame(features)
  train_df$.d <- factor(d)
  
  results <- purrr::map_dfr(mtry_grid, function(m) {
    rf <- ranger(
      formula       = .d ~ .,
      data          = train_df,
      mtry          = m,
      num.trees     = 500,
      probability   = TRUE,
      min.node.size = 10,
      seed          = 4231
    )
    
    ps_hat <- pmax(pmin(rf$predictions[, "1"], 0.99), 0.01)
    
    asmd <- compute_max_asmd(d, ps_hat, features)
    
    tibble::tibble(mtry = m, max_asmd = asmd)
  })
  
  print(results)
  
  best_mtry <- results$mtry[which.min(results$max_asmd)]
  cat("best mtry:", best_mtry, "\n\n")
  
  best_mtry
}

best_ps_x_female   <- tune_ps(female_data, "ps_x", "FEMALES")
best_ps_mx_female  <- tune_ps(female_data, "ps_mx", "FEMALES")
best_ps_wx_female  <- tune_ps(female_data, "ps_wx", "FEMALES")
best_ps_mwx_female <- tune_ps(female_data, "ps_mwx", "FEMALES")

best_ps_x_male  <- tune_ps(male_data, "ps_x", "MALES")
best_ps_mx_male  <- tune_ps(male_data, "ps_mx", "MALES")
best_ps_wx_male  <- tune_ps(male_data, "ps_wx", "MALES")
best_ps_mwx_male <- tune_ps(male_data, "ps_mwx", "MALES")

### Same as before, just replacing with RF for propensity score model ###
# Helper function to run propensity score models (RF)
propensity_models <- function(df, mtry_list = list(ps_x=4, ps_mx=4, ps_wx=7, ps_mwx=6)) {
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  # helper to fit RF
  fit_rf <- function(features, mtry_val) {
    data_rf <- data.frame(d = as.factor(d), features)
    
    model <- ranger(
      d ~ .,
      data = data_rf,
      probability = TRUE,
      mtry = mtry_val,
      num.trees = 500,
      seed = 4231
    )
    
    preds <- predict(model, data_rf)$predictions[, "1"]
    return(preds)
  }
  
  # PS models with custom mtry
  ps_x   <- fit_rf(as.data.frame(x),           mtry_list$ps_x)
  ps_mx  <- fit_rf(data.frame(m = m, x),       mtry_list$ps_mx)
  ps_wx  <- fit_rf(data.frame(w, x),           mtry_list$ps_wx)
  ps_mwx <- fit_rf(data.frame(m = m, w, x),    mtry_list$ps_mwx)
  
  list(ps_x = ps_x, ps_mx = ps_mx, ps_wx = ps_wx, ps_mwx = ps_mwx)
}

# Helper function to replicate Table VI (computation of each effects)
compute_effects <- function(df) {
  y <- df$y
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  # propensity scores
  ps <- propensity_models(df)
  ps_x <- ps$ps_x
  ps_mx <- ps$ps_mx
  ps_wx <- ps$ps_wx
  ps_mwx <- ps$ps_mwx
  
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
    ate                = ate,
    de_treat           = de_treat,
    de_control         = de_control,
    ie_total_treat     = ie_total_treat,
    ie_total_control   = ie_total_control,
    ie_partial_treat   = ie_partial_treat,
    ie_partial_control = ie_partial_control,
    ie_treat           = ie_treat,
    ie_control         = ie_control
  )
}

female_rf_est <- compute_effects(female_data)
male_rf_est <- compute_effects(male_data)
female_rf_est
male_rf_est

### --- Bootstrap for standard errors --- ###
# Same procedure as in original, just using parallel computation to speed up
# With 5 cores, it still takes an hour for each gender
# Results are saved in rf_female_results.rds and rf_male_results.rds

# Use parallelization
library(parallel)
cl <- makeCluster(5)
clusterExport(cl, varlist = c("compute_effects", "propensity_models"))
clusterEvalQ(cl, library(ranger))

run_bootstrap <- function(gender_df, R = 1999, seed = 4231) {
  
  compute_effects_boot <- function(data, indices, gender_data) {
    boot_df <- list(
      y = gender_data$y[indices],
      d = gender_data$d[indices],
      m = gender_data$m[indices],
      x = gender_data$x[indices, , drop = FALSE],
      w = gender_data$w[indices, , drop = FALSE]
    )
    est <- compute_effects(boot_df)
    return(unlist(est))
  }
  set.seed(seed)
  
  boot_result <- boot(
    data = 1:length(gender_df$y),
    statistic = compute_effects_boot,
    R = R,
    gender_data = gender_df,
    parallel = "snow",
    ncpus = 5,
    cl = cl
  )
  
  # Results table
  results_table <- data.frame(
    Effect = names(compute_effects(gender_df)),
    Estimate = boot_result$t0,
    SE = apply(boot_result$t, 2, sd, na.rm = TRUE),
    p_value = 2 * pnorm(-abs(boot_result$t0 / apply(boot_result$t, 2, sd, na.rm = TRUE)))
  )
  
  return(list(results = results_table, boot_object = boot_result))
}

# Run for females
female_rf_results <- run_bootstrap(female_data, R = 1999, seed = 4231)
#female_rf_results <- readRDS("rf_female_results.rds")
print(female_rf_results$results)
#saveRDS(female_rf_results, file = "rf_female_results.rds")

# Run for males
male_rf_results <- run_bootstrap(male_data, R = 1999, seed = 4231)
#male_rf_results <- readRDS("rf_male_results.rds")
print(male_rf_results$results)
#saveRDS(male_rf_results, file = "rf_male_results.rds")

stopCluster(cl)
