library(tidyverse)

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

female_df <- extract_components(females)
male_df <- extract_components(males)

# Helper function to run propensity score models (probit)
propensity_models <- function(df) {
  y <- df$y
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  ### Propensity scores (models) needed to obtain Proposition 1-4 (using probit)
  # P(D = 1|X)
  ps_x <- glm(d ~ x, family = binomial("probit"))$fitted
  # P(D = 1|M,X)
  ps_mx <- glm(d ~ cbind(m, x), family = binomial("probit"))$fitted
  # P(D = 1|W, X)
  ps_wx <- glm(d ~ cbind(w ,x), family = binomial("probit"))$fitted
  # P(D = 1|M,W,X)
  ps_mwx <- glm(d ~ cbind(m, w, x), family = binomial("probit"))$fitted
  
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

female_est <- compute_effects(female_df)
male_est <- compute_effects(male_df)
female_est
male_est

# Bootstrap procedure
run_bootstrap <- function(gender_df, R = 1999, seed = 123) {
  
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
  
  library(boot)
  set.seed(seed)
  
  boot_result <- boot(
    data = 1:length(gender_df$y),
    statistic = compute_effects_boot,
    R = R,
    gender_data = gender_df
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
female_results <- run_bootstrap(female_df, R = 1999, seed = 123)
print(female_results$results)
#saveRDS(female_results, file = "original_female_results.rds")

# Run for males
male_results <- run_bootstrap(male_df, R = 1999, seed = 123)  
print(male_results$results)
#saveRDS(male_results, file = "original_male_results.rds")
