library(tidyverse)
library(boot)

# import data
df <- read_csv("../data/causalmech.csv")

# split by gender (to replicate Table VI in report)
females <- df %>% filter(female == 1)
males <- df %>% filter(female == 0)

# helper function to extract Y (outcome), D (treatment), M (mediator),
# X (pre-treatment covariates), W (post-treatment, pre-mediator covariates)
extract_components <- function(df) {
  y <- df$exhealth30  # outcome Y: 1 if health is excellent 2.5 years after assignment
  d <- df$treat       # treatment indicator D: 1 = program, 0 = control
  m <- df$work2year2q # mediator M: 1 if employed in first half of second year
  
  # pre-treatment covariates X
  x <- df %>% select(schobef, trainyrbef, jobeverbef, jobyrbef, health012,
                     health0mis, pe_prb0, pe_prb0mis, everalc, alc12, everilldrugs,
                     age_cat, edumis, eduhigh, rwhite, everarr, hhsize, hhsizemis,
                     hhinc12, hhinc8, fdstamp, welf1, welf2, publicass) %>% as.matrix()
  
  # post-treatment, pre-mediator covariates W 
  w <- df %>% select(emplq4, emplq4full, pemplq4, pemplq4mis, vocq4, vocq4mis,
                     health1212, health123, pe_prb12, pe_prb12mis, narry1,
                     numkidhhf1zero, numkidhhf1onetwo, pubhse12, h_ins12a, h_ins12amis) %>% as.matrix()
  
  list(y=y, d=d, m=m, x=x, w=w)
}

female_df <- extract_components(females)
male_df   <- extract_components(males)

# Helper function to run propensity score models (probit)
propensity_models <- function(df) {
  y <- df$y 
  d <- df$d 
  m <- df$m
  x <- df$x
  w <- df$w
  
  ps_mwx <- glm(d ~ cbind(m, w, x), family = binomial("probit"))$fitted  
  ps_x <- glm(d ~ x, family = binomial("probit"))$fitted  
  ps_wx <- glm(d ~ cbind(w, x), family = binomial("probit"))$fitted  
  ps_mx <- glm(d ~ cbind(m, x), family = binomial("probit"))$fitted 
  
  list(ps_mwx=ps_mwx, ps_x=ps_x, ps_wx=ps_wx, ps_mx=ps_mx)
}

# Helper: normalised IPW mean  sum(num) / sum(den)
normalize <- function(numerator, denominator) sum(numerator) / sum(denominator)

# Compute all effects
compute_effects <- function(df) {
  y <- df$y
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  ps <- propensity_models(df)
  ps_mwx <- ps$ps_mwx  
  ps_x <- ps$ps_x 
  ps_wx <- ps$ps_wx 
  ps_mx <- ps$ps_mx  
  
  # --- Potential outcomes ---
  # E[Y(1, M(1))]
  y1m1 <- normalize(y * d/ps_x,  d/ps_x)
  # E[Y(0, M(0))]
  y0m0 <- normalize(y * (1-d) / (1-ps_x), (1-d)/(1-ps_x))
  # E[Y(1, M(0))]
  y1m0 <- normalize(y * d * (1-ps_mwx) / ((1-ps_x) * ps_mwx),
                    d * (1-ps_mwx) / ((1-ps_x) * ps_mwx))
  # E[Y(0, M(1))]
  y0m1 <- normalize(y * (1-d) * ps_mwx / (ps_x * (1-ps_mwx)),
                    (1-d) * ps_mwx / (ps_x * (1-ps_mwx)))
  
  # --- ATE ---
  ate <- mean(y[d==1]) - mean(y[d==0])
  
  # --- Average direct effects (Proposition 3 Equation 13) ---
  de_treat   <- y1m1 - y0m1
  de_control <- y1m0 - y0m0
  
  # --- Average indirect effects (Proposition 2 / equation 7) # THIS REFERS TO PRETREATMENT, SAME AS HUBER'S TABLE ---
  ie_treat   <- y1m1 - normalize(y * d * (1-ps_mx) / ps_mx,
                                 d * (1-ps_mx) / ps_mx)
  ie_control <- normalize(y * (1-d) * ps_mx / (ps_x * (1-ps_mx)),
                          (1-d) * ps_mx / (ps_x * (1-ps_mx))) - y0m0
  
  # --- Average partial indirect effects (Proposition 4) ---
  ie_partial_treat <- y1m1 -
    normalize(y * d / ps_mwx * (1-ps_mwx)/(1-ps_wx) * ps_wx/ps_x,
              d / ps_mwx * (1-ps_mwx)/(1-ps_wx) * ps_wx/ps_x)
  
  ie_partial_control <-
    normalize(y * (1-d) / (1-ps_mwx) * ps_mwx/ps_wx * (1-ps_wx)/(1-ps_x),
              (1-d) / (1-ps_mwx) * ps_mwx/ps_wx * (1-ps_wx)/(1-ps_x)) - y0m0
  
  # --- Average total indirect effects (Proposition 5) ---
  # Fit linear outcome model on treated
  lm_treat   <- lm(y[d==1] ~ cbind(m, w, x)[d==1, ])
  pred_treat  <- cbind(1, mean(m*(1-d)/(1-ps_x)), w, x) %*% lm_treat$coef
  ie_total_treat <- normalize((y - pred_treat) * d / ps_x, d / ps_x)
  
  # Fit linear outcome model on controls
  lm_control  <- lm(y[d==0] ~ cbind(m, w, x)[d==0, ])
  pred_control <- cbind(1, mean(m*d/ps_x), w, x) %*% lm_control$coef
  ie_total_control <- normalize((pred_control - y) * (1-d) / (1-ps_x),
                                (1-d) / (1-ps_x))
  
  list(
    ate = ate,
    de_treat = de_treat,
    de_control = de_control,
    ie_total_treat = ie_total_treat,
    ie_total_control = ie_total_control,
    ie_partial_treat = ie_partial_treat,
    ie_partial_control = ie_partial_control,
    ie_treat  = ie_treat,
    ie_control = ie_control
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
    unlist(compute_effects(boot_df))
  }
  
  set.seed(seed)
  boot_result <- boot(
    data = 1:length(gender_df$y),
    statistic = compute_effects_boot,
    R = R,
    gender_data = gender_df
  )
  
  estimates <- unlist(compute_effects(gender_df))
  ses <- apply(boot_result$t, 2, sd, na.rm = TRUE)
  
  results_table <- data.frame(
    Effect = names(estimates),
    Estimate = estimates,
    SE = ses,
    p_value = 2 * pnorm(-abs(estimates / ses))
  )
  
  list(results = results_table, boot_object = boot_result)
}

# Run for females
female_results <- run_bootstrap(female_df, R = 1999, seed = 123)
print(female_results$results)
# saveRDS(female_results, file = "results/original_female_results.rds")

# Run for males
male_results <- run_bootstrap(male_df, R = 1999, seed = 123)
print(male_results$results)
# saveRDS(male_results, file = "results/original_male_results.rds")


### --- Propensity score plotting to check common support --- ###
female_ps <- propensity_models(female_df)
male_ps <- propensity_models(male_df)

plot_data <- bind_rows(
  data.frame(ps = female_ps$ps_mwx,
             group = ifelse(female_df$d == 1, "Treated", "Control"),
             gender = "Females"),
  data.frame(ps = male_ps$ps_mwx,
             group = ifelse(male_df$d == 1, "Treated", "Control"),
             gender = "Males")
)

ggplot(plot_data, aes(x = ps, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ gender) +
  labs(title = "Propensity Score Overlap - P(D=1|M,W,X)",
       x = "Propensity Score", y = "Density", fill = "Group") +
  theme_minimal() +
  theme(panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold"))


### --- Trimming --- ###
trim_dataset <- function(df, trim) {
  y <- df$y
  d <- df$d
  m <- df$m
  x <- df$x
  w <- df$w
  
  ps_mwx <- glm(d ~ cbind(m, w, x), family = binomial("probit"))$fitted
  keep <- (ps_mwx > trim) & (ps_mwx < (1 - trim))
  cat("Dropped:", sum(!keep), "observations\n")
  
  list(y = y[keep], d = d[keep], m = m[keep],
       x = x[keep, , drop = FALSE],
       w = w[keep, , drop = FALSE])
}

female_df_trim <- trim_dataset(female_df, 0.05)
male_df_trim <- trim_dataset(male_df, 0.05)

female_est_trim <- compute_effects(female_df_trim)
male_est_trim <- compute_effects(male_df_trim)
female_est_trim
male_est_trim

female_trim_results <- run_bootstrap(female_df_trim, R = 1999)
male_trim_results <- run_bootstrap(male_df_trim,R = 1999)
female_trim_results <- readRDS("results/trim_female_results.rds")
male_trim_results <- readRDS("results/trim_male_results.rds")
print(female_trim_results$results)
print(male_trim_results$results)

#saveRDS(female_trim_results, file = "results/trim_female_results.rds")
#saveRDS(male_trim_results, file = "results/trim_male_results.rds")

female_ps_trim <- propensity_models(female_df_trim)
male_ps_trim <- propensity_models(male_df_trim)

plot_data_trim <- bind_rows(
  data.frame(ps = female_ps_trim$ps_mwx,
             group = ifelse(female_df_trim$d == 1, "Treated", "Control"),
             gender = "Females"),
  data.frame(ps = male_ps_trim$ps_mwx,
             group = ifelse(male_df_trim$d == 1, "Treated", "Control"),
             gender = "Males")
)

ggplot(plot_data_trim, aes(x = ps, fill = group)) +
  geom_density(alpha = 0.5, trim = T) +
  facet_wrap(~ gender) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "Propensity Score Overlap (Probit) - P(D=1|M,W,X)",
       x = "Propensity Score", y = "Density", fill = "Group") +
  theme_minimal() +
  theme(panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold"),
        strip.text   = element_text(face = "bold"))
#ggsave("results/probit_ps_overlap.png", width = 8, height = 5, dpi = 600)
