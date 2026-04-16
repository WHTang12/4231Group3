#######################################################################################
# This script generates the covariance balance (love plot) for the probit and rf models
########################################################################################

library(ranger)
library(tidyverse)

# import data
df <- read_csv("../data/causalmech.csv")

# split by gender
females <- df %>% filter(female == 1)

extract_components <- function(df) {
  y <- df$exhealth30
  d <- df$treat
  m <- df$work2year2q
  
  w <- df %>% select(emplq4, emplq4full, pemplq4, pemplq4mis, vocq4, vocq4mis, 
                     health1212, health123, pe_prb12, pe_prb12mis, narry1, 
                     numkidhhf1zero, numkidhhf1onetwo, pubhse12, h_ins12a, h_ins12amis)
  
  x <- df %>% select(schobef, trainyrbef, jobeverbef, jobyrbef, health012, 
                     health0mis, pe_prb0, pe_prb0mis, everalc, alc12, everilldrugs, 
                     age_cat, edumis, eduhigh, rwhite, everarr, hhsize, hhsizemis, 
                     hhinc12, hhinc8, fdstamp, welf1, welf2, publicass)
  
  list(y=y, d=d, m=m, x=x, w=w)
}

female_data <- extract_components(females)

#############################################
# extract variables
d <- female_data$d
x <- female_data$x
m <- female_data$m
w <- female_data$w

# put a prefix infront of covariates so easier to identify
colnames(x) <- paste0("x_", colnames(x))
colnames(w) <- paste0("w_", colnames(w))
covs <- data.frame(x, w, m = m)

ps_df <- data.frame(d = d, m = m, w, x)

# probit propensity model for P(D=1|M,W,X)
probit <- glm(d ~ ., data = ps_df, family = binomial("probit"))$fitted

###--- rf propensity model for P(D=1|M,W,X) (post-trimming) ---###
ps_df_rf <- data.frame(d = factor(d), m = m, w, x)

rf <- ranger(
  d ~ .,
  data = ps_df_rf,
  probability = TRUE,
  mtry = 20,
  min.node.size = 20,
  sample.fraction = 0.5,
  num.trees = 500,
  seed = 4231
)

rf_ps <- predict(rf, data = ps_df_rf)$predictions[, "1"]

# trimming
keep <- (rf_ps > 0.05) & (rf_ps < 0.95)

d_trim <- d[keep]
covs_trim <- covs[keep, ]
probit_trim <- probit[keep]
rf_ps_trim <- rf_ps[keep]

### --- Compute abs standardised mean difference --- ###
# helper function to compute (abs) standardized mean difference, exactly same as in previous scripts
compute_smd <- function(d, ps, covariates) {
  # which weights to use
  # Over here, i decide to use different weights to compute SMD depending on which propensity score model is estimated
  w <- ifelse(d == 1, 1/ps, 1/(1 - ps)) 
  
  # formula
  smd_vec <- apply(covariates, 2, function(x) {
    w_mean_t <- sum(w[d==1] * x[d==1]) / sum(w[d==1])
    w_mean_c <- sum(w[d==0] * x[d==0]) / sum(w[d==0])
    pooled_sd <- sqrt((var(x[d==1]) + var(x[d==0])) / 2)
    abs((w_mean_t - w_mean_c) / pooled_sd)
  })
  
  smd_vec
}

smd_probit <- compute_smd(d_trim, probit_trim, covs_trim)
smd_rf <- compute_smd(d_trim, rf_ps_trim, covs_trim)

### --- Plotting --- ###
plot_df <- tibble(
  variable = colnames(covs_trim),
  probit = smd_probit,
  rf = smd_rf
) %>%
  mutate(order_val = rf) %>%
  pivot_longer(cols = c(probit, rf),
               names_to = "model",
               values_to = "smd")

# this is just to swap the colors of the models, not important
plot_df <- plot_df %>%
  mutate(model = factor(model, levels = c("rf", "probit")))

# plot
ggplot(plot_df, aes(x = smd, y = reorder(variable, order_val), color = model)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey60") +
  scale_color_discrete(
    breaks = c("rf", "probit"),
    labels = c("rf" = "Random Forest", "probit" = "Probit")
  ) +
  labs(
    title = "Covariate Balance (Females): Probit vs Random Forest",
    x = "Absolute Standardized Mean Difference",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )
#ggsave("results/love_plot_females.png", width = 8, height = 10, dpi = 600)

### --- Males --- ###
males <- df %>% filter(female == 0)
male_data <- extract_components(males)
d <- male_data$d
x <- male_data$x
m <- male_data$m
w <- male_data$w

colnames(x) <- paste0("x_", colnames(x))
colnames(w) <- paste0("w_", colnames(w))
covs <- data.frame(x, w, m = m)

ps_df <- data.frame(d = d, m = m, w, x)
probit <- glm(d ~ ., data = ps_df, family = binomial("probit"))$fitted

ps_df_rf <- data.frame(d = factor(d), m = m, w, x)
rf <- ranger(
  d ~ .,
  data = ps_df_rf,
  probability = TRUE,
  mtry = 20,
  min.node.size = 20,
  sample.fraction = 0.5,
  num.trees = 500,
  seed = 4231
)

rf_ps <- predict(rf, data = ps_df_rf)$predictions[, "1"]

keep <- (rf_ps > 0.05) & (rf_ps < 0.95)

d_trim <- d[keep]
covs_trim <- covs[keep, ]
probit_trim <- probit[keep]
rf_ps_trim <- rf_ps[keep]

smd_probit <- compute_smd(d_trim, probit_trim, covs_trim)
smd_rf <- compute_smd(d_trim, rf_ps_trim, covs_trim)

plot_df <- tibble(
  variable = colnames(covs_trim),
  probit = smd_probit,
  rf = smd_rf
) %>%
  mutate(order_val = rf) %>%
  pivot_longer(cols = c(probit, rf),
               names_to = "model",
               values_to = "smd")

plot_df <- plot_df %>%
  mutate(model = factor(model, levels = c("rf", "probit")))

ggplot(plot_df, aes(x = smd, y = reorder(variable, order_val), color = model)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey60") +
  scale_color_discrete(
    breaks = c("rf", "probit"),
    labels = c("rf" = "Random Forest", "probit" = "Probit")
  ) +
  labs(
    title = "Covariate Balance (Males): Probit vs Random Forest",
    x = "Absolute Standardized Mean Difference",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )
#ggsave("results/love_plot_males.png", width = 8, height = 10, dpi = 600)
