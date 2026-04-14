library(tidyverse)
library(boot)
library(causalweight)

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

dml_f <- medDML(y = female_data$y,
                d = female_data$d, 
                m = female_data$m,
                x = female_data$x,
                trim = 0.05,
                normalized = T, # normalized weights for ipw weights, same as original ipw model
                MLmethod = "ensemble") # ensemble of lasso, rf, xgboost, svm

dml_f

dml_m <- medDML(y = male_data$y,
                d = male_data$d, 
                m = male_data$m,
                x = male_data$x,
                trim = 0.05,
                normalized = T, # normalized weights for ipw weights, same as original ipw model
                MLmethod = "ensemble") # ensemble of lasso, rf, xgboost, svm

dml_m
  
