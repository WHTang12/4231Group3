# import packages
library(tidyverse)

# import same dataset used in Huber(2014)
df <- read_csv("../data/causalmech.csv")

# check for NA's
any(is.na(df)) # none

# see summary statistics to spot any outliers
summary(df)

# important variables
table(df$treat) # treatment indicator D, 1 = program, 0 = control
table(df$exhealth30) # outcome Y, 1 if health is excellent 2.5 years after assignment, 0 otherwise
table(df$work2year2q) # mediator M, 1 if employed in first half of second year after assignment, 0 otherwise
table(df$female) # gender: 1 = female, 0 = male

# Data included is already filtered for observations where mediator and Y are observed

# check treatment rate by gender 
df %>% group_by(female) %>% summarise(treat_rate = mean(treat), n = n())

# check mediator rate by treatment and gender 
df %>% group_by(female, treat) %>% summarise(employ_rate = mean(work2year2q), n = n())
