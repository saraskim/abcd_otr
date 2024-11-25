# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Nuisance model code to estimate pseudo-outcomes for DAY 3 DIARRHEA outcome
#
#
# Sara Kim
# 28 January 2024
#
# ------------------------------------------------------------------------------------------

# ----- Make sure "0. wrappers.d3diar.R" script has been run

# ----- Load libraries
library(SuperLearner)
library(dplyr)
library(drotr)
library(ranger)
library(earth)
library(glmnet)
library(stats)
library(xgboost)


# ----- Load dataset 
# Deidentified dataset available upon reasonable request with data use agreement 
abcd_data <- read.csv("abcd_otr_deidentified.csv")[,-1] 
# Ensure categorical variables are factor format
abcd_data$site <- as.factor(abcd_data$site)
abcd_data$dy1_scrn_dehydr <- as.factor(abcd_data$dy1_scrn_dehydr)
abcd_data$dy1_ant_sex <- as.factor(abcd_data$dy1_ant_sex)
abcd_data$an_ses_quintile <- as.factor(abcd_data$an_ses_quintile)
abcd_data$dy1_scrn_vomitall <- as.factor(abcd_data$dy1_scrn_vomitall)
abcd_data$month_en <- as.factor(abcd_data$month_en)



# ----- Run nuisance models over five seeds to get pseudo-outcomes

# set seed
set.seed(NULL)
set.seed(150)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance150_day3diar_30apr2024.RData")

# set seed
set.seed(NULL)
set.seed(150)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance150_day3diar_30apr2024.RData")


# set seed
set.seed(NULL)
set.seed(150)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance150_day3diar_30apr2024.RData")

# set seed
set.seed(NULL)
set.seed(300)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance300_day3diar_30apr2024.RData")

# set seed
set.seed(NULL)
set.seed(450)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance450_day3diar_30apr2024.RData")

# set seed
set.seed(NULL)
set.seed(600)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance600_day3diar_30apr2024.RData")

# set seed
set.seed(NULL)
set.seed(750)

# Run nuisance model 
nuisance_output <- learn_nuisance(df = abcd_data,
                                  Y_name = "day3diar",
                                  A_name = "an_grp_01",
                                  W_list = W_list,
                                  id_name = "pid", 
                                  sl.library.outcome = sl.library.outcome,
                                  sl.library.treatment = sl.library.treatment,
                                  sl.library.missingness = sl.library.missingness,
                                  outcome_type = "binomial",
                                  k_folds = 10,
                                  ps_trunc_level = 0.01)

nuisance_models <- nuisance_output$nuisance_models

# Save nuisance output
saveRDS(nuisance_output, file="nuisance750_day3diar_30apr2024.RData")


