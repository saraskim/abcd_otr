# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Treatment effects for the Clinical Symptoms code for linear growth outcome
#
#
# Sara Kim
# 28 January 2024
#
# ------------------------------------------------------------------------------------------

# Please ensure "0. wrappers.laz.R" script has been run and outputs from "1. nuisance.laz.R" are in working directory

# Load libraries
library(SuperLearner)
library(dplyr)
library(drotr)
library(ranger)
library(earth)
library(glmnet)
library(stats)
library(xgboost)


# Load dataset 
abcd_data <- read.csv("abcd_otr_deidentified.csv")[,-1] 
abcd_data$site <- as.factor(abcd_data$site)
abcd_data$dy1_scrn_dehydr <- as.factor(abcd_data$dy1_scrn_dehydr)
abcd_data$dy1_ant_sex <- as.factor(abcd_data$dy1_ant_sex)
abcd_data$an_ses_quintile <- as.factor(abcd_data$an_ses_quintile)
abcd_data$dy1_scrn_vomitall <- as.factor(abcd_data$dy1_scrn_vomitall)
abcd_data$month_en <- as.factor(abcd_data$month_en)


# ----- Define CATE models to estimate child level expected benefit for the Clinical Symptoms rule


# clinical symptoms
SL.cate.clin <- function(Y, X, newX, family, ...){
  sl.cate.clin_fit <- glm(Y ~ dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + 
                            dy1_scrn_diardays + dy1_scrn_dehydr, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.clin_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.clin = sl.cate.clin_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.clin"
  # return the output
  return(out)
}
predict.SL.cate.clin <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.clin, newdata = newdata, type = 'response')
  return(pred)
}

# clinical symptoms with pairwise interaction
SL.cate.clin2 <- function(Y, X, newX, family, ...){
  sl.cate.clin2_fit <- glm(Y ~ (dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + 
                                  dy1_scrn_diardays + dy1_scrn_dehydr)^2, 
                           data = X,
                           family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.clin2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.clin2 = sl.cate.clin2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.clin2"
  # return the output
  return(out)
}
predict.SL.cate.clin2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.clin2, newdata = newdata, type = 'response')
  return(pred)
}


# ------ Define CATE model superlearner library 
sl.library.CATE <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost", "SL.glmnet.twoway",
                     "SL.cate.clin", "SL.cate.clin2")

# ------ Define variables in the clinical symptoms rule 
Z_lists <- c("dy1_scrn_vomitall", "dy1_scrn_lstools", "dy1_scrn_sstools", "dy1_scrn_diardays",
             "dy1_scrn_dehydr")


# ----- Run estimation for child level expected benefit over five seeds

# seed 150
nuisance_output <- readRDS("nuisance150_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(150)

results_laz_clin150 <- estimate_OTR(df = abcd_data,
                                       Y_name = "lazdiff",
                                       A_name = "an_grp_01",
                                       W_list = W_list,
                                       Z_list = Z_lists,
                                       id_name = "pid",
                                       sl.library.CATE = sl.library.CATE,
                                       nuisance_models = nuisance_models,
                                       validRows = validRows,
                                       k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                       threshold = 0.06,
                                       k_folds = 10,
                                       ps_trunc_level = 0.01,
                                       outcome_type = "gaussian")

saveRDS(results_laz_clin150$results, "results_laz_clin150.RData")
saveRDS(results_laz_clin150$CATE_models, "catemodels_laz_clin150.RData")

# seed 300
nuisance_output <- readRDS("nuisance300_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(300)

results_laz_clin300 <- estimate_OTR(df = abcd_data,
                                    Y_name = "lazdiff",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = 0.06,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "gaussian")

saveRDS(results_laz_clin300$results, "results_laz_clin300.RData")
saveRDS(results_laz_clin300$CATE_models, "catemodels_laz_clin300.RData")

# seed 450
nuisance_output <- readRDS("nuisance450_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(450)

results_laz_clin450 <- estimate_OTR(df = abcd_data,
                                    Y_name = "lazdiff",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = 0.06,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "gaussian")

saveRDS(results_laz_clin450$results, "results_laz_clin450.RData")
saveRDS(results_laz_clin450$CATE_models, "catemodels_laz_clin450.RData")

# seed 600
nuisance_output <- readRDS("nuisance600_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(600)

results_laz_clin600 <- estimate_OTR(df = abcd_data,
                                    Y_name = "lazdiff",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = 0.06,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "gaussian")

saveRDS(results_laz_clin600$results, "results_laz_clin600.RData")
saveRDS(results_laz_clin600$CATE_models, "catemodels_laz_clin600.RData")

# seed 750
nuisance_output <- readRDS("nuisance750_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(750)

results_laz_clin750 <- estimate_OTR(df = abcd_data,
                                    Y_name = "lazdiff",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = 0.06,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "gaussian")

saveRDS(results_laz_clin750$results, "results_laz_clin750.RData")
saveRDS(results_laz_clin750$CATE_models, "catemodels_laz_clin750.RData")

