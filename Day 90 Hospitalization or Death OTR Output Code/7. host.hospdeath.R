# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Treatment effects for the host (malnutrition + sociodemographics) rule code for DAY 90 REHOSPITALIZATION OR DEATH outcome
#
#
# Sara Kim
# 28 January 2024
#
# ------------------------------------------------------------------------------------------

# Please ensure "0. wrappers.hospdeath.R" script has been run and outputs from "1. nuisance.hospdeath.R" are in working directory

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


# ----- Define CATE models to estimate child level expected benefit for the host rule

# malnutrition and sociodemographics  
SL.cate.malsoc <- function(Y, X, newX, family, ...){
  sl.cate.malsoc_fit <- glm(Y ~ avemuac + wfazscore + lfazscore + wflzscore + 
                              dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5, 
                            data = X,
                            family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.malsoc_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.malsoc = sl.cate.malsoc_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.malsoc"
  # return the output
  return(out)
}
predict.SL.cate.malsoc <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.malsoc, newdata = newdata, type = 'response')
  return(pred)
}

# malnutrition and sociodemographic with pairwise interaction 
SL.cate.malsoc2 <- function(Y, X, newX, family, ...){
  sl.cate.malsoc2_fit <- glm(Y ~ (avemuac + wfazscore + lfazscore + wflzscore + 
                                    dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5)^2, 
                             data = X,
                             family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.malsoc2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.malsoc2 = sl.cate.malsoc2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.malsoc2"
  # return the output
  return(out)
}
predict.SL.cate.malsoc2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.malsoc2, newdata = newdata, type = 'response')
  return(pred)
}


# ------ Define CATE model superlearner library 
sl.library.CATE <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost", "SL.glmnet.twoway",
                     "SL.cate.malsoc", "SL.cate.malsoc2")

# ------ Define variables in the host rule 
Z_lists <- c("avemuac", "wfazscore", "lfazscore", "wflzscore", 
             "dy1_ant_sex", "agemchild", "an_ses_quintile", "an_tothhlt5")


# ----- Run estimation for child level expected benefit over five seeds

# seed 150
nuisance_output <- readRDS("nuisance150_hospdeath_24apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(150)

results_hosp_mal150 <- estimate_OTR(df = abcd_data,
                                    Y_name = "an_hosp90death",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = -0.02,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "binomial")

saveRDS(results_hosp_mal150$results, "results_hosp_mal150.RData")
saveRDS(results_hosp_mal150$CATE_models, "catemodels_hosp_mal150.RData")

# seed 300
nuisance_output <- readRDS("nuisance300_hospdeath_24apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(300)

results_hosp_mal300 <- estimate_OTR(df = abcd_data,
                                    Y_name = "an_hosp90death",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = -0.02,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "binomial")

saveRDS(results_hosp_mal300$results, "results_hosp_mal300.RData")
saveRDS(results_hosp_mal300$CATE_models, "catemodels_hosp_mal300.RData")

# seed 450
nuisance_output <- readRDS("nuisance450_hospdeath_24apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(450)

results_hosp_mal450 <- estimate_OTR(df = abcd_data,
                                    Y_name = "an_hosp90death",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = -0.02,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "binomial")

saveRDS(results_hosp_mal450$results, "results_hosp_mal450.RData")
saveRDS(results_hosp_mal450$CATE_models, "catemodels_hosp_mal450.RData")

# seed 600
nuisance_output <- readRDS("nuisance600_hospdeath_24apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(600)

results_hosp_mal600 <- estimate_OTR(df = abcd_data,
                                    Y_name = "an_hosp90death",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = -0.02,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "binomial")

saveRDS(results_hosp_mal600$results, "results_hosp_mal600.RData")
saveRDS(results_hosp_mal600$CATE_models, "catemodels_hosp_mal600.RData")

# seed 750
nuisance_output <- readRDS("nuisance750_hospeath_24apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(750)

results_hosp_mal750 <- estimate_OTR(df = abcd_data,
                                    Y_name = "an_hosp90death",
                                    A_name = "an_grp_01",
                                    W_list = W_list,
                                    Z_list = Z_lists,
                                    id_name = "pid",
                                    sl.library.CATE = sl.library.CATE,
                                    nuisance_models = nuisance_models,
                                    validRows = validRows,
                                    k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                    threshold = -0.02,
                                    k_folds = 10,
                                    ps_trunc_level = 0.01,
                                    outcome_type = "binomial")

saveRDS(results_hosp_mal750$results, "results_hosp_mal750.RData")
saveRDS(results_hosp_mal750$CATE_models, "catemodels_hosp_mal750.RData")

