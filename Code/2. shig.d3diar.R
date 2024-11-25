# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Treatment effects for the SHIGELLA RULE code for DAY 3 DIARRHEA outcome
#
#
# Sara Kim
# 28 January 2024
#
# ------------------------------------------------------------------------------------------

# Please ensure "0. wrappers.d3diar.R" script has been run and outputs from "1. nuisance.d3diar.R" are in working directory

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




# ----- Define CATE models to estimate child level expected benefit for the Shigella rule
# shigella pathogen quantity 
SL.cate.shig.glm <- function(Y, X, newX, family, ...){
  sl.cate.shig.glm_fit <- glm(Y ~ shigella_new + I(shigella_new > 0),
                              data = X,
                              family = family)
  pred <- predict(
    sl.cate.shig.glm_fit, newdata = newX, type = 'response'
  )
  fit <- list(fitted_model.cate.shig.glm = sl.cate.shig.glm_fit)
  out <- list(fit = fit, pred = pred)
  class(out$fit) <- "SL.cate.shig.glm"
  return(out)
}
predict.SL.cate.shig.glm <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.shig.glm, newdata = newdata, type = 'response')
  return(pred)
}


# ------ Define CATE model superlearner library
sl.library.CATE <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost", "SL.glmnet.twoway",
                     "SL.cate.shig.glm")

# ------ Define variables in the Shigella rule 
Z_lists <- c("shigella_new", "shigella_bin")


# ----- Run estimation for child level expected benefit over five seeds

# seed 150
nuisance_output <- readRDS("nuisance150_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(150)

results_diar_shig150 <- estimate_OTR(df = abcd_data,
                                       Y_name = "day3diar",
                                       A_name = "an_grp_01",
                                       W_list = W_list,
                                       Z_list = Z_lists,
                                       id_name = "pid",
                                       sl.library.CATE = sl.library.CATE,
                                       nuisance_models = nuisance_models,
                                       validRows = validRows,
                                       k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                       threshold = -0.07,
                                       k_folds = 10,
                                       ps_trunc_level = 0.01,
                                       outcome_type = "binomial")

saveRDS(results_diar_shig150$results, "results_diar_shig150.RData")
saveRDS(results_diar_shig150$CATE_models, "catemodels_diar_shig150.RData")

# seed 300
nuisance_output <- readRDS("nuisance300_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(300)

results_diar_shig300 <- estimate_OTR(df = abcd_data,
                                     Y_name = "day3diar",
                                     A_name = "an_grp_01",
                                     W_list = W_list,
                                     Z_list = Z_lists,
                                     id_name = "pid",
                                     sl.library.CATE = sl.library.CATE,
                                     nuisance_models = nuisance_models,
                                     validRows = validRows,
                                     k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                     threshold = -0.07,
                                     k_folds = 10,
                                     ps_trunc_level = 0.01,
                                     outcome_type = "binomial")

saveRDS(results_diar_shig300$results, "results_diar_shig300.RData")
saveRDS(results_diar_shig300$CATE_models, "catemodels_diar_shig300.RData")

# seed 450
nuisance_output <- readRDS("nuisance450_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(450)

results_diar_shig450 <- estimate_OTR(df = abcd_data,
                                     Y_name = "day3diar",
                                     A_name = "an_grp_01",
                                     W_list = W_list,
                                     Z_list = Z_lists,
                                     id_name = "pid",
                                     sl.library.CATE = sl.library.CATE,
                                     nuisance_models = nuisance_models,
                                     validRows = validRows,
                                     k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                     threshold = -0.07,
                                     k_folds = 10,
                                     ps_trunc_level = 0.01,
                                     outcome_type = "binomial")

saveRDS(results_diar_shig450$results, "results_diar_shig450.RData")
saveRDS(results_diar_shig450$CATE_models, "catemodels_diar_shig450.RData")

# seed 600
nuisance_output <- readRDS("nuisance600_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(600)

results_diar_shig600 <- estimate_OTR(df = abcd_data,
                                     Y_name = "day3diar",
                                     A_name = "an_grp_01",
                                     W_list = W_list,
                                     Z_list = Z_lists,
                                     id_name = "pid",
                                     sl.library.CATE = sl.library.CATE,
                                     nuisance_models = nuisance_models,
                                     validRows = validRows,
                                     k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                     threshold = -0.07,
                                     k_folds = 10,
                                     ps_trunc_level = 0.01,
                                     outcome_type = "binomial")

saveRDS(results_diar_shig600$results, "results_diar_shig600.RData")
saveRDS(results_diar_shig600$CATE_models, "catemodels_diar_shig600.RData")

# seed 750
nuisance_output <- readRDS("nuisance750_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(750)

results_diar_shig750 <- estimate_OTR(df = abcd_data,
                                     Y_name = "day3diar",
                                     A_name = "an_grp_01",
                                     W_list = W_list,
                                     Z_list = Z_lists,
                                     id_name = "pid",
                                     sl.library.CATE = sl.library.CATE,
                                     nuisance_models = nuisance_models,
                                     validRows = validRows,
                                     k_fold_assign_and_CATE = k_fold_assign_and_CATE,
                                     threshold = -0.07,
                                     k_folds = 10,
                                     ps_trunc_level = 0.01,
                                     outcome_type = "binomial")

saveRDS(results_diar_shig750$results, "results_diar_shig750.RData")
saveRDS(results_diar_shig750$CATE_models, "catemodels_diar_shig750.RData")


