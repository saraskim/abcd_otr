# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Treatment effects for the Pathogens and Symptoms code for linear growth outcome
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

# ----- Define CATE models to estimate child level expected benefit for the Pathogens + Symptoms rule
# --- Include "nested" models (Shigella rule, Rotavirus rule, Pathogen Quantities rule, Symptoms rule, etc) that contain variables in the Pathogen Quantities + Symptoms rule 

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

# rotavirus pathogen quantity 
SL.cate.rota.glm <- function(Y, X, newX, family, ...){
  sl.cate.rota.glm_fit <- glm(Y ~ rotavirus_new + I(rotavirus_new > 0),
                              data = X,
                              family = family)
  pred <- predict(
    sl.cate.rota.glm_fit, newdata = newX, type = 'response'
  )
  fit <- list(fitted_model.cate.rota.glm = sl.cate.rota.glm_fit)
  out <- list(fit = fit, pred = pred)
  class(out$fit) <- "SL.cate.rota.glm"
  return(out)
}
predict.SL.cate.rota.glm <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.rota.glm, newdata = newdata, type = 'response')
  return(pred)
}

# pathogen quantities
SL.cate.path <- function(Y, X, newX, family, ...){
  sl.cate.path_fit <- glm(Y ~ rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0), 
                          data = X,
                          family = family)
  pred <- predict(
    sl.cate.path_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.path = sl.cate.path_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.path"
  # return the output
  return(out)
}
predict.SL.cate.path <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.path, newdata = newdata, type = 'response')
  return(pred)
}

# pathogen quantities with pairwise interaction
SL.cate.path2 <- function(Y, X, newX, family, ...){
  sl.cate.path2_fit <- glm(Y ~ (rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                                  I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                                  sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                                  st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                                  campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                                  v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                                  cryptosporidium_new + I(cryptosporidium_new > 0))^2, 
                           data = X,
                           family = family)
  pred <- predict(
    sl.cate.path2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.path2 = sl.cate.path2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.path2"
  # return the output
  return(out)
}
predict.SL.cate.path2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.path2, newdata = newdata, type = 'response')
  return(pred)
}

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

# pathogen quantities and clinical symptoms
SL.cate.pathclin <- function(Y, X, newX, family, ...){
  sl.cate.pathclin_fit <- glm(Y ~ rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                                I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                                sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                                st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                                campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                                v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                                cryptosporidium_new + I(cryptosporidium_new > 0) + dy1_scrn_vomitall + 
                                dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr, 
                              data = X,
                              family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.pathclin_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.pathclin = sl.cate.pathclin_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.pathclin"
  # return the output
  return(out)
}
predict.SL.cate.pathclin <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.pathclin, newdata = newdata, type = 'response')
  return(pred)
}

# pathogen quantities and clinical symptoms with pairwise interaction
SL.cate.pathclin2 <- function(Y, X, newX, family, ...){
  sl.cate.pathclin2_fit <- glm(Y ~ (rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                                      I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                                      sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                                      st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                                      campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                                      v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                                      cryptosporidium_new + I(cryptosporidium_new > 0) + dy1_scrn_vomitall + 
                                      dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr)^2, 
                               data = X,
                               family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.pathclin2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.pathclin2 = sl.cate.pathclin2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.pathclin2"
  # return the output
  return(out)
}
predict.SL.cate.pathclin2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.pathclin2, newdata = newdata, type = 'response')
  return(pred)
}


# ------ Define CATE model superlearner library 
sl.library.CATE <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost", "SL.glmnet.twoway",
                     "SL.cate.shig.glm", "SL.cate.rota.glm",
                     "SL.cate.path", "SL.cate.path2",
                     "SL.cate.clin", "SL.cate.clin2",
                     "SL.cate.pathclin", "SL.cate.pathclin2")

# ------ Define variables in the pathogens + symptoms rule 
Z_lists <- c("rotavirus_new", "rotavirus_bin", "norovirus_new", "norovirus_bin", "adenovirus_new", "adenovirus_bin",
             "astrovirus_new", "astrovirus_bin", "sapovirus_new", "sapovirus_bin", "st_etec_new", "st_etec_bin",
             "shigella_new", "shigella_bin", "campylobacter_new", "campylobacter_bin", "tepec_new", "tepec_bin",
             "v_cholerae_new", "v_cholerae_bin", "salmonella_new", "salmonella_bin", "cryptosporidium_new", "cryptosporidium_bin",
             "dy1_scrn_vomitall", "dy1_scrn_lstools", "dy1_scrn_sstools", "dy1_scrn_diardays",
             "dy1_scrn_dehydr")


# ----- Run estimation for child level expected benefit over five seeds

# seed 150
nuisance_output <- readRDS("nuisance150_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(150)

results_laz_pc150 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_laz_pc150$results, "results_laz_pc150.RData")
saveRDS(results_laz_pc150$CATE_models, "catemodels_laz_pc150.RData")

# seed 300
nuisance_output <- readRDS("nuisance300_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(300)

results_laz_pc300 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_laz_pc300$results, "results_laz_pc300.RData")
saveRDS(results_laz_pc300$CATE_models, "catemodels_laz_pc300.RData")

# seed 450
nuisance_output <- readRDS("nuisance450_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(450)

results_laz_pc450 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_laz_pc450$results, "results_laz_pc450.RData")
saveRDS(results_laz_pc450$CATE_models, "catemodels_laz_pc450.RData")

# seed 600
nuisance_output <- readRDS("nuisance600_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(600)

results_laz_pc600 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_laz_pc600$results, "results_laz_pc600.RData")
saveRDS(results_laz_pc600$CATE_models, "catemodels_laz_pc600.RData")

# seed 750
nuisance_output <- readRDS("nuisance750_laz_19jun2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(750)

results_laz_pc750 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_laz_pc750$results, "results_laz_pc750.RData")
saveRDS(results_laz_pc750$CATE_models, "catemodels_laz_pc750.RData")

