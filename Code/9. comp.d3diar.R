# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Treatment effects for the comprehensive rule code for DAY 3 DIARRHEA outcome
#
#
# Sara Kim
# 28 January 2024
#
# ------------------------------------------------------------------------------------------

# Please ensure "0. wrappers.d3diar.R" script has been run and outputs from "1. nuisance.d3diar.R" are in working directory

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
# --- Include "nested" models (pathogen quantities rule, symptoms rule, host rule, etc) that contain variables in the comprehensive rule 

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

# clinical symptoms, malnutrition, and sociodemographic  
SL.cate.nopath <- function(Y, X, newX, family, ...){
  sl.cate.nopath_fit <- glm(Y ~ dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + 
                              dy1_scrn_diardays + dy1_scrn_dehydr + avemuac + wfazscore + lfazscore + 
                              wflzscore + dy1_ant_sex + agemchild + an_ses_quintile + 
                              an_tothhlt5, 
                            data = X,
                            family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.nopath_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.nopath = sl.cate.nopath_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.nopath"
  # return the output
  return(out)
}
predict.SL.cate.nopath <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.nopath, newdata = newdata, type = 'response')
  return(pred)
}

# illness characteristics, malnutrition, and sociodemographic with pairwise interaction 
SL.cate.nopath2 <- function(Y, X, newX, family, ...){
  sl.cate.nopath2_fit <- glm(Y ~ (dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + 
                                    dy1_scrn_diardays + dy1_scrn_dehydr + avemuac + wfazscore + lfazscore + 
                                    wflzscore + dy1_ant_sex + agemchild + an_ses_quintile + 
                                    an_tothhlt5)^2, 
                             data = X,
                             family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.nopath2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.nopath2 = sl.cate.nopath2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.nopath2"
  # return the output
  return(out)
}
predict.SL.cate.nopath2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.nopath2, newdata = newdata, type = 'response')
  return(pred)
}

# pathogen quantity, illness characteristics, malnutrition, and sociodemographic  
SL.cate.all <- function(Y, X, newX, family, ...){
  sl.cate.all_fit <- glm(Y ~ rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                           I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                           sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                           st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                           campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                           v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                           cryptosporidium_new + I(cryptosporidium_new > 0) + dy1_scrn_vomitall + 
                           dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr + 
                           avemuac + wfazscore + lfazscore + wflzscore + dy1_ant_sex + agemchild + 
                           an_ses_quintile + an_tothhlt5, 
                         data = X,
                         family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.all_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.all = sl.cate.all_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.all"
  # return the output
  return(out)
}
predict.SL.cate.all <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.all, newdata = newdata, type = 'response')
  return(pred)
}

# pathogen quantity, illness characteristics, malnutrition, and sociodemographic with pairwise interaction
SL.cate.all2 <- function(Y, X, newX, family, ...){
  sl.cate.all2_fit <- glm(Y ~ (rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                                 I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                                 sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                                 st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                                 campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                                 v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                                 cryptosporidium_new + I(cryptosporidium_new > 0) + dy1_scrn_vomitall + 
                                 dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr + 
                                 avemuac + wfazscore + lfazscore + wflzscore + dy1_ant_sex + agemchild + 
                                 an_ses_quintile + an_tothhlt5)^2, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.cate.all2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.cate.all2 = sl.cate.all2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.cate.all2"
  # return the output
  return(out)
}
predict.SL.cate.all2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.cate.all2, newdata = newdata, type = 'response')
  return(pred)
}


# ------ Define CATE model superlearner library
sl.library.CATE <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.earth", "SL.xgboost", "SL.glmnet.twoway",
                     "SL.cate.shig.glm", "SL.cate.rota.glm",
                     "SL.cate.clin", "SL.cate.clin2",
                     "SL.cate.pathclin", "SL.cate.pathclin2",
                     "SL.cate.malsoc", "SL.cate.malsoc2",
                     "SL.cate.nopath", "SL.cate.nopath2",
                     "SL.cate.all", "SL.cate.all2")

# ------ Define variables in the comprehensive rule 
Z_lists <- c("rotavirus_new", "rotavirus_bin", "norovirus_new", "norovirus_bin", "adenovirus_new", "adenovirus_bin",
             "astrovirus_new", "astrovirus_bin", "sapovirus_new", "sapovirus_bin", "st_etec_new", "st_etec_bin",
             "shigella_new", "shigella_bin", "campylobacter_new", "campylobacter_bin", "tepec_new", "tepec_bin",
             "v_cholerae_new", "v_cholerae_bin", "salmonella_new", "salmonella_bin", "cryptosporidium_new", "cryptosporidium_bin",
             "dy1_scrn_vomitall", "dy1_scrn_lstools", "dy1_scrn_sstools", "dy1_scrn_diardays",
             "dy1_scrn_dehydr", "avemuac", "wfazscore", "lfazscore", "wflzscore", 
             "dy1_ant_sex", "agemchild", "an_ses_quintile", "an_tothhlt5")

# ----- Run estimation for child level expected benefit over five seeds

# seed 150
nuisance_output <- readRDS("nuisance150_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(150)

results_diar_all150 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_diar_all150$results, "results_diar_all150.RData")
saveRDS(results_diar_all150$CATE_models, "catemodels_diar_all150.RData")

# seed 300
nuisance_output <- readRDS("nuisance300_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(300)

results_diar_all300 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_diar_all300$results, "results_diar_all300.RData")
saveRDS(results_diar_all300$CATE_models, "catemodels_diar_all300.RData")

# seed 450
nuisance_output <- readRDS("nuisance450_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(450)

results_diar_all450 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_diar_all450$results, "results_diar_all450.RData")
saveRDS(results_diar_all450$CATE_models, "catemodels_diar_all450.RData")

# seed 600
nuisance_output <- readRDS("nuisance600_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(600)

results_diar_all600 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_diar_all600$results, "results_diar_all600.RData")
saveRDS(results_diar_all600$CATE_models, "catemodels_diar_all600.RData")

# seed 750
nuisance_output <- readRDS("nuisance750_day3diar_30apr2024.RData")
nuisance_models <- nuisance_output$nuisance_models
k_fold_assign_and_CATE <- nuisance_output$k_fold_assign_and_CATE
validRows <- nuisance_output$validRows

set.seed(NULL)
set.seed(750)

results_diar_all750 <- estimate_OTR(df = abcd_data,
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

saveRDS(results_diar_all750$results, "results_diar_all750.RData")
saveRDS(results_diar_all750$CATE_models, "catemodels_diar_all750.RData")
