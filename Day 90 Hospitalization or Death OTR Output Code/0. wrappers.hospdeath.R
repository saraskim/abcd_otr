# -------------------------------------------------------------------------------------------
# 
# Code for ABCD Trial personalized treatment effects 
# Wrapper code and nuisance model superlearner libraries for DAY 90 RE-HOSPITALIZATION & DEATH outcome
# Make sure this code is run first prior to running analyses
#
#
# Sara Kim
# 28 January 2024
#
# ------------------------------------------------------------------------------------------

# ----- Load glmnet two way model 
.SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

SL.glmnet.twoway <- function (
    Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
    nlambda = 100, useMin = TRUE, loss = "deviance", ...
){
  .SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + .^2, X)
    newX <- model.matrix(~-1 + .^2, newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda,
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet.twoway"
  out <- list(pred = pred, fit = fit)
  return(out)
}
predict.SL.glmnet.twoway <- function (
    object, newdata, remove_extra_cols = TRUE,
    add_missing_cols = TRUE, ...
){
  .SL.require("glmnet")
  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~-1 + .^2, newdata)
  }
  original_cols = rownames(object$object$glmnet.fit$beta)
  if (remove_extra_cols) {
    extra_cols = setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste("Removing extra columns in prediction data:",
                    paste(extra_cols, collapse = ", ")))
      newdata = newdata[, !colnames(newdata) %in% extra_cols,
                        drop = FALSE]
    }
  }
  if (add_missing_cols) {
    missing_cols = setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste("Adding missing columns in prediction data:",
                    paste(missing_cols, collapse = ", ")))
      new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
      colnames(new_cols) = missing_cols
      newdata = cbind(newdata, new_cols)
      newdata = newdata[, original_cols, drop = FALSE]
    }
  }
  pred <- predict(object$object, newx = newdata, type = "response",
                  s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
  return(pred)
}


# ----- Define outcome models

# interaction for site and rotavirus 
SL.outcome.1 <- function(Y, X, newX, family, obsWeights, ...){
  sl.outcome.1_fit <- glm(Y ~ an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + 
                            dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site + an_grp_01*rotavirus_new +
                            + an_grp_01*I(rotavirus_new > 0), 
                          data = X,
                          family = family,
                          weights = obsWeights)
  pred <- predict(sl.outcome.1_fit, newdata = newX, type = 'response')
  fit <- list(fitted_model.outcome.1 = sl.outcome.1_fit)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- "SL.outcome.1"
  return(out)
}
predict.SL.outcome.1 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.1, newdata = newdata, type = 'response')
  return(pred)
}



# interaction for site and shigella
SL.outcome.2 <- function(Y, X, newX, family, ...){
  sl.outcome.2_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*shigella_new + an_grp_01*I(shigella_new > 0), 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.2 = sl.outcome.2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.2"
  # return the output
  return(out)
}
predict.SL.outcome.2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.2, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site and rotavirus SEASON
SL.outcome.3 <- function(Y, X, newX, family, ...){
  sl.outcome.3_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*rotaseason, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.3_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.3 = sl.outcome.3_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.3"
  # return the output
  return(out)
}
predict.SL.outcome.3 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.3, newdata = newdata, type = 'response')
  return(pred)
}
# rotaseason through interaction by month
SL.outcome.11 <- function(Y, X, newX, family, ...){
  sl.outcome.11_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                             I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                             sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                             st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                             campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                             v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                             cryptosporidium_new + I(cryptosporidium_new > 0) + 
                             dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                             avemuac + wfazscore + lfazscore + wflzscore +
                             site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                             month_en + rotaseason + an_grp_01*site +
                             an_grp_01*month_en, 
                           data = X,
                           family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.11_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.11 = sl.outcome.11_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.11"
  # return the output
  return(out)
}
predict.SL.outcome.11 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.11, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site, rotavirus, and shigella
SL.outcome.4 <- function(Y, X, newX, family, ...){
  sl.outcome.4_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*rotavirus_new + an_grp_01*I(rotavirus_new > 0) + 
                            an_grp_01*shigella_new + an_grp_01*I(shigella_new > 0), 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.4_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.4 = sl.outcome.4_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.4"
  # return the output
  return(out)
}
predict.SL.outcome.4 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.4, newdata = newdata, type = 'response')
  return(pred)
}



# interaction for site and 12 pathogen quantities 
SL.outcome.5 <- function(Y, X, newX, family, ...){
  sl.outcome.5_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*rotavirus_new + an_grp_01*I(rotavirus_new > 0) + an_grp_01*norovirus_new + 
                            an_grp_01*I(norovirus_new > 0) + an_grp_01*adenovirus_new + an_grp_01*I(adenovirus_new > 0) +
                            an_grp_01*astrovirus_new + an_grp_01*I(astrovirus_new > 0) + an_grp_01*sapovirus_new + 
                            an_grp_01*I(sapovirus_new > 0) + an_grp_01*st_etec_new + an_grp_01*I(st_etec_new > 0) +
                            an_grp_01*shigella_new + an_grp_01*I(shigella_new > 0) + an_grp_01*campylobacter_new + 
                            an_grp_01*I(campylobacter_new > 0) + an_grp_01*tepec_new + an_grp_01*I(tepec_new > 0) +
                            an_grp_01*v_cholerae_new + an_grp_01*I(v_cholerae_new > 0) + an_grp_01*salmonella_new + 
                            an_grp_01*I(salmonella_new > 0) + an_grp_01*cryptosporidium_new + an_grp_01*I(cryptosporidium_new > 0), 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.5_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.5 = sl.outcome.5_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.5"
  # return the output
  return(out)
}
predict.SL.outcome.5 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.5, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site and vomiting
SL.outcome.6 <- function(Y, X, newX, family, ...){
  sl.outcome.6_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*dy1_scrn_vomitall, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.6_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.6 = sl.outcome.6_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.6"
  # return the output
  return(out)
}
predict.SL.outcome.6 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.6, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site, 12 pathogens, and illness characteristics
SL.outcome.7 <- function(Y, X, newX, family, ...){
  sl.outcome.7_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*rotavirus_new + an_grp_01*I(rotavirus_new > 0) + an_grp_01*norovirus_new + 
                            an_grp_01*I(norovirus_new > 0) + an_grp_01*adenovirus_new + an_grp_01*I(adenovirus_new > 0) +
                            an_grp_01*astrovirus_new + an_grp_01*I(astrovirus_new > 0) + an_grp_01*sapovirus_new + 
                            an_grp_01*I(sapovirus_new > 0) + an_grp_01*st_etec_new + an_grp_01*I(st_etec_new > 0) +
                            an_grp_01*shigella_new + an_grp_01*I(shigella_new > 0) + an_grp_01*campylobacter_new + 
                            an_grp_01*I(campylobacter_new > 0) + an_grp_01*tepec_new + an_grp_01*I(tepec_new > 0) +
                            an_grp_01*v_cholerae_new + an_grp_01*I(v_cholerae_new > 0) + an_grp_01*salmonella_new + 
                            an_grp_01*I(salmonella_new > 0) + an_grp_01*cryptosporidium_new + an_grp_01*I(cryptosporidium_new > 0) +
                            an_grp_01*dy1_scrn_vomitall + an_grp_01*dy1_scrn_lstools + an_grp_01*dy1_scrn_sstools +
                            an_grp_01*dy1_scrn_diardays + an_grp_01*dy1_scrn_dehydr, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.7_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.7 = sl.outcome.7_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.7"
  # return the output
  return(out)
}
predict.SL.outcome.7 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.7, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site, malnutrition, and sociodemographics
SL.outcome.8 <- function(Y, X, newX, family, ...){
  sl.outcome.8_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*avemuac + an_grp_01*wfazscore + an_grp_01*lfazscore +
                            an_grp_01*wflzscore + an_grp_01*dy1_ant_sex + an_grp_01*agemchild +
                            an_grp_01*an_ses_quintile + an_grp_01*an_tothhlt5, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.8_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.8 = sl.outcome.8_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.8"
  # return the output
  return(out)
}
predict.SL.outcome.8 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.8, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site, illness characteristics, malnutrition, and sociodemographics
SL.outcome.9 <- function(Y, X, newX, family, ...){
  sl.outcome.9_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                            avemuac + wfazscore + lfazscore + wflzscore +
                            site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                            month_en + rotaseason + an_grp_01*site +
                            an_grp_01*dy1_scrn_vomitall + an_grp_01*dy1_scrn_lstools + 
                            an_grp_01*dy1_scrn_sstools + an_grp_01*dy1_scrn_diardays +
                            an_grp_01*dy1_scrn_dehydr +
                            an_grp_01*avemuac + an_grp_01*wfazscore + an_grp_01*lfazscore +
                            an_grp_01*wflzscore + an_grp_01*dy1_ant_sex + an_grp_01*agemchild +
                            an_grp_01*an_ses_quintile + an_grp_01*an_tothhlt5, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.9_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.9 = sl.outcome.9_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.9"
  # return the output
  return(out)
}
predict.SL.outcome.9 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.9, newdata = newdata, type = 'response')
  return(pred)
}


# interaction for site, 12 pathogens, illness characteristics, malnutrition, and sociodemographics
SL.outcome.10 <- function(Y, X, newX, family, ...){
  sl.outcome.10_fit <- glm(Y ~  an_grp_01 + rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                             I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                             sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                             st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                             campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                             v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                             cryptosporidium_new + I(cryptosporidium_new > 0) + 
                             dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + dy1_scrn_dehydr +
                             avemuac + wfazscore + lfazscore + wflzscore +
                             site + dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 +
                             month_en + rotaseason + an_grp_01*site +
                             an_grp_01*rotavirus_new + an_grp_01*I(rotavirus_new > 0) + an_grp_01*norovirus_new + 
                             an_grp_01*I(norovirus_new > 0) + an_grp_01*adenovirus_new + an_grp_01*I(adenovirus_new > 0) +
                             an_grp_01*astrovirus_new + an_grp_01*I(astrovirus_new > 0) + an_grp_01*sapovirus_new + 
                             an_grp_01*I(sapovirus_new > 0) + an_grp_01*st_etec_new + an_grp_01*I(st_etec_new > 0) +
                             an_grp_01*shigella_new + an_grp_01*I(shigella_new > 0) + an_grp_01*campylobacter_new + 
                             an_grp_01*I(campylobacter_new > 0) + an_grp_01*tepec_new + an_grp_01*I(tepec_new > 0) +
                             an_grp_01*v_cholerae_new + an_grp_01*I(v_cholerae_new > 0) + an_grp_01*salmonella_new + 
                             an_grp_01*I(salmonella_new > 0) + an_grp_01*cryptosporidium_new + an_grp_01*I(cryptosporidium_new > 0) +
                             an_grp_01*dy1_scrn_vomitall + an_grp_01*dy1_scrn_lstools + 
                             an_grp_01*dy1_scrn_sstools + an_grp_01*dy1_scrn_diardays +
                             an_grp_01*dy1_scrn_dehydr +
                             an_grp_01*avemuac + an_grp_01*wfazscore + an_grp_01*lfazscore +
                             an_grp_01*wflzscore + an_grp_01*dy1_ant_sex + an_grp_01*agemchild +
                             an_grp_01*an_ses_quintile + an_grp_01*an_tothhlt5, 
                           data = X,
                           family = family)
  # get predictions on newX
  pred <- predict(
    sl.outcome.10_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.outcome.10 = sl.outcome.10_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.outcome.10"
  # return the output
  return(out)
}
predict.SL.outcome.10 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.outcome.10, newdata = newdata, type = 'response')
  return(pred)
}



# ----- Define treatment models
SL.treatment <- function(Y, X, newX, family, ...){
  sl.treatment_fit <- glm(Y ~ rotavirus_new + I(rotavirus_new > 0) + norovirus_new + 
                            I(norovirus_new > 0) + adenovirus_new + I(adenovirus_new > 0) +
                            sapovirus_new + I(sapovirus_new > 0) + astrovirus_new + I(astrovirus_new > 0) +
                            st_etec_new + I(st_etec_new > 0) + shigella_new + I(shigella_new > 0) +
                            campylobacter_new + I(campylobacter_new > 0) + tepec_new + I(tepec_new > 0) +
                            v_cholerae_new + I(v_cholerae_new > 0) + salmonella_new + I(salmonella_new > 0) +
                            cryptosporidium_new + I(cryptosporidium_new > 0) + 
                            dy1_scrn_vomitall + dy1_scrn_lstools + dy1_scrn_sstools + dy1_scrn_diardays + 
                            dy1_scrn_dehydr + avemuac + wfazscore + lfazscore + wflzscore + site + 
                            dy1_ant_sex + agemchild + an_ses_quintile + an_tothhlt5 + month_en + rotaseason, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.treatment_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.treatment = sl.treatment_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.treatment"
  # return the output
  return(out)
}
predict.SL.treatment <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.treatment, newdata = newdata, type = 'response')
  return(pred)
}


# ----- define missingness models
SL.missing.1 <- function(Y, X, newX, family, ...){
  sl.missing.1_fit <- glm(Y ~ site + month_en, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.missing.1_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.missing.1 = sl.missing.1_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.missing.1"
  # return the output
  return(out)
}
predict.SL.missing.1 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.missing.1, newdata = newdata, type = 'response')
  return(pred)
}


SL.missing.2 <- function(Y, X, newX, family, ...){
  sl.missing.2_fit <- glm(Y ~ site + month_en + site*month_en, 
                          data = X,
                          family = family)
  # get predictions on newX
  pred <- predict(
    sl.missing.2_fit, newdata = newX, type = 'response'
  )
  # format the output as named list
  fit <- list(fitted_model.missing.2 = sl.missing.2_fit)
  out <- list(fit = fit, pred = pred)
  # give the object a class
  class(out$fit) <- "SL.missing.2"
  # return the output
  return(out)
}
predict.SL.missing.2 <- function(object, newdata, ...){
  pred <- predict(object$fitted_model.missing.2, newdata = newdata, type = 'response')
  return(pred)
}


# ----- Load W List of covariates
W_list <- c("rotavirus_new", "rotavirus_bin", "norovirus_new", "norovirus_bin", "adenovirus_new", 
            "adenovirus_bin", "sapovirus_new","sapovirus_bin", "astrovirus_new", "astrovirus_bin",
            "st_etec_new", "st_etec_bin", "shigella_new", "shigella_bin", "campylobacter_new", 
            "campylobacter_bin", "tepec_new", "tepec_bin", "v_cholerae_new", "v_cholerae_bin", 
            "salmonella_new", "salmonella_bin", "cryptosporidium_new", "cryptosporidium_bin",
            "dy1_scrn_vomitall", "dy1_scrn_lstools", "dy1_scrn_sstools", "dy1_scrn_diardays", 
            "dy1_scrn_dehydr", "avemuac", "wfazscore", "lfazscore", "wflzscore", "site", 
            "dy1_ant_sex", "agemchild", "an_ses_quintile", "an_tothhlt5", "month_en", "rotaseason")


# ----- define super learner libraries 
sl.library.outcome <- c("SL.glm", "SL.glmnet", "SL.glmnet.twoway", "SL.ranger", "SL.earth", 
                        "SL.xgboost", "SL.outcome.1", "SL.outcome.2", "SL.outcome.3", 
                        "SL.outcome.11", "SL.outcome.4", "SL.outcome.5", "SL.outcome.6", 
                        "SL.outcome.7", "SL.outcome.8", "SL.outcome.9", "SL.outcome.10")      
sl.library.treatment <- c("SL.mean", "SL.glm", "SL.treatment")  
sl.library.missingness <- c("SL.mean", "SL.missing.1", "SL.missing.2")
