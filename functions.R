library(tidyverse)
library(mice)
library(pROC)
library(caret)
library(cowplot)

create_cvfolds <- function(y, k=5, seed=12345){
  set.seed(seed)
  createFolds(y,k,list=TRUE,returnTrain=F)
}

table_cat <-  function(data, col, caption = "", level = c(), type) {
  col <- enquo(col)
  N <- dim(data)[1]
  N_event <- sum(data$y==1)
  smart %>% 
    mutate(col = !! col) %>%
    mutate(col = ifelse(is.na(col), 9, col)) %>%
    group_by(col, y) %>%
    summarise(pt = n()) %>%
    group_by(col) %>%
    mutate(N = sum(pt)) %>%
    mutate(pct = pt/N_event*100) %>%
    gather("metric", "value", -col, -y) %>%
    unite(metric,y,metric) %>%
    spread(metric, value) %>%
    mutate(all_pt = TRUE_N, all_pct = TRUE_N/N*100) %>%
    select(all_pt,all_pct, TRUE_pt, TRUE_pct) -> tmp
  if (type == "bi") {
    tmp %>%
      filter(col == 1) %>%
      mutate(col = caption)
  } else if (type == "multi") {
    tmp %>%
      filter(col != 9) %>%
      mutate(col = recode(col, !!! level)) %>% 
      ungroup %>%
      add_row(.before=T) %>%
      mutate(col=ifelse(is.na(col),caption, paste("\t- ", col, sep="")))
  }
}

table_num <- function(data, col, caption) {
  N <- dim(data)[1]
  col <- enquo(col)
  mean_col = summarise(data, mean(!! col, na.rm=T))
  sd_col = summarise(data, sd(!! col, na.rm=T))
  data %>%
    group_by(y) %>%
    summarise(pt=mean(!! col, na.rm=T), pct=sd(!! col, na.rm=T)) %>%
    gather(metric, value, -y) %>%
    unite(metric, y, metric) %>%
    spread(metric, value) %>%
    mutate(col=caption, all_pt=mean_col[1,1], all_pct=sd_col[1,1]) %>%
    select(col, all_pt, all_pct, TRUE_pt, TRUE_pct)
}

brier <- function(y,y_pred) {
  bs <- mean((y_pred-y)^2)
  return(bs)
}

calibration_glm <- function(y, y_pred) {
  y_lp <- log(y_pred/(1-y_pred))
  fit <- glm(y ~ y_lp, family="binomial")
  fit_offset <- glm(y ~ offset(y_lp), family="binomial")
  intercept <- fit_offset$coef
  slope <- fit$coef[2]
  results = c(intercept, slope)
  names(results) <- c("intercept", "slope")
  return(results)
}

original_model <- function(data, label) {
  if (label == "d.mort") {
    fit <- glm(d.mort~age+I(age^2)+d.motor_gr+d.pupil+hypoxia+hypotens+ctclass_gr+tsah, data=data, family="binomial")
  } else if (label == "d.unfav") {
    fit <- glm(d.unfav~age+I(age^2)+d.motor_gr+d.pupil+hypoxia+hypotens+ctclass_gr+tsah, data=data, family="binomial")
  }
  return(fit)
}

get_h2o_from_parameters <- function(parameters) {
  model_string <- str_split(parameters$model_id, pattern="_")[[1]][1]
  if (model_string == "GLM") {
    return(h2o.glm)
  } else if (model_string == "GBM") {
    return(h2o.gbm)
  } else if (model_string == "XGBoost") {
    return(h2o.xgboost)
  } else if (model_string == "DeepLearning") {
    return (h2o.deeplearning)
  } else if (model_string == "DRF") {
    return (h2o.randomForest) 
  } else if (model_string == "XRT") {
    return (h2o.randomForest)
  }
}

cv_h2o <- function(data, model_mort, model_unfav, cv_folds) {
  k <- length(cv_folds)
  aucs <- matrix(nrow=k, ncol=2)
  colnames(aucs) <- c("Death", "Unfavorable outcomes")
  calibrations_mort_h2o <- c()
  calibrations_unfav_h2o <- c()
  calibrations_glm_mort_h2o <- matrix(ncol=2, nrow=k)
  calibrations_glm_unfav_h2o <- matrix(ncol=2, nrow=k)
  briers <- matrix(ncol=2, nrow=k)
  
  for (i in 1:k) {
    print(i)
    test_index <- cv_folds[[i]]
    train <- as.h2o(data[-test_index, ])
    test <- data[test_index, ]
    model_mort$training_frame <- train
    model_unfav$training_frame <- train
    fit_mort_h2o <- do.call(get_h2o_from_parameters(model_mort), model_mort)
    fit_unfav_h2o <- do.call(get_h2o_from_parameters(model_unfav), model_unfav)
    mort_predict_h2o_df <- h2o.predict(fit_mort_h2o, as.h2o(test))
    mort_predict_h2o <- as.data.frame(mort_predict_h2o_df)[,"p1"]
    unfav_predict_h2o_df <- h2o.predict(fit_unfav_h2o, as.h2o(test))
    unfav_predict_h2o <- as.data.frame(unfav_predict_h2o_df)[,"p1"]
    
    data.frame(observed=test$d.mort, predicted=mort_predict_h2o) %>%
      mutate(cut=cut(predicted,10)) %>%
      group_by(cut) %>% summarise(observed=mean(observed), predicted=mean(predicted)) %>%
      mutate(fold=factor(i), ith=1:n()) %>%
      rbind(calibrations_mort_h2o) -> calibrations_mort_h2o
    
    data.frame(observed=test$d.unfav, predicted=unfav_predict_h2o) %>%
      mutate(cut=cut(predicted,10)) %>%
      group_by(cut) %>% summarise(observed=mean(observed), predicted=mean(predicted)) %>%
      mutate(fold=factor(i), ith=1:n()) %>%
      rbind(calibrations_unfav_h2o) -> calibrations_unfav_h2o
    
    calibrations_glm_mort_h2o[i,] <- calibration_glm(test$d.mort, mort_predict_h2o)
    calibrations_glm_unfav_h2o[i,] <- calibration_glm(test$d.unfav, unfav_predict_h2o)
    colnames(calibrations_glm_mort_h2o) <- c("intercept", "slope")
    colnames(calibrations_glm_unfav_h2o) <- c("intercept", "slope")
    
    aucs[i,1] <- auc(test$d.mort ~ mort_predict_h2o)
    aucs[i,2] <- auc(test$d.unfav ~ unfav_predict_h2o)
    
    briers[i,1] <- brier(test$d.mort, mort_predict_h2o)
    briers[i,2] <- brier(test$d.mort, mort_predict_h2o)
  }
  calibrations_mort_h2o <- calibrations_mort_h2o %>%
    group_by(ith) %>% summarise(observed=mean(observed), predicted=mean(predicted))
  calibrations_unfav_h2o <- calibrations_unfav_h2o %>%
    group_by(ith) %>% summarise(observed=mean(observed), predicted=mean(predicted))
  results <- list(
    auc=colMeans(aucs), 
    calibrations=list(mort=calibrations_mort_h2o, unfav=calibrations_unfav_h2o),
    calibrations_glm=list(mort=colMeans(calibrations_glm_mort_h2o), unfav=colMeans(calibrations_glm_unfav_h2o)),
    briers=colMeans(briers)
  )
  return(results)
}

original_study_lp <- list(
  mort=c(-3.267,-0.0198,0.000528,1.126,0.918,0.494,0.364,0.648,0.745,0.377,0.808,1.354,0.860,0.694),
  unfav=c(-2.842,0.00106,0.000391,1.690,1.275,0.675,0.440,0.938,0.449,0.740,0.637,0.628,0.619,0.657)
)


# AUC and calibration and brier with 5 cv
cv_original <- function(data, cv_folds) {
  k <- length(cv_folds)
  AUCs <- matrix(nrow=k, ncol=2)
  calibrations_mort <- c()
  calibrations_unfav <- c()
  calibrations_glm_mort <- matrix(ncol=2, nrow=k)
  calibrations_glm_unfav <- matrix(ncol=2, nrow=k)
  colnames(AUCs) <- c("Death", "Unfavorable outcomes")
  colnames(calibrations_glm_mort) <- c("Intercept", "Slope")
  colnames(calibrations_glm_unfav) <- c("Intercept", "Slope")
  briers <- matrix(nrow=k, ncol=2)
  
  for (i in 1:k) {
    test_index <- cv_folds[[i]]
    train <- data[-(data$.id %in% test_index), ]
    test <- data[data$.id %in% test_index, ]
    fit_mort <- original_model(train, "d.mort")
    fit_unfav <- original_model(train, "d.unfav")
    pred_mort <- predict(fit_mort, newdata=test, type="response")
    pred_unfav <- predict(fit_unfav, newdata=test, type="response")
    AUCs[i,1] <- auc(test$d.mort ~ pred_mort)
    AUCs[i,2] <- auc(test$d.unfav ~ pred_unfav)
    
    calibration_mort <- data.frame(observed=test$d.mort, predicted=pred_mort) %>%
      mutate(cut=cut(predicted,10)) %>%
      group_by(cut) %>% summarise(observed=mean(observed), predicted=mean(predicted)) %>%
      mutate(fold=factor(i), ith=1:10)
    
    calibration_unfav <- data.frame(observed=test$d.unfav, predicted=pred_unfav) %>%
      mutate(cut=cut(predicted,10)) %>%
      group_by(cut) %>% summarise(observed=mean(observed), predicted=mean(predicted)) %>%
      mutate(fold=factor(i), ith=1:10)
    
    calibrations_mort <- rbind(calibrations_mort, calibration_mort)
    calibrations_unfav <- rbind(calibrations_unfav, calibration_unfav)
    
    calibrations_glm_mort[i,] <- calibration_glm(test$d.mort, pred_mort)
    calibrations_glm_unfav[i,] <- calibration_glm(test$d.unfav, pred_unfav)
    
    briers[i,] <- brier(test$d.mort, pred_mort)
    briers[i,] <- brier(test$d.unfav, pred_unfav)
  }
  calibrations_mort <-  calibrations_mort %>%
    group_by(ith) %>% summarise(observed=mean(observed), predicted=mean(predicted))
  calibrations_unfav <- calibrations_unfav %>%
    group_by(ith) %>% summarise(observed=mean(observed), predicted=mean(predicted))
  
  results <- list(
    auc=colMeans(AUCs),
    calibrations=list(mort=calibrations_mort, unfav=calibration_unfav),
    calibrations_glm=list(mort=colMeans(calibrations_glm_mort), unfav=colMeans(calibrations_glm_unfav)),
    briers=colMeans(briers)
  )
  return(results)
}

createCalibrationPlotForH2oOriginal_mort <- function(results_h2o, results_original, method="glm") {
  fig <- ggplot() +
    geom_point(results_h2o$calibrations$mort, mapping=aes(x=predicted, y=observed), color='blue') + 
    geom_smooth(results_h2o$calibrations$mort, mapping=aes(x=predicted, y=observed, color="H2O"), method=method, se=F) + 
    geom_point(results_original$calibrations$mort, mapping=aes(x=predicted, y=observed), color='red') + 
    geom_smooth(results_original$calibrations$mort, mapping=aes(x=predicted, y=observed, color="Original"), method=method, se=F) + 
    geom_abline(intercept=0, slope=1) + 
    scale_colour_manual(name="Legend", values=cols) +
    labs(title="Mortality") +
    xlim(c(0,1)) + ylim(c(0,1))
}

createCalibrationPlotForH2oOriginal_unfav <- function(results_h2o, results_original, method="glm") {
  fig <- ggplot() +
    geom_point(results_h2o$calibrations$unfav, mapping=aes(x=predicted, y=observed), color='blue') + 
    geom_smooth(results_h2o$calibrations$unfav, mapping=aes(x=predicted, y=observed, color="H2O"), method=method, se=F) + 
    geom_point(results_original$calibrations$unfav, mapping=aes(x=predicted, y=observed), color='red') + 
    geom_smooth(results_original$calibrations$unfav, mapping=aes(x=predicted, y=observed, color="Original"), method=method, se=F) + 
    geom_abline(intercept=0, slope=1) + 
    scale_colour_manual(name="Legend", values=cols) +
    labs(title="Unfavorable outcomes") +
    xlim(c(0,1)) + ylim(c(0,1))
}

get_model_from_leaderboard <- function(leaderboard, model_path, idx, y) {
  model <- h2o.loadModel(paste(model_path, leaderboard[idx,1], sep="/"))
  hyperparameters <- model@parameters
  x <- hyperparameters$x
  hyperparameters$monotone_constraints <- NULL
  hyperparameters$y <- y
  hyperparameters$x <- x
  return(hyperparameters)
}

process_tbi_data <- function(tbi_path) {
  columns <- c("trial","age","hypoxia","hypotens","cisterns","shift","tsah","edh","d.pupil","d.motor","ctclass","d.sysbpt","hb","glucose","d.unfav","d.mort")
  columns_dict <- c("Trial","Age","Hypoxia","Hypotension","Compressed Cistern","Midline Shift","Traumatic Subarachnoid Hemorrhage", "Epidural Hematoma", "Pupils", "Motor", "CT brain classification","SBP","Hemoglobin","Blood glucose","Unfavarable outcomes", "Death")
  tbi <- foreign::read.spss(tbi_path, to.data.frame=T)[,columns]
  tbi <- tbi %>% mutate(trial=case_when(str_detect(trial, "US") ~ "US", 
                                        str_detect(trial, "Inter") ~ "International")) 
  
  tbi$ctclass_gr <- with(tbi, factor(case_when(
    ctclass %in% c(1,2) ~ "Grade 1-2",
    ctclass == 3 ~ "Grade 3",
    ctclass == 4 ~ "Grade 4",
    ctclass %in% c(5,6) ~ "Grade 5-6"
  )))
  tbi$d.motor_gr <- with(tbi, relevel(factor(case_when(
    d.motor %in% c(1,2) ~ "Grade 1-2",
    d.motor == 3 ~ "Grade 3",
    d.motor == 4 ~ "Grade 4",
    d.motor %in% c(5,6) ~ "Grade 5-6"
  )), ref="Grade 5-6"))
  
  
  tbi$d.mort.factor <- factor(tbi$d.mort)
  tbi$d.unfav.factor <- factor(tbi$d.unfav)
  
  return(tbi)
}

original_surv <- function(data, survmodel="cox") {
  S <- Surv(data$TEVENT, data$y)
  if (survmodel == "cox") {
    f <- coxph(S ~ AGE2+SEX+factor(SMOKING)+factor(alcohol)+BMIO+SYSTH+HDLO+DIABETES+HISTCAR2+HOMOCO+log(CREATO)+factor(albumin)+STENOSIS+IMTO, data=data)
  } else if (survmodel == "exp") {
    f <- survreg(S ~ AGE2+SEX+factor(SMOKING)+factor(alcohol)+BMIO+SYSTH+HDLO+DIABETES+HISTCAR2+HOMOCO+log(CREATO)+factor(albumin)+STENOSIS+IMTO, data=data, dist="exponential")
  }
  return(f)
}