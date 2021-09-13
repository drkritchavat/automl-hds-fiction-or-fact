source("report/functions.R")
library(h2o)
h2o.init(ice_root="/home/kp/tmp_h2o")

# load and process tbi_data
tbi <- process_tbi_data("/home/kp/hds/Dissertation/Datasets/cpm/tbi/TBI.sav")
# set label
y <- list(mort="d.mort.factor",unfav="d.unfav.factor")

# set model path
path_model_mort <- "/home/kp/hds/Dissertation/AutoML/AutoML_tbi_mort@@d_mort"
path_model_unfav <- "/home/kp/hds/Dissertation/AutoML/AutoML_tbi_mort_2@@d_unfav"
# loading leaderboard
leaderboard_mort <- read.csv(paste(path_model_mort, "leaderboard.csv", sep="/"))
leaderboard_unfav <- read.csv(paste(path_model_unfav, "leaderboard.csv", sep="/"))


h2o_model_mort <- get_model_from_leaderboard(leaderboard_mort, path_model_mort, 1, y$mort)
h2o_model_unfav <- get_model_from_leaderboard(leaderboard_unfav, path_model_unfav, 1, y$unfav)

# 5 fold CV
k <- 5
set.seed(123)
cv_folds <- createFolds(1:dim(tbi)[1], k = k, list = TRUE, returnTrain = FALSE)
set.seed(123)
cv_folds_inter <- createFolds(1:sum(tbi$trial == "International"), k = k, list = TRUE, returnTrain = FALSE)
set.seed(123)
cv_folds_us <- createFolds(1:sum(tbi$trial == "US"), k = k, list = TRUE, returnTrain = FALSE)

# Pooled
results_h2o <- cv_h2o(tbi, h2o_model_mort, h2o_model_unfav, cv_folds)
# Subgroups
results_h2o_inter <- cv_h2o(tbi[tbi$trial=="International", ], h2o_model_mort, h2o_model_unfav, cv_folds_inter)
results_h2o_us <- cv_h2o(tbi[tbi$trial=="US", ], h2o_model_mort, h2o_model_unfav, cv_folds_us)

# cv for 2nd+ models
results_others <- list()
for (i in 2:10) {
  h2o_model_mort <- get_model_from_leaderboard(leaderboard_mort, path_model_mort, i, y$mort)
  h2o_model_unfav <- get_model_from_leaderboard(leaderboard_unfav, path_model_unfav, 1, y$unfav)
  results_others[[i]] <- cv_h2o(tbi, h2o_model_mort, h2o_model_unfav, cv_folds)
  print(i)
}

# save cv results
saveRDS(results_h2o, "/home/kp/hds/Dissertation/report/h2o_results.Rdata")

saveRDS(results_h2o_inter, "/home/kp/hds/Dissertation/report/h2o_results_inter.Rdata")
saveRDS(results_h2o_us, "/home/kp/hds/Dissertation/report/h2o_results_us.Rdata")

saveRDS(results_others, "/home/kp/hds/Dissertation/report/h2o_results_others.Rdata")

##################################
# sensitivity
##################################

best_mort <- h2o.getModel("GBM_grid__1_AutoML_20210903_011310_model_21")@parameters
best_unfav <- h2o.getModel("DRF_1_AutoML_20210903_012128")@parameters
best_mort$y <- "d.mort.factor"
best_mort$y <- "d.unfav.factor"
best_mort$monotone_constraints <- NULL
best_unfav$monotone_constraints <- NULL

k <- length(cv_folds)
aucs <- matrix(nrow=k, ncol=2)
colnames(aucs) <- c("Death", "Unfavorable outcomes")
calibrations_mort_h2o <- c()
calibrations_unfav_h2o <- c()
calibrations_glm_mort_h2o <- matrix(ncol=2, nrow=k)
calibrations_glm_unfav_h2o <- matrix(ncol=2, nrow=k)
briers <- matrix(ncol=2, nrow=k)
imp.tbi %>%
  mutate(d.unfav.factor=factor(d.unfav), d.mort.factor=factor(d.mort)) %>%
  select(.id,d.mort,d.unfav,d.mort.factor,d.unfav.factor,age,d.motor_gr,d.pupil,hypoxia,hypotens,ctclass_gr,tsah) %>%
  mutate(sq.age=age^2) -> data
for (i in 1:k) {
  test_index <- cv_folds[[i]]
  train <- as.h2o(data[-(data$.id %in% test_index), ])
  test <- data[data$.id %in% test_index, ]
  best_mort$training_frame <- train
  best_unfav$training_frame <- train
  fit_mort_h2o <- do.call(get_h2o_from_parameters(best_mort), best_mort)
  fit_unfav_h2o <- do.call(get_h2o_from_parameters(best_unfav), best_unfav)
  mort_predict_h2o_df <- h2o.predict(fit_mort_h2o, as.h2o(test))
  mort_predict_h2o <- as.data.frame(mort_predict_h2o_df)[,"p1"]
  unfav_predict_h2o_df <- h2o.predict(fit_unfav_h2o, as.h2o(test))
  unfav_predict_h2o <- as.data.frame(unfav_predict_h2o_df)[,"predict"]
  
  data.frame(observed=test$d.mort, predicted=mort_predict_h2o) %>%
    mutate(cut=cut(predicted,10)) %>%
    group_by(cut) %>% summarise(observed=mean(observed), predicted=mean(predicted)) %>%
    mutate(fold=factor(i), ith=1:n()) %>%
    rbind(calibrations_mort_h2o) -> calibrations_mort_h2o
  
  data.frame(observed=test$d.mort, predicted=mort_predict_h2o) %>%
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


h2o.shutdown(prompt = FALSE)
