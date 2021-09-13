library(mice)
library(tidyverse)

impute_cols <- c("trial","age","hypoxia","hypotens","cisterns","shift","tsah","edh","d.pupil","d.motor","ctclass","d.sysbpt","hb","glucose","d.unfav","d.mort")
impute_cols_dict <- c("Trial","Age","Hypoxia","Hypotension","Compressed Cistern","Midline Shift","Traumatic Subarachnoid Hemorrhage", "Epidural Hematoma", "Pupils", "Motor", "CT brain classification","SBP","Hemoglobin","Blood glucose","Unfavarable outcomes", "Death")
tbi <- foreign::read.spss("/home/kp/hds/Dissertation/Datasets/cpm/tbi/TBI.sav", to.data.frame=T)[,impute_cols]
tbi <- tbi %>% mutate(trial=case_when(str_detect(trial, "US") ~ "US", 
                                      str_detect(trial, "Inter") ~ "International")) 
N <- dim(tbi)[1]

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

imp <- mice(tbi, m=10)
# save imputed object

saveRDS(imp, "/home/kp/hds/Dissertation/report/tbi_imputed.Rdata")


#library(tidyverse)
#imp.tbi %>%
# select(d.mort,d.unfav,age,d.motor_gr,d.pupil,hypoxia,hypotens,ctclass_gr,tsah) %>%
#  mutate(sq.age=age^2) %>%
#  write.csv( "/home/kp/hds/Dissertation/report/tbi_imputed.csv", row.names = F)

#library("h2o")
#h2o.init(ice_root="/home/kp/tmp_h2o")
