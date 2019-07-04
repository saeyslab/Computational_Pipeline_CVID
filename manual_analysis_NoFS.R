library(tidyverse)
library(randomForest)
library(readxl)
library(e1071)



# ---------------------------- Class to use -----------------------------------

HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency")

group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)

seed <- 1


# --------------------------- Read meta-data ----------------------------------

meta_data_file <- "../../PhD/CVID/Clinical_features_data.txt"
meta_data <- read.delim(meta_data_file, header = TRUE, sep = "\t", 
                        check.names = FALSE) %>% 
  dplyr::filter(diagnosis_group %in% group_to_use) %>%
  dplyr::filter(SAMPLE_ID != "PIDHC081")

if(length(group_to_use) == length(PID_HC)){
  meta_data$diagnosis_group <- as.character(meta_data$diagnosis_group)
  meta_data[!(meta_data$diagnosis_group %in% HC_CVID),"diagnosis_group"] <- "PAD"
  meta_data$diagnosis_group <- as.factor(meta_data$diagnosis_group)
}

# PIDHC081 is removed because he has too few cells for tube 2
rownames(meta_data) <- meta_data$SAMPLE_ID



# ------------------------------ Random Forests -------------------------------


# ------------------------- read in correct matrix ----------------------------

tube <- 2

if (tube != "Full"){
  giant_matrix <- readRDS(file = paste0("Manual_gating/Data/tube",tube,".rds"))
  cols_not_to_convert <- c(1,2,3)
  giant_matrix[,-c(1,2,3)] <- apply(giant_matrix[,-c(1,2,3)], 2, type.convert)
} else {
  giant_matrix <- readRDS(file = "Manual_gating/Data/tube1til5.rds")
  giant_matrix[,-c(1,2,3)] <- apply(giant_matrix[,-c(1,2,3)], 2, type.convert)
  cols_not_to_convert <- c(1,2,3)
}




# Leave out next patients because they are missing tube 3  values
giant_matrix <- giant_matrix %>% 
  filter(!(`SAMPLE ID` %in% c("PID057", "PID217", "PIDHC089")))


logFile = paste0("Manual_gating/LogFiles/log_file_for_manual_gating_tube",tube,
                 "_FlowSOM_NoFS_RF.txt")
cat(paste0("This is a log file for the patients that were predicted wrong from the RF without a certain experiment day"), file=logFile, append=FALSE, sep = "\n")


giant_matrix$Diagnosis_2 <- as.character(giant_matrix$Diagnosis)

giant_matrix$Diagnosis_2[which(giant_matrix$Diagnosis == "PID")] <- "No_CVID"
giant_matrix$Diagnosis_2[which(giant_matrix$Diagnosis == "HC")] <- "No_CVID"

giant_matrix$Diagnosis_2 <- as.factor(giant_matrix$Diagnosis_2)

giant_matrix <- giant_matrix %>% select(Diagnosis_2, everything())

cols_not_to_convert <- c(1,2,3,4)


meta_data <- meta_data %>% filter(`SAMPLE_ID` %in% giant_matrix$`SAMPLE ID`)
experiment <- 13
set.seed(seed)

results_RF_2 <- list()
results_RF_3 <- list()

source("R_Scripts/Z_score.R")

for(experiment in c(13:33)){
  set.seed(1)
  
  #print(paste0("Experiment day ", experiment))
  allfiles <- list.files(paste0("../../PhD/CVID/CVID_Delfien_orig/Experiment CVID_",
                                experiment), 
                         pattern = ".*PID.*.fcs$", recursive = TRUE)
  
  train_ID <- gsub('.*(PID.*)_.*_.*','\\1',allfiles) %>% unique
  test_data <- giant_matrix %>% dplyr::filter(`SAMPLE ID` %in% train_ID)
  train_data <- giant_matrix %>% dplyr::filter(!(`SAMPLE ID` %in% test_data$`SAMPLE ID`))
  
  
  z_scored_train_data_full <- suppressWarnings(UseZScore(to_be_scored = train_data, 
                                                         cols_not_to_convert,
                                                         meta_data, 
                                                         verbose = FALSE, 
                                                         start = 5, 
                                                         stop = dim(train_data)[2]))
  z_scored_train_data <- z_scored_train_data_full$z_score_matrix
  
  
  z_scored_test_data <- suppressWarnings(UseZScoreTestData(test_data, 
                                                           cols_not_to_convert, 
                                                           meta_data, 
                                                           z_scored_train_data_full$HC_scores,
                                                           start = 5, 
                                                           stop = dim(test_data)[2]))
  
  
  
  # For the three variables
  res3 <- randomForest(x = z_scored_train_data[,-c(1:4)], 
                       y = z_scored_train_data$Diagnosis)
  
  res3$prediction <- predict(res3, z_scored_test_data[,-c(1:4)])
  res3$testset <- z_scored_test_data$`SAMPLE ID`
  res3$testsetpheno <- z_scored_test_data$Diagnosis
  results_RF_3[[paste0("Experiment",experiment)]] <- res3
  
  
  # For the two variables
  res2 <- randomForest(x = z_scored_train_data[,-c(1:4)], 
                       y = as.factor(z_scored_train_data$Diagnosis_2))
  
  res2$prediction <- predict(res2, z_scored_test_data[,-c(1:4)])
  res2$testset <- z_scored_test_data$`SAMPLE ID`
  res2$testsetpheno <- z_scored_test_data$Diagnosis_2
  results_RF_2[[paste0("Experiment",experiment)]] <- res2
  
  wrongly_predicted_3 <- which(res3$prediction != z_scored_test_data$Diagnosis)
  
  if (length(wrongly_predicted_3) != 0){
    wrongly_3 <- data.frame("Patient" = z_scored_test_data$`SAMPLE ID`[wrongly_predicted_3],
                            "Diagnosis" = z_scored_test_data$Diagnosis[wrongly_predicted_3],
                            "Wrong_diagnosis" = res3$prediction[wrongly_predicted_3])
    
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for three classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output(print(wrongly_3, print.gap=3), append = TRUE, file=logFile)
  }
  
  wrongly_predicted_2 <- which(res2$prediction != z_scored_test_data$Diagnosis_2)
  
  if (length(wrongly_predicted_2) != 0){
    wrongly_2 <- data.frame("Patient" = z_scored_test_data$`SAMPLE ID`[wrongly_predicted_2],
                            "Diagnosis" = z_scored_test_data$Diagnosis_2[wrongly_predicted_2],
                            "Wrong_diagnosis" = res2$prediction[wrongly_predicted_2])
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for two classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output( print(wrongly_2, print.gap=3), append = TRUE, file=logFile)
    cat("\n", file = logFile, append = TRUE)
    
  }
}



save(results_RF_2, results_RF_3, file = paste0("Manual_gating/Data/Randomforests_Manual_Gating_tube",tube , ".Rdata"))


final_2 <- matrix(0, nrow = 2, ncol = 2)
for (experiment in results_RF_2){
  t_2 <- table(experiment$prediction, experiment$testsetpheno)
  if (("CVID" %in% colnames(t_2) == FALSE)){
    t_2 <- cbind("CVID"=c(0,0), t_2)
  }
  if (("No_CVID" %in% colnames(t_2) == FALSE)){
    t_2 <- cbind("CVID"=c(0,0), t_2)
  }
  final_2 <- final_2 + t_2
}
final_3 <- matrix(0, nrow = 3, ncol = 3)
for (experiment in results_RF_3){
  t_3 <- table(experiment$prediction, experiment$testsetpheno)
  final_3 <- final_3 + t_3
}

cat("\n", file = logFile, append = TRUE)
cat(paste0("Result for three classes: "), file = logFile, append = TRUE, sep = "\n")
capture.output( print(final_3, print.gap=3), append = TRUE, file=logFile)


cat("\n", file = logFile, append = TRUE)
cat(paste0("Result for two classes: "), file = logFile, append = TRUE, sep = "\n")

capture.output( print(final_2, print.gap=3), append = TRUE, file=logFile)



# ------------------------- Support vector machine ----------------------------

# ------------------------- read in correct matrix ----------------------------

tube <- 2

if (tube != "Full"){
  giant_matrix <- readRDS(file = paste0("Manual_gating/Data/tube",tube,".rds"))
  cols_not_to_convert <- c(1,2,3)
  giant_matrix[,-c(1,2,3)] <- apply(giant_matrix[,-c(1,2,3)], 2, type.convert)
} else {
  giant_matrix <- readRDS(file = "Manual_gating/Data/tube1til5.rds")
  giant_matrix[,-c(1,2,3)] <- apply(giant_matrix[,-c(1,2,3)], 2, type.convert)
  cols_not_to_convert <- c(1,2,3)
}




# Leave out next patients because they are missing tube 3  values
giant_matrix <- giant_matrix %>% 
  filter(!(`SAMPLE ID` %in% c("PID057", "PID217", "PIDHC089")))


logFile = paste0("Manual_gating/LogFiles/log_file_for_manual_gating_tube",tube,
                 "_FlowSOM_NoFS_SVM.txt")
cat(paste0("This is a log file for the patients that were predicted wrong from the RF without a certain experiment day"), file=logFile, append=FALSE, sep = "\n")


giant_matrix$Diagnosis_2 <- as.character(giant_matrix$Diagnosis)

giant_matrix$Diagnosis_2[which(giant_matrix$Diagnosis == "PID")] <- "No_CVID"
giant_matrix$Diagnosis_2[which(giant_matrix$Diagnosis == "HC")] <- "No_CVID"

giant_matrix$Diagnosis_2 <- as.factor(giant_matrix$Diagnosis_2)

giant_matrix <- giant_matrix %>% select(Diagnosis_2, everything())

cols_not_to_convert <- c(1,2,3,4)


meta_data <- meta_data %>% filter(`SAMPLE_ID` %in% giant_matrix$`SAMPLE ID`)
experiment <- 13
set.seed(seed)

results_SVM_2 <- list()
results_SVM_3 <- list()

source("R_Scripts/Z_score.R")

experiment <- 13
for(experiment in c(13:33)){
  set.seed(1)
  
  #print(paste0("Experiment day ", experiment))
  allfiles <- list.files(paste0("../../PhD/CVID/CVID_Delfien_orig/Experiment CVID_",
                                experiment), 
                         pattern = ".*PID.*.fcs$", recursive = TRUE)
  
  train_ID <- gsub('.*(PID.*)_.*_.*','\\1',allfiles) %>% unique
  test_data <- giant_matrix %>% dplyr::filter(`SAMPLE ID` %in% train_ID)
  train_data <- giant_matrix %>% dplyr::filter(!(`SAMPLE ID` %in% test_data$`SAMPLE ID`))
  
  
  z_scored_train_data_full <- suppressWarnings(UseZScore(to_be_scored = train_data, 
                                                         cols_not_to_convert,
                                                         meta_data, 
                                                         verbose = FALSE, 
                                                         start = 5, 
                                                         stop = dim(train_data)[2]))
  z_scored_train_data <- z_scored_train_data_full$z_score_matrix
  
  
  z_scored_test_data <- suppressWarnings(UseZScoreTestData(test_data, 
                                                           cols_not_to_convert, 
                                                           meta_data, 
                                                           z_scored_train_data_full$HC_scores,
                                                           start = 5, 
                                                           stop = dim(test_data)[2]))
  
  
  
  # For the three variables
  res3 <- svm(Diagnosis ~ ., data = z_scored_train_data[,-c(1,3,4)], kernel = "linear", cost =1, scale = TRUE)
  
  
  res3$prediction <- predict(res3, z_scored_test_data[,-c(1:4)])
  res3$testset <- z_scored_test_data$`SAMPLE ID`
  res3$testsetpheno <- z_scored_test_data$Diagnosis
  results_SVM_3[[paste0("Experiment",experiment)]] <- res3
  
  
  # For the two variables
  res2 <- svm(Diagnosis_2 ~ ., data = z_scored_train_data[,-c(2,3,4)], kernel = "linear", cost =1, scale = TRUE)
  
  
  res2$prediction <- predict(res2, z_scored_test_data[,-c(1:4)])
  res2$testset <- z_scored_test_data$`SAMPLE ID`
  res2$testsetpheno <- z_scored_test_data$Diagnosis_2
  results_SVM_2[[paste0("Experiment",experiment)]] <- res2
  
  wrongly_predicted_3 <- which(res3$prediction != z_scored_test_data$Diagnosis)
  
  if (length(wrongly_predicted_3) != 0){
    wrongly_3 <- data.frame("Patient" = z_scored_test_data$`SAMPLE ID`[wrongly_predicted_3],
                            "Diagnosis" = z_scored_test_data$Diagnosis[wrongly_predicted_3],
                            "Wrong_diagnosis" = res3$prediction[wrongly_predicted_3])
    
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for three classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output(print(wrongly_3, print.gap=3), append = TRUE, file=logFile)
  }
  
  wrongly_predicted_2 <- which(res2$prediction != z_scored_test_data$Diagnosis_2)
  
  if (length(wrongly_predicted_2) != 0){
    wrongly_2 <- data.frame("Patient" = z_scored_test_data$`SAMPLE ID`[wrongly_predicted_2],
                            "Diagnosis" = z_scored_test_data$Diagnosis_2[wrongly_predicted_2],
                            "Wrong_diagnosis" = res2$prediction[wrongly_predicted_2])
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for two classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output( print(wrongly_2, print.gap=3), append = TRUE, file=logFile)
    cat("\n", file = logFile, append = TRUE)
    
  }
}



save(results_SVM_2, results_SVM_3, file = paste0("Manual_gating/Data/SVM_Manual_Gating_tube",tube , ".Rdata"))


final_2 <- matrix(0, nrow = 2, ncol = 2)
for (experiment in results_SVM_2){
  t_2 <- table(experiment$prediction, experiment$testsetpheno)
  if (("CVID" %in% colnames(t_2) == FALSE)){
    t_2 <- cbind("CVID"=c(0,0), t_2)
  }
  if (("No_CVID" %in% colnames(t_2) == FALSE)){
    t_2 <- cbind("CVID"=c(0,0), t_2)
  }
  final_2 <- final_2 + t_2
}
final_3 <- matrix(0, nrow = 3, ncol = 3)
for (experiment in results_SVM_3){
  t_3 <- table(experiment$prediction, experiment$testsetpheno)
  final_3 <- final_3 + t_3
}

cat("\n", file = logFile, append = TRUE)
cat(paste0("Result for three classes: "), file = logFile, append = TRUE, sep = "\n")
capture.output( print(final_3, print.gap=3), append = TRUE, file=logFile)


cat("\n", file = logFile, append = TRUE)
cat(paste0("Result for two classes: "), file = logFile, append = TRUE, sep = "\n")
capture.output( print(final_2, print.gap=3), append = TRUE, file=logFile)


