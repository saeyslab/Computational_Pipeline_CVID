library(tidyverse)
library(randomForest)
library(readxl)
library(e1071)
library(qsub)

# ---------------------------- Class to use -----------------------------------

HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency")

group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)

seed <- 1


# ----------------- load in meta-data -----------------------------------------                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ---------- Read meta-data ----------------------------------

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



# ------------------------------- RF with FS ----------------------------------


# ------------------------- read in correct matrix ----------------------------

tube <- 1

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
                 "_FlowSOM_WithFS_RF.txt")
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

results_RF_3 <- list()
results_RF_2 <- list()
selected_featureslist_2 <- list()
selected_featureslist_3 <- list()

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
  
  # three class model (Kruskal Wallis tests)
  
  handle_3 <- qsub_lapply(
    X = c(5:(dim(z_scored_train_data)[2])),
    qsub_config = override_qsub_config(
      name = "Kruskal_Wallis_test", 
      verbose = F,
      wait = T,
      memory = "4G"
    ),
    qsub_environment = c("z_scored_train_data"),
    FUN = function(x){
      set.seed(1)
      library(tidyverse)
      train_data_new <-  z_scored_train_data %>% select(x,'Diagnosis')
      colnames(train_data_new) <- gsub("-", "neg", colnames(train_data_new))
      colnames(train_data_new) <- gsub("\\+", "pos", colnames(train_data_new))
      colnames(train_data_new) <- gsub("/", "_", colnames(train_data_new))
      colnames(train_data_new) <- gsub("%", "Perc_", colnames(train_data_new))
      colnames(train_data_new) <- gsub("\\([^)]*\\)", "", colnames(train_data_new))
      colnames(train_data_new) <- gsub(" ", "", colnames(train_data_new))
      name <- colnames(train_data_new)[1]
      test <- kruskal.test(formula(paste0(name ,"~Diagnosis")), data = train_data_new )
      
      
      p_value <- test$p.value
      names(p_value) <- colnames(z_scored_train_data)[x]
      return(p_value)
    }  )
  
  
  # Correlation check to remove correlated features
  
  na_removed_3 <- unlist(handle_3[!is.na(handle_3)])
  adj_p_3 <- stats::p.adjust(na_removed_3, "BH")
  sorted_3 <- sort(adj_p_3)
  z_scored_train_data <- z_scored_train_data %>% select("SAMPLE ID", "Diagnosis", "Diagnosis_2",  names(sorted_3))
  selected_features_3 <- c(4,5)
  for(i_3 in 6:length(sorted_3)){
    max_3 <- 0
    for(f_3 in selected_features_3){
      if(sd(z_scored_train_data %>% .[[f_3]]) > 0 & sd(z_scored_train_data %>% .[[i_3]]) >0 ){
        tmp_3 <- abs(cor(z_scored_train_data %>% .[[i_3]],z_scored_train_data %>% .[[f_3]]))
        if(tmp_3 > max_3) max_3 <- tmp_3
      } else {
        max_3 <- 1
      }
    }
    if(max_3 < 0.2){
      selected_features_3 <- c(selected_features_3,i_3)
      print(selected_features_3)
    } 
  }
  
  
  selected_featurenames_3 <- colnames(z_scored_train_data)[selected_features_3]
  selected_featureslist_3[[paste0("Experiment",experiment)]] <- selected_featurenames_3
  
  selected_train_data_3 <- z_scored_train_data %>% select("SAMPLE ID", "Diagnosis", "Diagnosis_2", selected_featurenames_3)
  selected_test_data_3 <- z_scored_test_data %>% select("SAMPLE ID", "Diagnosis", "Diagnosis_2",selected_featurenames_3)    
  
  
  
  res3 <- randomForest(x = selected_train_data_3[,-c(1:4)], 
                       y = selected_train_data_3$Diagnosis)
  
  res3$prediction <- predict(res3, selected_test_data_3[,-c(1:4)])
  res3$testset <- selected_test_data_3$`SAMPLE ID`
  res3$testsetpheno <- selected_test_data_3$Diagnosis
  results_RF_3[[paste0("Experiment",experiment)]] <- res3
  
  
  # For the two classes (Wilcoxon tests)
  handle_2 <- qsub_lapply(
    X = c(5:dim(z_scored_train_data)[2]),
    qsub_config = override_qsub_config(
      name = "Wilcox_test", 
      verbose = F,
      wait = T,
      memory = "4G",
      remove_tmp_folder = TRUE
    ),
    qsub_environment = c("z_scored_train_data"),
    FUN = function(x){
      set.seed(1)
      library(tidyverse)
      warning(sessionInfo())
      test <- wilcox.test(z_scored_train_data %>% tbl_df() %>% dplyr::filter(Diagnosis_2 == "CVID") %>% 
                            .[[x]],
                          z_scored_train_data %>% tbl_df() %>% dplyr::filter(Diagnosis_2 == "No_CVID") %>% 
                            .[[x]])
      p_value <- test$p.value
      names(p_value) <- colnames(z_scored_train_data)[x]
      return(p_value)
    } )
  
  
  # Remove correlated features
  na_removed_2 <- unlist(handle_2[!is.na(handle_2)])
  na_removed_2 <- na.omit(na_removed_2)
  adj_p_2 <- stats::p.adjust(na_removed_2, "BH")
  sorted_2 <- sort(adj_p_2)
  z_scored_train_data <- z_scored_train_data %>% select("SAMPLE ID", "Diagnosis", "Diagnosis_2", names(sorted_2))
  selected_features_2 <- c(4,5)
  for(i_2 in 6:length(sorted_2)){
    max_2 <- 0
    for(f_2 in selected_features_2){
      if(sd(z_scored_train_data %>% .[[f_2]]) > 0 & sd(z_scored_train_data %>% .[[i_2]]) >0 ){
        tmp_2 <- abs(cor(z_scored_train_data %>% .[[i_2]],z_scored_train_data %>% .[[f_2]]))
        if(tmp_2 > max_2) max_2 <- tmp_2
      } else {
        max_2 <- 1
      }
    }
    if(max_2 < 0.2){
      selected_features_2 <- c(selected_features_2,i_2)
      print(selected_features_2)
    } 
  }
  
  selected_featurenames_2 <- colnames(z_scored_train_data)[selected_features_2]
  selected_featureslist_2[[paste0("Experiment",experiment)]] <- selected_featurenames_2
  
  selected_train_data_2 <- z_scored_train_data %>% select("SAMPLE ID", "Diagnosis", "Diagnosis_2", selected_featurenames_2)
  selected_test_data_2 <- z_scored_test_data %>% select("SAMPLE ID", "Diagnosis", "Diagnosis_2", selected_featurenames_2)    
  
  
  res2 <- randomForest(x = selected_train_data_2[,-c(1:4)], 
                       y = as.factor(selected_train_data_2$Diagnosis_2))
  
  res2$prediction <- predict(res2, selected_test_data_2[,-c(1:4)])
  res2$testset <- selected_test_data_2$`SAMPLE ID`
  res2$testsetpheno <- selected_test_data_2$Diagnosis_2
  results_RF_2[[paste0("Experiment",experiment)]] <- res2
  
  
  
  wrongly_predicted_3 <- which(res3$prediction != selected_test_data_3$Diagnosis)
  
  if (length(wrongly_predicted_3) != 0){
    wrongly_3 <- data.frame("Patient" = selected_test_data_3$`SAMPLE ID`[wrongly_predicted_3],
                            "Diagnosis" = selected_test_data_3$Diagnosis[wrongly_predicted_3],
                            "Wrong_diagnosis" = res3$prediction[wrongly_predicted_3])
    
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for three classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output(print(wrongly_3, print.gap=3), append = TRUE, file=logFile)
  }
  
  wrongly_predicted_2 <- which(res2$prediction != selected_test_data_2$Diagnosis_2)
  
  if (length(wrongly_predicted_2) != 0){
    wrongly_2 <- data.frame("Patient" = selected_test_data_2$`SAMPLE ID`[wrongly_predicted_2],
                            "Diagnosis" = selected_test_data_2$Diagnosis_2[wrongly_predicted_2],
                            "Wrong_diagnosis" = res2$prediction[wrongly_predicted_2])
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for two classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output( print(wrongly_2, print.gap=3), append = TRUE, file=logFile)
    cat("\n", file = logFile, append = TRUE)
    
  }
}



save(results_RF_2, results_RF_3, 
     selected_featureslist_2, selected_featureslist_3, 
     file = paste0("Manual_gating/Data/Randomforests_With_FS_Manual_Gating_tube_",tube , ".Rdata"))


#load(file = paste0("Manual_gating/Data/Randomforests_With_FS_Manual_Gating_tube_",tube , ".Rdata"))

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




# ------------------------------- SVM model -----------------------------------
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
  filter(!(`SAMPLE ID` %in% c("PIDHC089")))


logFile = paste0("Manual_gating/LogFiles/log_file_for_manual_gating_tube",tube,
                 "_WihFS_SVM.txt")
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
  
  # Read in selected features from feature selection 
  load(paste0("Manual_gating/Data/Randomforests_With_FS_Manual_Gating_tube_",tube , ".Rdata"))
  
  selected_featurenames_2 <- selected_featureslist_2[[(experiment-12)]]
  selected_featurenames_3 <- selected_featureslist_3[[(experiment-12)]]
  selected_train_data_2 <- z_scored_train_data %>% select( "Diagnosis_2", "Diagnosis",  "SAMPLE ID", selected_featurenames_2)
  selected_test_data_2 <- z_scored_test_data %>% select(  "Diagnosis_2","Diagnosis", "SAMPLE ID", selected_featurenames_2)    
  selected_train_data_3 <- z_scored_train_data %>% select("Diagnosis_2", "Diagnosis",  "SAMPLE ID", selected_featurenames_3)
  selected_test_data_3 <- z_scored_test_data %>% select("Diagnosis_2", "Diagnosis",  "SAMPLE ID", selected_featurenames_3)    
  
  
  
  # For the three variables
  res3 <- svm(Diagnosis ~ ., data = selected_train_data_3[,-c(1,3,4)], kernel = "linear", cost =1, scale = TRUE)
  
  
  res3$prediction <- predict(res3, selected_test_data_3[,-c(1:4)])
  res3$testset <- selected_test_data_3$`SAMPLE ID`
  res3$testsetpheno <- selected_test_data_3$Diagnosis
  results_SVM_3[[paste0("Experiment",experiment)]] <- res3
  
  
  # For the two variables
  res2 <- svm(Diagnosis_2 ~ ., data = selected_train_data_2[,-c(2,3,4)], kernel = "linear", cost =1, scale = TRUE)
  
  
  res2$prediction <- predict(res2, selected_test_data_2[,-c(1:4)])
  res2$testset <-selected_test_data_2$`SAMPLE ID`
  res2$testsetpheno <- selected_test_data_2$Diagnosis_2
  results_SVM_2[[paste0("Experiment",experiment)]] <- res2
  
  wrongly_predicted_3 <- which(res3$prediction != selected_test_data_3$Diagnosis)
  
  if (length(wrongly_predicted_3) != 0){
    wrongly_3 <- data.frame("Patient" = selected_test_data_3$`SAMPLE ID`[wrongly_predicted_3],
                            "Diagnosis" = selected_test_data_3$Diagnosis[wrongly_predicted_3],
                            "Wrong_diagnosis" = res3$prediction[wrongly_predicted_3])
    
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for three classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output(print(wrongly_3, print.gap=3), append = TRUE, file=logFile)
  }
  
  wrongly_predicted_2 <- which(res2$prediction != selected_test_data_2$Diagnosis_2)
  
  if (length(wrongly_predicted_2) != 0){
    wrongly_2 <- data.frame("Patient" = selected_test_data_2$`SAMPLE ID`[wrongly_predicted_2],
                            "Diagnosis" = selected_test_data_2$Diagnosis_2[wrongly_predicted_2],
                            "Wrong_diagnosis" = res2$prediction[wrongly_predicted_2])
    cat(paste0("For experiment day ", experiment, " patients wrongly predicted for two classes: "), file = logFile, append = TRUE, sep = "\n")
    capture.output( print(wrongly_2, print.gap=3), append = TRUE, file=logFile)
    cat("\n", file = logFile, append = TRUE)
    
  }
}



save(results_SVM_2, results_SVM_3, file = paste0("Manual_gating/Data/SVM_WithFS_Manual_Gating_tube",tube , ".Rdata"))


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


