library(tidyverse)
library(e1071)
library(crayon)

tube <- "1_2_5" #(1 = PBMCs, 2 = B cells, 5 = T cells, 1_2_5 = Combination of all tubes)
seed <- 1
kernel <- "linear"


# --------------- With what group will we be working on? ----------------------
HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency")

group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID_New_preprocessing" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)


# ------ List directories and the fsom objects with the feature matrices ------
full_directory <- paste0("FlowSOM/Tube",tube, "/Data",PID,"/Crossvalidation")
working_directory <- file.path(full_directory, paste0("SVM_", kernel, "_noFS"))
Logfile_directory <- file.path(full_directory, "Logfiles")
matrices_directory <- file.path(full_directory, "Matrices")

suppressWarnings(dir.create(working_directory))
suppressWarnings(dir.create(Logfile_directory))


# --------------------------- Open meta_data----------------------------------
meta_data_file <- "../../PhD/CVID/Clinical_features_data.txt"
meta_data <- read.delim(meta_data_file, header = TRUE, sep = "\t", 
                        check.names = FALSE) %>% 
  dplyr::filter(diagnosis_group %in% group_to_use) %>%
  dplyr::filter(SAMPLE_ID != "PIDHC081")

if(length(group_to_use) == length(PID_HC)){
  meta_data$diagnosis_group <- as.character(meta_data$diagnosis_group)
  meta_data[!(meta_data$diagnosis_group %in% HC_CVID),"diagnosis_group"] <- "PID"
  meta_data$diagnosis_group <- as.factor(meta_data$diagnosis_group)
}

# PIDHC081 is removed because he has too few cells for tube 2
rownames(meta_data) <- meta_data$SAMPLE_ID



# ---------------- Define the number of matrix combinations -------------------

all <- c("meta_clusters", "clusters", "percentages_clusters", "percentages_metaclusters","total", "percentages_to")
to_use_data <- all[4]

# ------------------------- z-score + RF---------------------------------------
#Load in data object
source("R_Scripts/Z_score.R")

for (to_use_data in all){
  
  
  
  #set up writing
  logFile = paste0(Logfile_directory,"/log_file_for_",
                   to_use_data,
                   "_FlowSOM_NoFS_SVM_", kernel, ".txt")
  cat(paste0("This is a log file for the patients that were predicted wrong from the RF without a certain experiment day where following features were used: ", to_use_data), file=logFile, append=FALSE, sep = "\n")
  
  
  results_SVM_3 <- list()
  results_SVM_2 <- list()
  
  
  experiment_day <- 33
  for (experiment_day in c(13:33)){
    set.seed(seed)
    print(paste0("Now calculating experiment ",experiment_day))
    
    # Make correct featureset
    
    if (to_use_data == "meta_clusters"){
      percentages_metaclusters <- readRDS(file.path(matrices_directory, paste0("percentages_metaclusters_experiment_day_",experiment_day, ".RDS")))
      MFI_metaclusters <- readRDS(file.path(matrices_directory, paste0("MFI_metaclusters_experiment_day_",experiment_day, ".RDS")))
      
      combined_data <- cbind("SAMPLE ID" = rownames(percentages_metaclusters), 
                             percentages_metaclusters, 
                             MFI_metaclusters) %>% 
        tbl_df %>% 
        arrange(`SAMPLE ID`)
    } else if (to_use_data == "clusters"){
      percentages_clusters <- readRDS(file.path(matrices_directory, paste0("percentages_clusters_experiment_day_",experiment_day, ".RDS")))
      MFI_clusters <- readRDS(file.path(matrices_directory, paste0("MFI_clusters_experiment_day_",experiment_day, ".RDS")))
      
      combined_data <- cbind("SAMPLE ID" = rownames(percentages_clusters), 
                             percentages_clusters, 
                             MFI_clusters) %>% 
        tbl_df %>% 
        arrange(`SAMPLE ID`)
    } else if (to_use_data == "total"){
      percentages_metaclusters <- readRDS(file.path(matrices_directory, paste0("percentages_metaclusters_experiment_day_",experiment_day, ".RDS")))
      MFI_metaclusters <- readRDS(file.path(matrices_directory, paste0("MFI_metaclusters_experiment_day_",experiment_day, ".RDS")))
      percentages_clusters <- readRDS(file.path(matrices_directory, paste0("percentages_clusters_experiment_day_",experiment_day, ".RDS")))
      MFI_clusters <- readRDS(file.path(matrices_directory, paste0("MFI_clusters_experiment_day_",experiment_day, ".RDS")))
      percentages_to_clusters_metaclusters <- readRDS(file.path(matrices_directory, paste0("percentages_to_clusters_metaclusters_experiment_day_",experiment_day, ".RDS")))
      
      colnames(percentages_clusters) <- paste(colnames(percentages_clusters), "c", sep = "_")
      colnames(percentages_metaclusters) <- paste(colnames(percentages_metaclusters), "m", sep = "_")
      colnames(MFI_clusters) <- paste(colnames(MFI_clusters), "c", sep = "_")
      colnames(MFI_metaclusters) <- paste(colnames(MFI_metaclusters), "m", sep = "_")
      colnames(percentages_to_clusters_metaclusters) <- paste(colnames(percentages_to_clusters_metaclusters), "c_to_mc", sep = "_")
      
      
      
      combined_data <- cbind("SAMPLE ID" = rownames(percentages_clusters), 
                             percentages_clusters, 
                             MFI_clusters,
                             percentages_metaclusters,
                             MFI_metaclusters,
                             percentages_to_clusters_metaclusters) %>% 
        tbl_df %>% 
        arrange(`SAMPLE ID`)
    } else if (to_use_data == "percentages_clusters"){
      percentages_clusters <- readRDS(file.path(matrices_directory, paste0("percentages_clusters_experiment_day_",experiment_day, ".RDS")))
      combined_data <- cbind("SAMPLE ID" = rownames(percentages_clusters), 
                             percentages_clusters) %>% 
        tbl_df %>% 
        arrange(`SAMPLE ID`)
    } else if (to_use_data == "percentages_metaclusters"){
      percentages_metaclusters <- readRDS(file.path(matrices_directory, paste0("percentages_metaclusters_experiment_day_",experiment_day, ".RDS")))
      combined_data <- cbind("SAMPLE ID" = rownames(percentages_metaclusters), 
                             percentages_metaclusters) %>% 
        tbl_df %>% 
        arrange(`SAMPLE ID`)
    } else if (to_use_data == "percentages_to"){
      percentages_to_clusters_metaclusters <- readRDS(file.path(matrices_directory, paste0("percentages_to_clusters_metaclusters_experiment_day_",experiment_day, ".RDS")))
      
      combined_data <- cbind("SAMPLE ID" = rownames(percentages_to_clusters_metaclusters), 
                             percentages_to_clusters_metaclusters) %>% 
        tbl_df %>% 
        arrange(`SAMPLE ID`)
    } 
    
    combined_data <- combined_data %>% dplyr::filter(`SAMPLE ID` %in% meta_data$SAMPLE_ID)
    
    combined_data[,-1] <- apply(combined_data[,-1], 2, type.convert)
    
    suppressWarnings(meta_data <- meta_data %>% 
                       dplyr::filter(`SAMPLE_ID` %in% 
                                       combined_data$`SAMPLE ID`) %>% 
                       arrange(`SAMPLE_ID`))
    
    
    combined_data <- combined_data %>% 
      .[!duplicated(.[,"SAMPLE ID"]),] %>%  
      mutate("diagnosis" = meta_data$diagnosis_group) 
    
    combined_data <- combined_data %>% mutate("diagnosis_2" = ifelse(combined_data$diagnosis == "CVID", "CVID", "No_CVID"))
    combined_data <- combined_data %>% mutate("diagnosis_2" = as.factor(diagnosis_2))
    
    to_compare_tube <- tube
    if (tube == "1_2" || tube == "1_2_5"){to_compare_tube <- 1}
    
    files <- list.files(paste0("Preprocessing/Tube",to_compare_tube,"/"), 
                        recursive = TRUE, 
                        pattern = ".*_pp.fcs")
    
    
    test_data_names <- files[grep(paste0("Exp_",experiment_day),files)] %>% 
      gsub('.*/(.*)_.*_.*_.*','\\1',.)
    
    test_data_names <- test_data_names[which(test_data_names %in% meta_data$SAMPLE_ID)]
    train_data_names <- combined_data$`SAMPLE ID`[-which(test_data_names %in% combined_data$`SAMPLE ID`)]  
    train_data_names <- train_data_names[which(train_data_names %in% meta_data$SAMPLE_ID)]
    
    test_data <- combined_data %>% dplyr::filter(`SAMPLE ID` %in% test_data_names)
    train_data <- combined_data %>% dplyr::filter(`SAMPLE ID` %in% train_data_names)
    
    dimension <- train_data %>% colnames %>% length
    
    train <- suppressWarnings(UseZScore(to_be_scored = train_data, 
                                        cols_not_to_convert = c(1,dimension -1, dimension), 
                                        meta_data = meta_data, 
                                        verbose = TRUE,
                                        start = 2,
                                        stop = dimension-2))
    
    z_scored_train_data <- train$z_score_matrix
    
    
    z_scored_test_data <- suppressWarnings(UseZScoreTestData(test_data, 
                                                             cols_not_to_convert = c(1, dimension - 1, dimension), 
                                                             meta_data, 
                                                             HC_scores = train$HC_scores,
                                                             start = 2, 
                                                             stop = dimension-2))
    
    # For the three classes
    res3 <- svm(diagnosis ~ ., data = z_scored_train_data[,-c(1,3)], kernel = "linear", cost =1, scale = TRUE)
    res3$prediction <- predict(res3, newdata = z_scored_test_data[,-c(1:3)])
    res3$testset <- z_scored_test_data$`SAMPLE ID`
    res3$testsetpheno <- z_scored_test_data$diagnosis
    results_SVM_3[[paste0("Experiment",experiment_day)]] <- res3
    
    # For the two classes
    res2 <- svm(diagnosis_2 ~ ., data = z_scored_train_data[,-c(1,2)], kernel = "linear", cost =1, scale = TRUE)
    
    res2$prediction <- predict(res2, z_scored_test_data[,-c(1:3)])
    res2$testset <- z_scored_test_data$`SAMPLE ID`
    res2$testsetpheno <- z_scored_test_data$diagnosis_2
    results_SVM_2[[paste0("Experiment",experiment_day)]] <- res2
    
    wrongly_predicted_3 <- which(res3$prediction != z_scored_test_data$diagnosis)
    
    if (length(wrongly_predicted_3) != 0){
      wrongly_3 <- data.frame("Patient" = z_scored_test_data$`SAMPLE ID`[wrongly_predicted_3],
                              "Diagnosis" = z_scored_test_data$diagnosis[wrongly_predicted_3],
                              "Wrong_diagnosis" = res3$prediction[wrongly_predicted_3])
      
      cat(paste0("For experiment day ", experiment_day, " patients wrongly predicted for three classes: "), file = logFile, append = TRUE, sep = "\n")
      capture.output( print(wrongly_3, print.gap=3), append = TRUE, file=logFile)
    }
    
    wrongly_predicted_2 <- which(res2$prediction != test_data$diagnosis_2)
    
    if (length(wrongly_predicted_2) != 0){
      wrongly_2 <- data.frame("Patient" = z_scored_test_data$`SAMPLE ID`[wrongly_predicted_2],
                              "Diagnosis" = z_scored_test_data$diagnosis_2[wrongly_predicted_2],
                              "Wrong_diagnosis" = res2$prediction[wrongly_predicted_2])
      cat(paste0("For experiment day ", experiment_day, " patients wrongly predicted for two classes: "), file = logFile, append = TRUE, sep = "\n")
      capture.output( print(wrongly_2, print.gap=3), append = TRUE, file=logFile)
      cat("\n", file = logFile, append = TRUE)
      
    }
    
    
  }
  
  save(results_SVM_2, results_SVM_3, file = file.path(working_directory, paste0("SVM_", kernel, "_NoFS_", to_use_data, ".Rdata")))
  
  # Build confusion matrix
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
  
}
