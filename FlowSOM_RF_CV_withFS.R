library(tidyverse)
library(randomForest)
library(crayon)
library(qsub)

tube <- "1_2_5" #(1 = PBMCs, 2 = B cells, 5 = T cells, 1_2_5 = Combination of all tubes)
seed <- 1

# --------------- With what group will we be working on? ----------------------
HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency")

group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID_New_preprocessing" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)


# ------ List directories and the fsom objects with the feature matrices ------
full_directory <- paste0("FlowSOM/Tube",tube, "/Data",PID,"/Crossvalidation")
working_directory <- file.path(full_directory, "RF_WithFS")
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

meta_clusters <- cbind("percentages_metaclusters", "MFI_metaclusters")
clusters <- cbind("percentages_clusters", "MFI_clusters")
total <- cbind("percentages_clusters", "percentages_metaclusters", 
               "percentages_to_clusters_metaclusters", "MFI_clusters", 
               "MFI_metaclusters")

all <- c("meta_clusters", "clusters", "percentages_clusters", "percentages_metaclusters","total", "percentages_to")


to_use_data <- all[1]


# Define what should be done in the file. Should it run on the cluster, extract the results from the cluster or built the models (calculate)?
cluster <- "extract"

# Define for what experiments the code should run
experiment_list <- c(13:33)

# ------------------------- z-score + RF---------------------------------------
#Load in data object
experiment_day <- 13
source("R_Scripts/Z_score.R")

for (to_use_data in all){
  
  
  
  #set up log script
  logFile = paste0(Logfile_directory,"/log_file_for_",
                   to_use_data,
                   "_FlowSOM_WithFS_RF.txt")
  cat(paste0("This is a log file for the patients that were predicted wrong from the RF without a certain experiment day where following features were used: ", to_use_data), file=logFile, append=FALSE, sep = "\n")
  
  
  results_RF_3 <- list()
  results_RF_2 <- list()
  selected_featureslist_2 <- list()
  selected_featureslist_3 <- list()
  
  
  for (experiment_day in experiment_list){
    set.seed(seed)
    print(paste0("Now calculating experiment ",experiment_day))
    
    # Make the correct featureset
    
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
    
    combined_data[,-1] <- apply(combined_data[,-1], 2, as.numeric)
    
    suppressWarnings( meta_data <- meta_data %>% 
      dplyr::filter(`SAMPLE_ID` %in% 
                      combined_data$`SAMPLE ID`) %>% 
      arrange(`SAMPLE_ID`))
    
    
    combined_data <- combined_data %>% 
      .[!duplicated(.[,"SAMPLE ID"]),] %>%  
      mutate("diagnosis" = meta_data$diagnosis_group) 
    
    combined_data <- combined_data %>% mutate("diagnosis_2" = ifelse(combined_data$diagnosis == "CVID", "CVID", "No_CVID"))
    
    
    to_compare_tube <- tube
    
    # This line is added in order to have the original filenames to compare to.
    if (tube == "1_2_5"){to_compare_tube <- 1}
    
    files <- list.files(paste0("Preprocessing/Tube",to_compare_tube,"/"), 
                        recursive = TRUE, 
                        pattern = ".*_pp.fcs")
    
 
    # Determine test and training data -- only perform wlicoxon or Kruskal wallis test on train data   
    test_data_names <- files[grep(paste0("Exp_",experiment_day),files)] %>% 
      gsub('.*/(.*)_.*_.*_.*','\\1',.)
    
    test_data_names <- test_data_names[which(test_data_names %in% meta_data$SAMPLE_ID)]
    train_data_names <- combined_data$`SAMPLE ID`[-which(test_data_names %in% combined_data$`SAMPLE ID`)]  
    train_data_names <- train_data_names[which(train_data_names %in% meta_data$SAMPLE_ID)]
    
    test_data <- combined_data %>% dplyr::filter(`SAMPLE ID` %in% test_data_names)
    train_data <- combined_data %>% dplyr::filter(`SAMPLE ID` %in% train_data_names)
    
  # Now the there will be a seperation for the predicion of two (CVID and no_CVID) and three classes (CVID, HCs, other PAD)
    
    dimension <- train_data %>% colnames %>% length
    
    train <- suppressWarnings(UseZScore(to_be_scored = train_data, 
                                        cols_not_to_convert = c(1,dimension -1, dimension), 
                                        meta_data = meta_data, 
                                        verbose = FALSE,
                                        start = 2,
                                        stop = dimension-2))
    
    z_scored_train_data <- train$z_score_matrix
    
    
    z_scored_test_data <- suppressWarnings(UseZScoreTestData(test_data, 
                                                             cols_not_to_convert = c(1, dimension - 1, dimension), 
                                                             meta_data, 
                                                             HC_scores = train$HC_scores,
                                                             start = 2, 
                                                             stop = dimension-2, 
                                                             verbose = FALSE))
    
    if (cluster == "run"){
      
      # Perform Wilcoxon test for every extracted features for the two-class prediction
      
      if (file.exists(file.path(working_directory, paste0("Randomforests_WithFS_handle2_exp_",experiment_day,"_", to_use_data, ".RDS"))) == FALSE){
    handle_2 <- qsub_lapply(
      X = c(4:dim(z_scored_train_data)[2]),
      qsub_config = override_qsub_config(
        name = "Wilcox_test",
        verbose = F,
        wait = F,
        memory = "4G"
      ),
      qsub_environment = c("z_scored_train_data"),
      FUN = function(x){
        set.seed(1)
        library(tidyverse)
        test <- wilcox.test(z_scored_train_data %>% filter(diagnosis_2 == "CVID") %>%
                              .[[x]],
                            z_scored_train_data %>% filter(diagnosis_2 == "No_CVID") %>%
                              .[[x]])
        p_value <- test$p.value
        names(p_value) <- colnames(z_scored_train_data)[x]
        return(p_value)
      } )
    saveRDS(handle_2, file = file.path(working_directory, paste0("Randomforests_WithFS_handle2_exp_",experiment_day,"_", to_use_data, ".RDS")))
      } else { print("You already have this file ;)")}
      
    }
    
    
    if (cluster == "extract"){
      # Extract and remove stored file
      
      if (file.exists(file.path(working_directory, paste0("Randomforests_WithFS_handle2_exp_",experiment_day,"_", to_use_data, ".RDS")))){
        
      qsub_obj_2 <- readRDS(file = file.path(working_directory, paste0("Randomforests_WithFS_handle2_exp_",experiment_day,"_", to_use_data, ".RDS")))
      
      if (file.exists(qsub_obj_2$src_dir) == TRUE){ 
      
      handle_2 <- qsub_retrieve(qsub_obj_2)
      
      } else {
        print (paste0( "Did you really run this handle 2 file for experiment day ", experiment_day, " for ", to_use_data, "?"))
        next()}
      
      saveRDS(handle_2, file = file.path(working_directory, paste0("Randomforests_WithFS_handle2_retrieved_exp_", experiment_day, "_", to_use_data, ".RDS")))
      unlink(file.path(working_directory, paste0("Randomforests_WithFS_handle2_exp_",experiment_day,"_", to_use_data, ".RDS")))
      }  else { print (paste0("file from experiment day ", experiment_day, " for " , to_use_data, " is not here for the handle 2"))} 
      
    }
    
    if (cluster == "calculate"){
      # Remove correlated features by feature selection
      handle_2 <- readRDS(file = file.path(working_directory, paste0("Randomforests_WithFS_handle2_retrieved_exp_", experiment_day, "_", to_use_data, ".RDS")))
      na_removed_2 <- unlist(handle_2[!is.na(handle_2)])
      na_removed_2 <- na.omit(na_removed_2)
      adj_p_2 <- stats::p.adjust(na_removed_2, "BH")
      sorted_2 <- sort(adj_p_2)
      z_scored_train_data <- z_scored_train_data %>% select("SAMPLE ID", names(sorted_2), "diagnosis", "diagnosis_2")
      selected_features_2 <- c(2,3)
      for(i_2 in (4:length(sorted_2))-2){
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
      selected_featureslist_2[[paste0("Experiment",experiment_day)]] <- selected_featurenames_2
  
      selected_train_data_2 <- z_scored_train_data %>% select("SAMPLE ID", "diagnosis", "diagnosis_2", selected_featurenames_2)
      selected_test_data_2 <- z_scored_test_data %>% select("SAMPLE ID", "diagnosis", "diagnosis_2", selected_featurenames_2)
  
  
      res2 <- randomForest(x = selected_train_data_2[,-c(1:3)],
                           y = as.factor(selected_train_data_2$diagnosis_2))
  
      res2$prediction <- predict(res2, selected_test_data_2[,-c(1:3)])
      res2$testset <- selected_test_data_2$`SAMPLE ID`
      res2$testsetpheno <- selected_test_data_2$diagnosis_2
      results_RF_2[[paste0("Experiment",experiment_day)]] <- res2

    }
    
    
    ## Perfrom Kruskal wallis tests for the three class prediction
  if (cluster == "run"){
    if (file.exists(file.path(working_directory, paste0("Randomforests_WithFS_handle3_exp_",experiment_day,"_", to_use_data, ".RDS"))) == FALSE){
      
    handle_3 <- qsub_lapply(
      X = c(4:(dim(z_scored_train_data)[2])),
      qsub_config = override_qsub_config(
        name = "Kruskal_Wallis_test",
        verbose = F,
        wait = F,
        memory = "4G"
      ),
      qsub_environment = c("z_scored_train_data"),
      FUN = function(x){
        set.seed(1)
        library(tidyverse)
        train_data_new <-  z_scored_train_data %>% select(x,'diagnosis')
        colnames(train_data_new) <- sub("-", "_", colnames(train_data_new))
        colnames(train_data_new) <- sub("/", "_", colnames(train_data_new))

        name <- colnames(train_data_new)[1]
        test <- kruskal.test(formula(paste(name ,"~diagnosis")), data = train_data_new )


        p_value <- test$p.value
        names(p_value) <- colnames(z_scored_train_data)[x]
        return(p_value)
      } )

    saveRDS(handle_3, file = file.path(working_directory, paste0("Randomforests_WithFS_handle3_exp_",experiment_day,"_", to_use_data, ".RDS")))
    } else { print("You already have this file ;)")}
  }
    if (cluster == "extract"){
      if (file.exists(file.path(working_directory, paste0("Randomforests_WithFS_handle3_exp_",experiment_day,"_", to_use_data, ".RDS")))){
      
      qsub_obj_3 <- readRDS(file = file.path(working_directory, paste0("Randomforests_WithFS_handle3_exp_",experiment_day,"_", to_use_data, ".RDS")))
      if (file.exists(qsub_obj_3$src_dir) == TRUE){ 
        
        handle_3 <- qsub_retrieve(qsub_obj_3)
        
      } else {
        print (paste0( "Did you really run this handle 3 file for experiment day ", experiment_day, " for ", to_use_data, "?"))
        next()}
      
      saveRDS(handle_3, file = file.path(working_directory, paste0("Randomforests_WithFS_handle3_retrieved_exp_", experiment_day, "_", to_use_data, ".RDS")))
      unlink(file.path(working_directory, paste0("Randomforests_WithFS_handle3_exp_",experiment_day,"_", to_use_data, ".RDS"))) 
      } else { print (paste0("file from experiment day ", experiment_day, " for " , to_use_data, " is not here for the handle 3"))} 
      
    }
    if (cluster == "calculate") {
    
      # Remove correlated features by feature selection
      
    handle_3 <- readRDS(file = file.path(working_directory, paste0("Randomforests_WithFS_handle3_retrieved_exp_", experiment_day, "_", to_use_data, ".RDS")))
    na_removed_3 <- unlist(handle_3[!is.na(handle_3)])
    adj_p_3 <- stats::p.adjust(na_removed_3, "BH")
    sorted_3 <- sort(adj_p_3)
    z_scored_train_data <- z_scored_train_data %>% select("SAMPLE ID", names(sorted_3), "diagnosis", "diagnosis_2")
    selected_features_3 <- c(2,3)
    for(i_3 in (4:length(sorted_3))-2){
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
    selected_featureslist_3[[paste0("Experiment",experiment_day)]] <- selected_featurenames_3

    selected_train_data_3 <- z_scored_train_data %>% select("SAMPLE ID", "diagnosis", "diagnosis_2", selected_featurenames_3)
    selected_test_data_3 <- z_scored_test_data %>% select("SAMPLE ID", "diagnosis", "diagnosis_2",selected_featurenames_3)



    res3 <- randomForest(x = selected_train_data_3[,-c(1:3)],
                         y = selected_train_data_3$diagnosis)

    res3$prediction <- predict(res3, selected_test_data_3[,-c(1:3)])
    res3$testset <- selected_test_data_3$`SAMPLE ID`
    res3$testsetpheno <- selected_test_data_3$diagnosis
    results_RF_3[[paste0("Experiment",experiment_day)]] <- res3


    wrongly_predicted_3 <- which(res3$prediction != selected_test_data_3$diagnosis)

    if (length(wrongly_predicted_3) != 0){
      wrongly_3 <- data.frame("Patient" = selected_test_data_3$`SAMPLE ID`[wrongly_predicted_3],
                              "Diagnosis" = selected_test_data_3$diagnosis[wrongly_predicted_3],
                              "Wrong_diagnosis" = res3$prediction[wrongly_predicted_3])

      cat(paste0("For experiment day ", experiment_day, " patients wrongly predicted for three classes: "), file = logFile, append = TRUE, sep = "\n")
      capture.output( print(wrongly_3, print.gap=3), append = TRUE, file=logFile)
    }

    wrongly_predicted_2 <- which(res2$prediction != selected_test_data_2$diagnosis_2)

    if (length(wrongly_predicted_2) != 0){
      wrongly_2 <- data.frame("Patient" = selected_test_data_2$`SAMPLE ID`[wrongly_predicted_2],
                              "Diagnosis" = selected_test_data_2$diagnosis[wrongly_predicted_2],
                              "Wrong_diagnosis" = res2$prediction[wrongly_predicted_2])
      cat(paste0("For experiment day ", experiment_day, " patients wrongly predicted for two classes: "), file = logFile, append = TRUE, sep = "\n")
      capture.output( print(wrongly_2, print.gap=3), append = TRUE, file=logFile)
      cat("\n", file = logFile, append = TRUE)

    }
    }
    

  }
  
  if (cluster == "calculate"){
    save(results_RF_2, results_RF_3, selected_featureslist_2, selected_featureslist_3, file = file.path(working_directory, paste0("Randomforests_WithFS_", to_use_data, ".Rdata")))
  
  
    # Make confusion matrix
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
  }
  
}

