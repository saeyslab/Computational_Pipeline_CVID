library(tidyverse)

tube <- "1_2_5"
seed <- 1

# ------- Load in name matrices ----------------------------------------------

#names <- readRDS("FlowSOM/Tube1/Data/Crossvalidation/Cluster_Metacluster_Labels.Rds")


# --------------- With what group will we be working on? ----------------------
HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency", "HFM")

group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID_New_preprocessing" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)

# ------ List directories and the fsom objects with the feature matrices ------
full_directory <- paste0("FlowSOM/Tube",tube, "/Data", PID,"/Crossvalidation")

working_directory <- file.path(full_directory, "fsom_objects")

storing_directory <- file.path(full_directory, "Matrices")


suppressWarnings(dir.create(storing_directory))

files_fsoms <- list.files(working_directory, pattern = "fsom_objects_experiment_day...Rds")

fcs_dir = paste0("Preprocessing/Tube",tube)

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

allfiles <- list.files(fcs_dir, pattern = "^P.*_pp.fcs$", recursive = TRUE)
files <- allfiles[grep(paste(rownames(meta_data), collapse="|"), allfiles)]

# --------------------------- Make orderly matrices ---------------------------


experiment_day <- 13

for (experiment_day in c(14:33)){
  
  # For entire tube
   # full_directory <- paste0("FlowSOM/Tube",tube, "/Data", PID)
   # file1 <- "fsom_objects.Rds"
   # results <- readRDS(file.path(full_directory, file1))
   # 
  
  # For different experiments
  file1 <- files_fsoms[grep(experiment_day, files_fsoms)]
  results <- readRDS(file.path(working_directory, file1))
  
  patient_names <- sapply(results, function(x)x[[1]])

  
  percentages_clusters <- matrix(0, nrow = length(patient_names), 
                                 ncol = length(results[[1]]$percentages_clusters), 
                                 dimnames = list(patient_names,
                                                 paste("c",seq_len(length(results[[1]]$percentages_clusters)), sep = "_")))
  

  percentages_metaclusters <- matrix(0, nrow = length(patient_names), 
                                     ncol = length(results[[1]]$percentages_metaclusters), 
                                     dimnames = list(patient_names,
                                                     paste("mc",seq_len(length(results[[1]]$percentages_metaclusters)), sep = "_")))
   
  percentages_to_clusters_metaclusters <- matrix(0, nrow = length(patient_names), 
                                                 ncol = length(results[[1]]$percentages_cluster2metacluster), 
                                                 dimnames = list(patient_names,
                                                                 paste("ctomc",seq_len(length(results[[1]]$percentages_cluster2metacluster)),  sep = "_")))
  
  names_MFI_clusters <- lapply(colnames(results[[1]]$MFI_clusters),
                               function(x) paste("c",x, 1:dim(percentages_clusters)[2], sep="_")) %>% 
    unlist %>% 
    unname
  
  MFI_clusters <- matrix(0, nrow = length(patient_names), 
                                 ncol = length(names_MFI_clusters), 
                                 dimnames = list(patient_names,
                                                 names_MFI_clusters))
  
   
  names_MFI_metaclusters <- lapply(colnames(results[[1]]$MFI_metaclusters),
                                   function(x) paste("mc",x, 1:dim(percentages_metaclusters)[2], sep="_")) %>% 
    unlist %>% 
    unname
  
  MFI_metaclusters <- matrix(0, nrow = length(patient_names), 
                         ncol = length(names_MFI_metaclusters), 
                         dimnames = list(patient_names,
                                         names_MFI_metaclusters))
  
  patient <- results[[32]]
  for (patient in results){
    percentages_clusters[patient$Patient,] <- patient$percentages_clusters
    percentages_metaclusters[patient$Patient,] <- patient$percentages_metaclusters
    percentages_to_clusters_metaclusters[patient$Patient,] <- patient$percentages_cluster2metacluster
    
    
    print(patient$Patient)
    for (colour in colnames(results[[1]]$MFI_clusters)){
      MFI_clusters[patient$Patient, grep(colour, colnames(MFI_clusters))] <- patient$MFI_clusters[,colour] %>% ifelse(is.na(.), 0, .)
      MFI_metaclusters[patient$Patient,grep(colour, colnames(MFI_metaclusters))] <- patient$MFI_metaclusters[,colour] %>% ifelse(is.na(.), 0, .)
    }
    
  }

  
  # For entire tube
  # saveRDS(percentages_clusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/percentages_clusters.RDS"))
  # saveRDS(percentages_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/percentages_metaclusters.RDS"))
  # saveRDS(percentages_to_clusters_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/percentages_to_clusters_metaclusters.RDS"))
  # saveRDS(MFI_clusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/MFI_clusters.RDS"))
  # saveRDS(MFI_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/MFI_metaclusters.RDS"))
  # 
  # For individual tubes
  saveRDS(percentages_clusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/percentages_clusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(percentages_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/percentages_metaclusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(percentages_to_clusters_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/percentages_to_clusters_metaclusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(MFI_clusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/MFI_clusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(MFI_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/MFI_metaclusters_experiment_day_", experiment_day, ".RDS"))

}



# --------------------- Make percentages for tube 1 and 2 and 5 ---------------

tube <- 1
full_directory_1 <- paste0("FlowSOM/Tube",tube, "/Data", PID,"/Crossvalidation")
working_directory_1 <- file.path(full_directory_1, "Matrices")

tube <- 2
full_directory_2 <- paste0("FlowSOM/Tube",tube, "/Data", PID,"/Crossvalidation")
working_directory_2 <- file.path(full_directory_2, "Matrices")

tube <- 5
full_directory_5 <- paste0("FlowSOM/Tube",tube, "/Data", PID,"/Crossvalidation")
working_directory_5 <- file.path(full_directory_5, "Matrices")


tube <- "1_2"
full_directory_3 <- paste0("FlowSOM/Tube",tube, "/Data", PID,"/Crossvalidation")
storing_directory <- file.path(full_directory_3, "Matrices")

tube <- "1_2_5"
full_directory_full <- paste0("FlowSOM/Tube",tube, "/Data", PID,"/Crossvalidation")
storing_directory <- file.path(full_directory_full, "Matrices")



experiment_day <- 13
for (experiment_day in c(13:33)){
  print(experiment_day)
  percentages_metaclusters_1 <- readRDS(file.path(working_directory_1, paste0("percentages_metaclusters_experiment_day_",experiment_day, ".RDS")))
  MFI_metaclusters_1 <- readRDS(file.path(working_directory_1, paste0("MFI_metaclusters_experiment_day_",experiment_day, ".RDS")))
  percentages_clusters_1 <- readRDS(file.path(working_directory_1, paste0("percentages_clusters_experiment_day_",experiment_day, ".RDS")))
  MFI_clusters_1 <- readRDS(file.path(working_directory_1, paste0("MFI_clusters_experiment_day_",experiment_day, ".RDS")))
  percentages_to_clusters_metaclusters_1 <- readRDS(file.path(working_directory_1, paste0("percentages_to_clusters_metaclusters_experiment_day_",experiment_day, ".RDS")))
  percentages_metaclusters_2 <- readRDS(file.path(working_directory_2, paste0("percentages_metaclusters_experiment_day_",experiment_day, ".RDS")))
  MFI_metaclusters_2 <- readRDS(file.path(working_directory_2, paste0("MFI_metaclusters_experiment_day_",experiment_day, ".RDS")))
  percentages_clusters_2 <- readRDS(file.path(working_directory_2, paste0("percentages_clusters_experiment_day_",experiment_day, ".RDS")))
  MFI_clusters_2 <- readRDS(file.path(working_directory_2, paste0("MFI_clusters_experiment_day_",experiment_day, ".RDS")))
  percentages_to_clusters_metaclusters_2 <- readRDS(file.path(working_directory_2, paste0("percentages_to_clusters_metaclusters_experiment_day_",experiment_day, ".RDS")))
  percentages_metaclusters_5 <- readRDS(file.path(working_directory_5, paste0("percentages_metaclusters_experiment_day_",experiment_day, ".RDS")))
  MFI_metaclusters_5 <- readRDS(file.path(working_directory_5, paste0("MFI_metaclusters_experiment_day_",experiment_day, ".RDS")))
  percentages_clusters_5 <- readRDS(file.path(working_directory_5, paste0("percentages_clusters_experiment_day_",experiment_day, ".RDS")))
  MFI_clusters_5 <- readRDS(file.path(working_directory_5, paste0("MFI_clusters_experiment_day_",experiment_day, ".RDS")))
  percentages_to_clusters_metaclusters_5 <- readRDS(file.path(working_directory_5, paste0("percentages_to_clusters_metaclusters_experiment_day_",experiment_day, ".RDS")))
  
  
  MFI_metaclusters_1 <- MFI_metaclusters_1[rownames(MFI_metaclusters_2),]
  MFI_clusters_1 <- MFI_clusters_1[rownames(MFI_clusters_2),]
  percentages_clusters_1 <- percentages_clusters_1[rownames(percentages_clusters_2),]
  percentages_metaclusters_1 <- percentages_metaclusters_1[rownames(percentages_metaclusters_2),]
  percentages_to_clusters_metaclusters_1 <- percentages_to_clusters_metaclusters_1[rownames(percentages_to_clusters_metaclusters_2),]
  
  
  colnames(percentages_metaclusters_1) <- paste0("Tube1_", colnames(percentages_metaclusters_1))
  colnames(percentages_clusters_1) <- paste0("Tube1_", colnames(percentages_clusters_1))
  colnames(percentages_to_clusters_metaclusters_1) <- paste0("Tube1_", colnames(percentages_to_clusters_metaclusters_1))
  colnames(MFI_metaclusters_1) <- paste0("Tube1_", colnames(MFI_metaclusters_1))
  colnames(MFI_clusters_1) <- paste0("Tube1_", colnames(MFI_clusters_1))
  
  
  colnames(percentages_metaclusters_2) <- paste0("Tube2_", colnames(percentages_metaclusters_2))
  colnames(percentages_clusters_2) <- paste0("Tube2_", colnames(percentages_clusters_2))
  colnames(percentages_to_clusters_metaclusters_2) <- paste0("Tube2_", colnames(percentages_to_clusters_metaclusters_2))
  colnames(MFI_metaclusters_2) <- paste0("Tube2_", colnames(MFI_metaclusters_2))
  colnames(MFI_clusters_2) <- paste0("Tube2_", colnames(MFI_clusters_2))
  
  
  colnames(percentages_metaclusters_5) <- paste0("Tube5_", colnames(percentages_metaclusters_5))
  colnames(percentages_clusters_5) <- paste0("Tube5_", colnames(percentages_clusters_5))
  colnames(percentages_to_clusters_metaclusters_5) <- paste0("Tube5_", colnames(percentages_to_clusters_metaclusters_5))
  colnames(MFI_metaclusters_5) <- paste0("Tube5_", colnames(MFI_metaclusters_5))
  colnames(MFI_clusters_5) <- paste0("Tube5_", colnames(MFI_clusters_5))
  
  
  percentages_clusters1_2 <- cbind(percentages_clusters_1, percentages_clusters_2)
  percentages_metaclusters1_2 <- cbind(percentages_metaclusters_1, percentages_metaclusters_2)
  percentages_to_clusters_metaclusters1_2 <-  cbind(percentages_to_clusters_metaclusters_1, percentages_to_clusters_metaclusters_2)
  MFI_clusters1_2 <- cbind(MFI_clusters_1, MFI_clusters_2)
  MFI_metaclusters1_2 <- cbind(MFI_metaclusters_1, MFI_metaclusters_2)
  
  percentages_clusters1_5 <- cbind(percentages_clusters_1, percentages_clusters_2, percentages_clusters_5)
  percentages_metaclusters1_5 <- cbind(percentages_metaclusters_1, percentages_metaclusters_2, percentages_metaclusters_5)
  percentages_to_clusters_metaclusters1_5 <-  cbind(percentages_to_clusters_metaclusters_1, percentages_to_clusters_metaclusters_2, percentages_to_clusters_metaclusters_5)
  MFI_clusters1_5 <- cbind(MFI_clusters_1, MFI_clusters_2, MFI_clusters_5)
  MFI_metaclusters1_5 <- cbind(MFI_metaclusters_1, MFI_metaclusters_2, MFI_metaclusters_5)
  
  tube <- "1_2_5"
  saveRDS(percentages_clusters1_5, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/percentages_clusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(percentages_metaclusters1_5, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/percentages_metaclusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(percentages_to_clusters_metaclusters1_5, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/percentages_to_clusters_metaclusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(MFI_clusters1_5, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/MFI_clusters_experiment_day_", experiment_day, ".RDS"))
  saveRDS(MFI_metaclusters1_5, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/Matrices/MFI_metaclusters_experiment_day_", experiment_day, ".RDS"))
  
  
}



# --------------------------- Make percentages for tube 1 and 2 for entire tube ---------------

tube <- 1
working_directory_1 <- paste0("FlowSOM/Tube",tube, "/Data", PID)

tube <- 2
working_directory_2 <- paste0("FlowSOM/Tube",tube, "/Data", PID)

tube <- 5
working_directory_5 <- paste0("FlowSOM/Tube",tube, "/Data", PID)


tube <- "1_2_5"
working_directory_3 <- paste0("FlowSOM/Tube",tube, "/Data", PID)


  percentages_metaclusters_1 <- readRDS(file.path(working_directory_1, paste0("percentages_metaclusters.RDS")))
  MFI_metaclusters_1 <- readRDS(file.path(working_directory_1, paste0("MFI_metaclusters.RDS")))
  percentages_clusters_1 <- readRDS(file.path(working_directory_1, paste0("percentages_clusters.RDS")))
  MFI_clusters_1 <- readRDS(file.path(working_directory_1, paste0("MFI_clusters.RDS")))
  percentages_to_clusters_metaclusters_1 <- readRDS(file.path(working_directory_1, paste0("percentages_to_clusters_metaclusters.RDS")))
  percentages_metaclusters_2 <- readRDS(file.path(working_directory_2, paste0("percentages_metaclusters.RDS")))
  MFI_metaclusters_2 <- readRDS(file.path(working_directory_2, paste0("MFI_metaclusters.RDS")))
  percentages_clusters_2 <- readRDS(file.path(working_directory_2, paste0("percentages_clusters.RDS")))
  MFI_clusters_2 <- readRDS(file.path(working_directory_2, paste0("MFI_clusters.RDS")))
  percentages_to_clusters_metaclusters_2 <- readRDS(file.path(working_directory_2, paste0("percentages_to_clusters_metaclusters.RDS")))
  percentages_metaclusters_5 <- readRDS(file.path(working_directory_5, paste0("percentages_metaclusters.RDS")))
  MFI_metaclusters_5 <- readRDS(file.path(working_directory_5, paste0("MFI_metaclusters.RDS")))
  percentages_clusters_5 <- readRDS(file.path(working_directory_5, paste0("percentages_clusters.RDS")))
  MFI_clusters_5 <- readRDS(file.path(working_directory_5, paste0("MFI_clusters.RDS")))
  percentages_to_clusters_metaclusters_5 <- readRDS(file.path(working_directory_5, paste0("percentages_to_clusters_metaclusters.RDS")))
  
  colnames(percentages_metaclusters_1) <- paste0("Tube1_", colnames(percentages_metaclusters_1))
  colnames(percentages_clusters_1) <- paste0("Tube1_", colnames(percentages_clusters_1))
  colnames(percentages_to_clusters_metaclusters_1) <- paste0("Tube1_", colnames(percentages_to_clusters_metaclusters_1))
  colnames(MFI_metaclusters_1) <- paste0("Tube1_", colnames(MFI_metaclusters_1))
  colnames(MFI_clusters_1) <- paste0("Tube1_", colnames(MFI_clusters_1))
  
  
  colnames(percentages_metaclusters_2) <- paste0("Tube2_", colnames(percentages_metaclusters_2))
  colnames(percentages_clusters_2) <- paste0("Tube2_", colnames(percentages_clusters_2))
  colnames(percentages_to_clusters_metaclusters_2) <- paste0("Tube2_", colnames(percentages_to_clusters_metaclusters_2))
  colnames(MFI_metaclusters_2) <- paste0("Tube2_", colnames(MFI_metaclusters_2))
  colnames(MFI_clusters_2) <- paste0("Tube2_", colnames(MFI_clusters_2))
  
  colnames(percentages_metaclusters_5) <- paste0("Tube5_", colnames(percentages_metaclusters_5))
  colnames(percentages_clusters_5) <- paste0("Tube5_", colnames(percentages_clusters_5))
  colnames(percentages_to_clusters_metaclusters_5) <- paste0("Tube5_", colnames(percentages_to_clusters_metaclusters_5))
  colnames(MFI_metaclusters_5) <- paste0("Tube5_", colnames(MFI_metaclusters_5))
  colnames(MFI_clusters_5) <- paste0("Tube5_", colnames(MFI_clusters_5))
  
  MFI_clusters_1 <- MFI_clusters_1[rownames(MFI_clusters_2),]
  MFI_metaclusters_1 <- MFI_metaclusters_1[rownames(MFI_metaclusters_2),]
  percentages_clusters_1 <-  percentages_clusters_1[rownames( percentages_clusters_2),]
  percentages_metaclusters_1 <-  percentages_metaclusters_1[rownames( percentages_metaclusters_2),]
  percentages_to_clusters_metaclusters_1 <- percentages_to_clusters_metaclusters_1[rownames(percentages_to_clusters_metaclusters_2),]
  
  
  percentages_clusters <- cbind(percentages_clusters_1, percentages_clusters_2, percentages_clusters_5)
  percentages_metaclusters <- cbind(percentages_metaclusters_1, percentages_metaclusters_2,percentages_metaclusters_5)
  percentages_to_clusters_metaclusters <-  cbind(percentages_to_clusters_metaclusters_1, percentages_to_clusters_metaclusters_2, percentages_to_clusters_metaclusters_5)
  MFI_clusters <- cbind(MFI_clusters_1, MFI_clusters_2, MFI_clusters_5)
  MFI_metaclusters <- cbind(MFI_metaclusters_1, MFI_metaclusters_2, MFI_metaclusters_5)
  
  tube <- "1_2_5"
  saveRDS(percentages_clusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/percentages_clusters.RDS"))
  saveRDS(percentages_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/percentages_metaclusters.RDS"))
  saveRDS(percentages_to_clusters_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/percentages_to_clusters_metaclusters.RDS"))
  saveRDS(MFI_clusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/MFI_clusters.RDS"))
  saveRDS(MFI_metaclusters, file = paste0("FlowSOM/Tube", tube, "/Data", PID, "/MFI_metaclusters.RDS"))




