library(tidyverse)
library(readxl)

# ------------------------- reading in tube files -----------------------------

allfiles <- list.files("../../PhD/CVID/Excelsheets_data_verwerking", pattern = ".*_all_data" ,recursive = TRUE, full.names = TRUE)

tube1 <- read_excel( allfiles[1], sheet = "Tube 1 edited") #PBMCs
tube1 <- tube1[!duplicated(tube1$`SAMPLE ID`),]
tube2 <- read_excel( allfiles[2], sheet = "Tube 2 edited") #B cells
tube2 <- tube2[!duplicated(tube2$`SAMPLE ID`),]
tube5 <- read_excel( allfiles[5], sheet = "Tube 5 edited") #T cells
tube5 <- tube5[!duplicated(tube5$`SAMPLE ID`),]


big_manual_gating_matrix <- cbind(tube1,tube2, tube5)
big_manual_gating_matrix <- big_manual_gating_matrix[,!duplicated(colnames(big_manual_gating_matrix))]

# ---------------------- Make big matrix --------------------------------------

big_manual_gating_matrix <- big_manual_gating_matrix %>% 
  dplyr::filter(big_manual_gating_matrix$`SAMPLE ID` %in% 
                  meta_data$SAMPLE_ID)
meta_data <- meta_data[big_manual_gating_matrix$`SAMPLE ID`,]

big_manual_gating_matrix <- cbind("Diagnosis" = meta_data$diagnosis_group,
                                  big_manual_gating_matrix)
big_manual_gating_matrix <- big_manual_gating_matrix[,-c(3,4)]

big_manual_gating_matrix <- big_manual_gating_matrix %>% select(Diagnosis, 
                                                                `SAMPLE ID`,  
                                                                `EXPERIMENT NAME`, 
                                                                everything())

saveRDS(big_manual_gating_matrix, "Manual_gating/Data/tube1til5.rds")



# ------------------- Make individual tube matrices ---------------------------



for (tube in c(1,2,5)){
  if (tube == 1){ tube_matrix <- tube1
  } else if (tube == 2) {tube_matrix <- tube2
  } else{tube_matrix <- tube5}
  
  big_manual_gating_matrix <- tube_matrix %>% 
    dplyr::filter(tube5$`SAMPLE ID` %in% meta_data$SAMPLE_ID)
  meta_data <- meta_data[big_manual_gating_matrix$`SAMPLE ID`,]
  
  big_manual_gating_matrix <- cbind("Diagnosis" = meta_data$diagnosis_group,
                                    big_manual_gating_matrix)
  big_manual_gating_matrix <- big_manual_gating_matrix[,-c(3,4)]
  big_manual_gating_matrix <- big_manual_gating_matrix %>% select(Diagnosis, 
                                                                  `SAMPLE ID`,  
                                                                  `EXPERIMENT NAME`, 
                                                                  everything())
  
  
  saveRDS(big_manual_gating_matrix, paste0("Manual_gating/Data/tube",tube,".rds"))
  
}

# ---------------------------- Class to use -----------------------------------

HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency")

group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)

seed <- 1
