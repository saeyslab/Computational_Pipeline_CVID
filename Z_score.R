CalculateZScore <- function(HC_matrix,cols_not_to_convert){
  mean <- apply(HC_matrix[,-cols_not_to_convert], 2, mean)
  sd <- apply(HC_matrix[,-cols_not_to_convert], 2, sd)
  return(list("mean" = mean, "sd" = sd))
}


## This function was created to try to eliminate the age-related immunechanges between patient groups
# The mean and standard deviation are calculated per age group for the healthy controls 
# All the patients and healthy controls belonging to this age group are then corrected by using the z-score based on these two measurements

UseZScore <- function(to_be_scored, cols_not_to_convert, meta_data, verbose = FALSE, start, stop){
  library(dplyr)
  ## No age HC so no z-score could be applied --> to remove from the z-scored to_be_scored
  end <- list()
  
  HC_scores <- list()
  
  to_be_scored[,-cols_not_to_convert] <- apply(to_be_scored[-cols_not_to_convert], 2, as.numeric)
  
  agegroup <- 3
  for (agegroup in c(2:10)){
    
    if (verbose == TRUE){
      print(paste0("Calculating following agegroup: ",agegroup))}
    
    Patients_age_names <- meta_data %>% dplyr::filter(age_group_Z_scores == agegroup) %>% .$`SAMPLE_ID`
    
    Patients_this_age <- to_be_scored %>% dplyr::filter(`SAMPLE ID` %in% Patients_age_names)
    
    HC_this_age_names <- meta_data %>% dplyr::filter(diagnosis_group == "HC" & age_group_Z_scores == agegroup) %>% .$`SAMPLE_ID`
    
    HC_this_age <- to_be_scored %>% dplyr::filter(`SAMPLE ID` %in% HC_this_age_names)
    
    if (dim(HC_this_age)[1] < 2){
      next}
    
    HC_calculations <- CalculateZScore(HC_this_age, cols_not_to_convert)
    HC_scores[[as.character(agegroup)]] <- HC_calculations
    
    z_score_matrix <- sapply(c(start:stop), function(x)((Patients_this_age %>% .[[x]])- HC_calculations$mean[x - (start -1)])/HC_calculations$sd[x- (start - 1)])
    
    colnames(z_score_matrix) <- colnames(Patients_this_age)[-cols_not_to_convert]
    
    z_score_matrix <-  tbl_df(z_score_matrix)
    
    if (length(cols_not_to_convert) == 1){
    final_z_score_matrix <- bind_cols("SAMPLE ID" = Patients_this_age[,cols_not_to_convert], z_score_matrix)
    }else{final_z_score_matrix <- bind_cols(Patients_this_age[,cols_not_to_convert], z_score_matrix)}
    
    
    
    end[[as.character(agegroup)]] <- final_z_score_matrix
    
    
  }
  z_scored <- bind_rows(end)
  z_scored <- z_scored %>% arrange(`SAMPLE ID`)
  
  
  z_scored <- z_scored %>% replace(., is.na(.), 0)
  z_scored[which(z_scored == Inf, arr.ind = TRUE)] <- 0
  z_scored[which(z_scored == -Inf, arr.ind = TRUE)] <- 0
  
  
  return(list("z_score_matrix" = z_scored, "HC_scores" = HC_scores))
  
}


UseZScoreTestData <- function(test_data, cols_not_to_convert, meta_data, HC_scores, start, stop, verbose = TRUE){
  
  end <- list()
  test_data[,-cols_not_to_convert] <- apply(test_data[-cols_not_to_convert], 2, as.numeric)
  
  ## No age HC so no z-score could be applied --> to remove from the z-scored giant_matrix
  start_agegroup <- as.numeric(names(HC_scores)[1])
  
  agegroup <- 6
  for (agegroup in c(start_agegroup:10)){
    
    if (verbose == TRUE){
      print(paste0("Calculating following agegroup: ",agegroup))}
    
    Patients_age_names <- meta_data %>% dplyr::filter(age_group_Z_scores == agegroup) %>% .$`SAMPLE_ID`
    
    Patients_this_age <- test_data %>% dplyr::filter(`SAMPLE ID` %in% Patients_age_names)
    
    if(dim(Patients_this_age)[1] == 0 ){
      next
    }
  
    
    HC_calculations <- HC_scores[[as.character(agegroup)]]
    
    z_score_matrix <- sapply(c(start:stop), function(x)((Patients_this_age %>% .[[x]])- HC_calculations$mean[x - (start -1)])/HC_calculations$sd[x- (start - 1)])
    
    if (dim(Patients_this_age)[1] >1){
    colnames(z_score_matrix) <- colnames(Patients_this_age)[-cols_not_to_convert]
    
    z_score_matrix <-  tbl_df(z_score_matrix)
    
    if (length(cols_not_to_convert) == 1){
      final_z_score_matrix <- bind_cols("SAMPLE ID" = Patients_this_age[,cols_not_to_convert], z_score_matrix)
    }else{final_z_score_matrix <- bind_cols(Patients_this_age[,cols_not_to_convert], z_score_matrix)}
    } else { 
      names(z_score_matrix) <- colnames(Patients_this_age)[-cols_not_to_convert]
      final_z_score_matrix <- c(Patients_this_age[,cols_not_to_convert], z_score_matrix)
    }
    
    end[[as.character(agegroup)]] <- final_z_score_matrix
    
    
  }
  z_scored <- bind_rows(end)
  z_scored <- z_scored %>% arrange(`SAMPLE ID`)
  
  z_scored <- z_scored %>% replace(., is.na(.), 0)
  z_scored[which(z_scored == Inf, arr.ind = TRUE)] <- 0
  return(z_scored)
  
}
