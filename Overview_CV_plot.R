library(ggplot2)
library(dplyr)
library(randomForest)


HC_CVID <- c("HC", "CVID")
PID_HC <- c(HC_CVID, "IgG_deficiency", "IgG_subclass_deficiency", "IgM_deficiency",
            "IgA_deficiency", "Syndromic_hypogamma")


group_to_use <- PID_HC  # Choose between HC_CVID and PID_HC
PID <- "_PID_New_preprocessing" #Set to empty character if you are only working with HC en CVID (otherwise: _PID)


# --------------------------- Figure Dir --------------------------------------

figure_dir <- paste0("FlowSOM/Tube",tube,"/Figures",PID,"/Crossvalidation/")
suppressWarnings(dir.create(figure_dir))
specific_figure_dir <- file.path(figure_dir, FS)
suppressWarnings(dir.create(specific_figure_dir))
randomforest_result_dir <- paste0("FlowSOM/Tube", tube, "/Data", PID, "/Crossvalidation/RF_", FS ) 


# ------------------------- Make nice overview of results ---------------------



type <- "Manual_gating"

working_dir <- file.path("Manual_gating/Data")

tube <- "Full"

name_classifier <- "Randomforests"

final_results <- data.frame("Tube" =NA, "Type" =NA, "FS"=NA, "Nr_Classes"=NA, "Mistakes"=NA, "Specificity"=NA, "Sensitivity"=NA, "Classifier" = NA, "Precision" = NA, "Recall" = NA, "F_measure" = NA, "Balanced_accuracy" = NA) 
#colnames(final_results) <- c("Tube", "Type", "FS", "Nr_Classes", "Mistakes", "Specificity", "Sensitivity", "Classifier")


for(name_classifier in c("Randomforests", "SVM")){
for (FS in c("WithFS", "NoFS")){
  for (tube in c(1,2,5, "Full")){
   
    load(file.path(working_dir,paste0(name_classifier,"_",FS,"_Manual_Gating_tube",tube,".Rdata")))
    mistakes_2Classes <- c()
    mistakes_3Classes <- c()
    
    
    if (name_classifier == "Randomforests"){
      classifier <- "RF"
      experiment2 <- results_RF_2
      experiment3 <- results_RF_3
    } else { classifier <-  "SVM_linear"
    experiment2 <- results_SVM_2
    experiment3 <- results_SVM_3
    }
    final_2 <- matrix(0, nrow = 2, ncol = 2)
    
    for (experiment in experiment2){
      t_2 <- table(experiment$prediction, experiment$testsetpheno)
      if (("CVID" %in% colnames(t_2) == FALSE)){
        t_2 <- cbind("CVID"=c(0,0), t_2)
      }
      if (("No_CVID" %in% colnames(t_2) == FALSE)){
        t_2 <- cbind("CVID"=c(0,0), t_2)
      }
      wrongly_predicted_2 <- which(experiment$prediction != experiment$testsetpheno)
      if (length(wrongly_predicted_2) >0) {mistakes_2Classes <- c(mistakes_2Classes, experiment$testset[wrongly_predicted_2])}
      
      final_2 <- final_2 + t_2
    }
    final_3 <- matrix(0, nrow = 3, ncol = 3)
    for (experiment in experiment3){
      t_3 <- table(experiment$prediction, experiment$testsetpheno)
      final_3 <- final_3 + t_3
      wrongly_predicted_3 <- which(experiment$prediction != experiment$testsetpheno)
      if (length(wrongly_predicted_3) >0) {mistakes_3Classes <- c(mistakes_3Classes, experiment$testset[wrongly_predicted_3])}
      
    
      }
    
    sensitivity_2 <- final_2[1] / (final_2[1] + final_2[2])
    specificity_2 <- final_2[4] / (final_2[3] + final_2[4])
    
    
    sensitivity_3 <- final_3[1] / (final_3[1] + final_3[2] + final_3[3])
    specificity_3 <- (final_3[5] + final_3[6] + final_3[8] + final_3[9]) / (final_3[5] + final_3[6] + final_3[8] + final_3[9] + final_3[4] + final_3[7])
    
    precision_2 <- final_2[1] / (final_2[1] + final_2[3])
    recall_2 <- final_2[1] / (final_2[1] + final_2[2])
    
    f_2 <- 2*((precision_2*recall_2) / (precision_2 +recall_2))
    
    precision_3 <- final_3[1] / (final_3[1] + final_3[4] + final_3[7])
    recall_3 <- (final_3[1]) / (final_3[1] + final_3[2] + final_3[3])
    
    f_3 <- 2*((precision_3*recall_3) / (precision_3 +recall_3))
    
    
    balanced_2 <- ((final_2[1]/ (final_2[1] + final_2[2])) + (final_2[4]/ (final_2[3] + final_2[4])))/2
    balanced_3 <- ((final_3[1]/ (final_3[1] + final_3[2] + final_3[3])) + (final_3[5]/ (final_3[4] + final_3[5] + final_3[6])) + (final_3[9]/ (final_3[7] + final_3[8] + final_3[9])))/3
    
    if (tube == "Full"){
      tube <- "1_2_5"
    }
    
    results_2Classes <- data.frame("Tube" = tube, 
                                   "Type" = "Manual", 
                                   "FS"= FS, 
                                   "Nr_Classes"=2, 
                                   "Specificity"= specificity_2, 
                                   "Sensitivity"= sensitivity_2,
                                   "Precision" = precision_2,
                                   "Recall" = recall_2,
                                   "F_measure" = f_2, 
                                   "Classifier" = classifier,
                                   "Balanced_accuracy" = balanced_2 )
    results_2Classes$Mistakes <- list(mistakes_2Classes)
    results_3Classes <- data.frame("Tube" = tube, 
                                   "Type" = "Manual", 
                                   "FS"= FS, 
                                   "Nr_Classes"=3, 
                                   "Specificity"= specificity_3, 
                                   "Sensitivity"= sensitivity_3,
                                   "Precision" = precision_3,
                                   "Recall" = recall_3,
                                   "F_measure" = f_3, 
                                   "Classifier" = classifier,
                                   "Balanced_accuracy" = balanced_3)  
    results_3Classes$Mistakes <- list(mistakes_3Classes)
    
    final_results <- rbind(final_results, results_3Classes, results_2Classes)
  }
}
}
working_dir <- file.path("FlowSOM")

name_classifier <- "RF"
for(name_classifier in c("RF", "SVM_linear")){
  FS <- "NoFS"
  for (FS in c("WithFS", "NoFS")){
    tube <- 1
    for (tube in c(1,2, 5,"1_2_5")){
      
      working_dir_tube <- file.path(working_dir, paste0("Tube",tube), "Data_PID_New_preprocessing", "Crossvalidation", paste0(name_classifier, "_", FS))
      for (to_use_data in c("clusters", "meta_clusters", "percentages_clusters", "percentages_metaclusters", "percentages_to", "total")){
      if (name_classifier == "RF"){
        classifier <- "Randomforests"
      } else {classifier <- "SVM_linear"}
        
      load(file.path(working_dir_tube,paste0(classifier,"_",FS,"_", to_use_data, ".Rdata")))
      mistakes_2Classes <- c()
      mistakes_3Classes <- c()
      
      
      if (name_classifier == "RF"){
        experiment2 <- results_RF_2
        experiment3 <- results_RF_3
      } else {
      experiment2 <- results_SVM_2
      experiment3 <- results_SVM_3
      }
      final_2 <- matrix(0, nrow = 2, ncol = 2)
      
      for (experiment in experiment2){
        t_2 <- table(experiment$prediction, experiment$testsetpheno)
        if (("CVID" %in% colnames(t_2) == FALSE)){
          t_2 <- cbind("CVID"=c(0,0), t_2)
        }
        if (("No_CVID" %in% colnames(t_2) == FALSE)){
          t_2 <- cbind("CVID"=c(0,0), t_2)
        }
        wrongly_predicted_2 <- which(experiment$prediction != experiment$testsetpheno)
        if (length(wrongly_predicted_2) >0) {mistakes_2Classes <- c(mistakes_2Classes, experiment$testset[wrongly_predicted_2])}
        
        final_2 <- final_2 + t_2
      }
      final_3 <- matrix(0, nrow = 3, ncol = 3)
      for (experiment in experiment3){
        t_3 <- table(experiment$prediction, experiment$testsetpheno)
        final_3 <- final_3 + t_3
        wrongly_predicted_3 <- which(experiment$prediction != experiment$testsetpheno)
        if (length(wrongly_predicted_3) >0) {mistakes_3Classes <- c(mistakes_3Classes, experiment$testset[wrongly_predicted_3])}
        
        
      }
      
      sensitivity_2 <- final_2[1] / (final_2[1] + final_2[2])
      specificity_2 <- final_2[4] / (final_2[3] + final_2[4])
      
      
      sensitivity_3 <- final_3[1] / (final_3[1] + final_3[2] + final_3[3])
      specificity_3 <- (final_3[5] + final_3[6] + final_3[8] + final_3[9]) / (final_3[5] + final_3[6] + final_3[8] + final_3[9] + final_3[4] + final_3[7])
      
      precision_2 <- final_2[1] / (final_2[1] + final_2[3])
      recall_2 <- final_2[1] / (final_2[1] + final_2[2])
      
      f_2 <- 2*((precision_2*recall_2) / (precision_2 +recall_2))
      
      precision_3 <- final_3[1] / (final_3[1] + final_3[4] + final_3[7])
      recall_3 <- (final_3[1]) / (final_3[1] + final_3[2] + final_3[3])
      
      f_3 <- 2*((precision_3*recall_3) / (precision_3 +recall_3))
      
      balanced_2 <- ((final_2[1]/ (final_2[1] + final_2[2])) + (final_2[4]/ (final_2[3] + final_2[4])))/2
      balanced_3 <- ((final_3[1]/ (final_3[1] + final_3[2] + final_3[3])) + (final_3[5]/ (final_3[4] + final_3[5] + final_3[6])) + (final_3[9]/ (final_3[7] + final_3[8] + final_3[9])))/3
      
      results_2Classes <- data.frame("Tube" = tube, 
                                     "Type" = paste0("FlowSOM_", to_use_data), 
                                     "FS"= FS, 
                                     "Nr_Classes"=2, 
                                     "Specificity"= specificity_2, 
                                     "Sensitivity"= sensitivity_2,
                                     "Precision" = precision_2,
                                     "Recall" = recall_2,
                                     "F_measure" = f_2, 
                                     "Classifier" = name_classifier,
                                     "Balanced_accuracy" = balanced_2)
      results_2Classes$Mistakes <- list(mistakes_2Classes)
      results_3Classes <- data.frame("Tube" = tube, 
                                     "Type" = paste0("FlowSOM_", to_use_data), 
                                     "FS"= FS, 
                                     "Nr_Classes"=3, 
                                     "Specificity"= specificity_3, 
                                     "Sensitivity"= sensitivity_3,
                                     "Precision" = precision_3,
                                     "Recall" = recall_3,
                                     "F_measure" = f_3, 
                                     "Classifier" = name_classifier,
                                     "Balanced_accuracy" = balanced_3)  
      results_3Classes$Mistakes <- list(mistakes_3Classes)
      
      final_results <- rbind(final_results, results_3Classes, results_2Classes)
    }
  }
  }
}


final_results <- final_results[-1,]
final_results$Type <- as.factor(final_results$Type)
final_results$Tube <- as.factor(final_results$Tube)
final_results$FS <- as.factor(final_results$FS)
final_results$Nr_Classes <- as.factor(final_results$Nr_Classes)
final_results$Classifier <- as.factor(final_results$Classifier)

classes2 <- final_results[which(final_results$Nr_Classes == 2),]
classes3 <- final_results[which(final_results$Nr_Classes == 3),]



table(classes2$Mistakes %>% unlist) %>% sort
table(classes3$Mistakes %>% unlist) %>% sort


#2 classes
mean(classes2$Balanced_accuracy[which(classes2$Type == "Manual" & classes2$FS == "WithFS")])
sd(classes2$Balanced_accuracy[which(classes2$Type == "Manual" & classes2$FS == "WithFS")])

mean(classes2$Balanced_accuracy[which(classes2$Type != "Manual" & classes2$FS == "WithFS")])
sd(classes2$Balanced_accuracy[which(classes2$Type != "Manual" & classes2$FS == "WithFS")])

mean(classes2$Balanced_accuracy[which(final_results$Type == "Manual" & final_results$FS == "NoFS")])
sd(classes2$Balanced_accuracy[which(final_results$Type == "Manual" & final_results$FS == "NoFS")])

mean(classes2$Balanced_accuracy[which(classes2$Type != "Manual" & classes2$FS == "NoFS")])
sd(classes2$Balanced_accuracy[which(classes2$Type != "Manual" & classes2$FS == "NoFS")])

#3 classes
mean(classes3$Balanced_accuracy[which(classes3$Type == "Manual" & classes3$FS == "WithFS")])
sd(classes3$Balanced_accuracy[which(classes3$Type == "Manual" & classes3$FS == "WithFS")])

mean(classes3$Balanced_accuracy[which(classes3$Type != "Manual" & classes3$FS == "WithFS")])
sd(classes3$Balanced_accuracy[which(classes3$Type != "Manual" & classes3$FS == "WithFS")])

mean(classes3$Balanced_accuracy[which(final_results$Type == "Manual" & final_results$FS == "NoFS")])
sd(classes3$Balanced_accuracy[which(final_results$Type == "Manual" & final_results$FS == "NoFS")])

mean(classes3$Balanced_accuracy[which(classes3$Type != "Manual" & classes3$FS == "NoFS")])
sd(classes3$Balanced_accuracy[which(classes3$Type != "Manual" & classes3$FS == "NoFS")])



# -------------------- Make overview figure -----------------------------------
library(ggplot2)
library(reshape2)
library(cowplot)

new_results <- final_results
new_results.long <- melt(new_results, measure = c("Specificity", "Sensitivity"))


# p_1 <- ggplot(data = new_results.long, aes(x = Type, y = value, shape = Classifier, color = variable )) +
#   geom_point(size = 5, aes(fill = variable, alpha = FS))+ 
#   geom_point(size = 5) +
#   scale_shape_manual(values=c(21,22,23)) +
#   scale_alpha_manual(values=c("WithFS"=0, "NoFS"=1)) +
#   theme(panel.spacing.x = unit(0, "lines"))
# p_1 <- p_1 + theme(legend.title=element_blank()) 
# p_1 <- p_1 + facet_grid(Nr_Classes ~ Tube, labeller = labeller(Nr_Classes =  c("2" = "2 classes", "3" = "3 classes"), Tube = c("1" = "Panel 1", "2" = "Panel 2", "1_2" = "Panel 1 & 2")))
# p_1 <- p_1 + theme(axis.text.x = element_text(angle = 270, hjust = 0)) +ggtitle("Sensitivity/specificity")
# plot(p_1)
# 
# ggsave(
#   "Overview_Figures/sensitivity.png",
#   plot = p_1,
#   width = 12,
#   height = 12,
#   dpi = 300
# )

new_results$Type <- as.character(new_results$Type)
new_results$Type[which(new_results$Type == "FlowSOM_percentages_to")] <- "FlowSOM_percentages_clusters_to_metaclusters"
new_results$Type[which(new_results$Type == "FlowSOM_meta_clusters")] <- "FlowSOM_metaclusters"
new_results$Type <- factor(new_results$Type, levels = c("Manual", "FlowSOM_percentages_metaclusters", "FlowSOM_percentages_clusters", "FlowSOM_percentages_clusters_to_metaclusters", "FlowSOM_metaclusters", "FlowSOM_clusters", "FlowSOM_total"))
new_results$Tube <- factor(new_results$Tube, levels = c("1", "2", "5", "1_2_5"))

new_results_test <- new_results
new_results_test$FS <- factor(new_results$FS, levels(new_results$FS)[c(2,1)])
new_results_test <- new_results_test %>% arrange(FS)



colorado <- function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}


p_2 <- ggplot(data = new_results_test, aes(x = Type, y = Balanced_accuracy, shape = Classifier )) +
  
  geom_point(aes(shape=Classifier, fill=FS), size=5, alpha = 0.6) +
  scale_shape_manual(values=c(21,22), name = "Classifier", labels = c("Random Forest", "Support Vector Machine")) +
  scale_fill_manual(values=c(NoFS='orangered2',WithFS='blue4', None = "darkred"), name = "Feature Selection", labels = c("Without", "With", "Overlapping results"), limits = c("NoFS", "WithFS", "None")) + 
  
  #scale_alpha_manual(values=c("WithFS"=0, "NoFS"=1)) +
  theme(panel.spacing.x = unit(0, "lines"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
p_2 <- p_2 + theme(legend.title=element_blank())  + theme_bw()
p_2 <- p_2 + facet_grid(Nr_Classes ~ Tube, labeller = labeller(Nr_Classes =  c("2" = "2 classes", "3" = "3 classes"), Tube = c("1" = "Panel 1: PBMC's", "2" = "Panel 2: B cells", "5" = "Panel 3: T cells", "1_2_5" = "All panels combined")))
p_2 <- p_2 + theme(axis.text.x = element_text(angle = 270, hjust = 0, face = colorado(new_results_test$Type, "Manual")), axis.title.x = element_blank())+ ylab( "Balanced accuracy")
p_2 <- p_2 + geom_point(data = new_results_test[which(new_results_test$Type == "Manual"),], size = 4, stroke = 2, show.legend = FALSE, alpha = 0.6)

plot(p_2)

ggsave(
  "Overview_Figures/Balanced_accuracy_new.png",
  plot = p_2,
  width = 13,
  height = 10,
  dpi = 400
)



new_results_test_3 <- new_results_test[which(new_results_test$Nr_Classes == 3),]



p_3 <- ggplot(data = new_results_test_3, aes(x = Type, y = Balanced_accuracy, shape = Classifier )) +
  
  geom_point(aes(shape=Classifier, fill=FS), size=5, alpha = 0.6) +
  scale_shape_manual(values=c(21,22), name = "Classifier", labels = c("Random Forest", "Support Vector Machine")) +
  scale_fill_manual(values=c(NoFS='orangered2',WithFS='blue4', None = "darkred"), name = "Feature Selection", labels = c("Without", "With", "Overlapping results"), limits = c("NoFS", "WithFS", "None")) + 
  
  #scale_alpha_manual(values=c("WithFS"=0, "NoFS"=1)) +
  theme(panel.spacing.x = unit(0, "lines"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
p_3 <- p_3 + theme(legend.title=element_blank())  + theme_bw()
p_3 <- p_3 + facet_grid( ~ Tube, labeller = labeller(Tube = c("1" = "Panel 1: PBMC's", "2" = "Panel 2: B cells", "5" = "Panel 3: T cells", "1_2_5" = "All panels combined")))
p_3 <- p_3 + theme(axis.text.x = element_text(angle = 270, hjust = 0, face = colorado(new_results_test$Type, "Manual")), axis.title.x = element_blank())+ ylab( "Balanced accuracy")
p_3 <- p_3 + geom_point(data = new_results_test_3[which(new_results_test_3$Type == "Manual"),], size = 4, stroke = 2, show.legend = FALSE, alpha = 0.6)


plot(p_3)

ggsave(
  "Overview_Figures/Balanced_accuracy_new_three_classes.png",
  plot = p_3,
  width = 13,
  height = 6,
  dpi = 400
)


new_results_test_3_nofs <- new_results_test_3[which(new_results_test_3$FS == "NoFS"),]


p_3 <- ggplot(data = new_results_test_3_nofs, aes(x = Type, y = Balanced_accuracy, shape = Classifier )) +
  
  geom_point(aes(shape=Classifier, fill=FS), size=5, alpha = 0.6) +
  scale_shape_manual(values=c(21,22), name = "Classifier", labels = c("Random Forest", "Support Vector Machine")) +
  scale_fill_manual(values=c(NoFS='orangered2',WithFS='blue4', None = "darkred"), name = "Feature Selection", labels = c("Without", "With", "Overlapping results"), limits = c("NoFS", "WithFS", "None")) + 
  
  #scale_alpha_manual(values=c("WithFS"=0, "NoFS"=1)) +
  theme(panel.spacing.x = unit(0, "lines"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
p_3 <- p_3 + theme(legend.title=element_blank())  + theme_bw()
p_3 <- p_3 + facet_grid( ~ Tube, labeller = labeller(Tube = c("1" = "Panel 1: PBMC's", "2" = "Panel 2: B cells", "5" = "Panel 3: T cells", "1_2_5" = "All panels combined")))
p_3 <- p_3 + theme(axis.text.x = element_text(angle = 270, hjust = 0, face = colorado(new_results_test$Type, "Manual")), axis.title.x = element_blank())+ ylab( "Balanced accuracy")
p_3 <- p_3 + geom_point(data = new_results_test_3_nofs[which(new_results_test_3_nofs$Type == "Manual"),], size = 4, stroke = 2, show.legend = FALSE, alpha = 0.6)


plot(p_3)

ggsave(
  "Overview_Figures/Balanced_accuracy_new_three_classes_noFS.png",
  plot = p_3,
  width = 13,
  height = 6,
  dpi = 400
)


