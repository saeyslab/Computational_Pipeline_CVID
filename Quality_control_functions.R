# Function to remove events on the margins
removeMargins <- function(flowFrame,dimensions){
  # Look up the accepted ranges for the dimensions
  meta <- pData(flowFrame@parameters)
  rownames(meta) <- meta[,"name"]
  
  # Initialize variables
  selection <- rep(TRUE,times=dim(flowFrame)[1])
  e <- flowFrame@exprs
  
  # Make selection
  for(d in dimensions){
    selection <- selection & 
      e[,d] > max(min(meta[d,"minRange"],0),min(e[,d])) & ## first: the min between the minRange and 0 (minRange could be wrong), than the max between this and the min found in the expression matrix
      e[,d] < min(meta[d,"maxRange"],max(e[,d]))
  }
  return(selection)
}

gating_subset_adapted <- function(flowjo_res, gate){

  if(!is.null(flowjo_res$flowSet)){
    if (length(gate) > 1){
      res <- lapply(seq_len(length(flowjo_res$flowSet)),
                    function(i){
                      inbetween <- flowjo_res$gates[[i]][,gate]
                      flowjo_res$flowSet[[i]][which(apply(inbetween, 1, any)),]
                      
                    })
      
      names(res) <- sampleNames(flowjo_res$flowSet)
      return(list(
        flowSet = flowCore::flowSet(res),
        gates = lapply(flowjo_res$gates, function(x){
          inbetween <- x[,gate]
          x[which(apply(inbetween,1,any)),]})
      ))
    } else {
      res <- lapply(seq_len(length(flowjo_res$flowSet)),
                    function(i){
                      flowjo_res$flowSet[[i]][flowjo_res$gates[[i]][,gate],]
                    })
      names(res) <- sampleNames(flowjo_res$flowSet)
      return(list(
        flowSet = flowCore::flowSet(res),
        gates = lapply(flowjo_res$gates, function(x)x[x[,gate], ])
      ))
    }
    
  } else {
    return(list(flowFrame = flowjo_res$flowFrame[flowjo_res$gates[,gate], ],
                gates = flowjo_res$gates[flowjo_res$gates[,gate], ]))
  }
}

parse_flowjo_compensation <- function(files,
                         wsp_file,
                         experiment,
                         fcsDir,
                         tube, 
                         group = "All Samples",
                         plot = FALSE,
                         html_plot = "_QC") {
  
  # Load in FlowJo workspace
  wsp <- flowWorkspace::openWorkspace(wsp_file)
  o <- capture.output(
    gates <- suppressMessages(flowWorkspace::parseWorkspace(wsp, group))
  )
  files_in_wsp <- gates@data@origSampleVector
  counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
  result <- list()
  file <- files[11]
  for(file in files){
    print(paste0("Processing ", file))
    file_id <- grep(gsub(".*/", "", file), files_in_wsp)
    if(length(file_id) == 0) {stop("File not found. Files available: ",
                                   gsub("_[0-9]*$", "\n", files_in_wsp))}
    gate_names <- flowWorkspace::getNodes(gates[[file_id]], path = "auto")
    gatingMatrix <- matrix(FALSE,
                           nrow = counts[file_id],
                           ncol = length(gate_names),
                           dimnames = list(NULL, gate_names))
    for (gate in gate_names) {
      gatingMatrix[, gate] <- flowWorkspace::getIndiceMat(gates[[file_id]],
                                                          gate)
    }
    
    
    if (tube %in% c(1,2)){
      if (experiment %in% c(13,14,15)){
        ff <- read.FCS(file.path(fcsDir, "Tube 1-2-6", file))
      } else {
        ff <- read.FCS(file.path(fcsDir,"CVID plate 1-2-6", file))
      }
    } else if (tube == 5){
      if (experiment %in% c(13,14,15)){
        ff <- read.FCS(file.path(fcsDir, "Tube 5", file))
      } else {
        ff <- read.FCS(file.path(fcsDir,"CVID plate 5", file))
      }
    }
    
    # Remove NA files
    selection_NA <- complete.cases(ff@exprs)
    ff <- ff[complete.cases(ff@exprs),]
    
    # Remove events that have measures outside of the detection range
    selection_margins <- removeMargins(ff, dimensions = colnames(ff))
    ff <- ff[selection_margins, ]
    
    # Compensation with compensation matrix of workspace
    compensation <- getCompensationMatrices(gates[[1]])
    
    ff <- compensate(ff, compensation)
    
    # Transformation
    ff <- transform(ff, estimateLogicle(ff, colnames(compensation@spillover)))

    # Quality Control
    ff_AI <- flow_auto_qc(ff, folder_results =  paste0("Preprocessing/Tube",tube,"/Exp_",experiment,"/QC/"), 
                          output =2, outlier_binsFS = TRUE, fcs_QC = FALSE, html_report = html_plot)
    
    selected_flowAI <- ff_AI@exprs[,"QCvector"] <= 10000
    
    ff <- ff[selected_flowAI, ]
    
    gatingMatrix <- gatingMatrix[selection_NA,]
    gatingMatrix <- gatingMatrix[selection_margins,]
    gatingMatrix <- gatingMatrix[selected_flowAI,]
    
    result[[file]] <- list("flowFrame" = ff,
                           "gates" = gatingMatrix)
    
    if (plot) {
      flowWorkspace::plot(gates[[file_id]])
    }
  }
  if (length(files) == 1){
    result <- result[[1]]
  } else {
    result <- list(flowSet = flowCore::flowSet(lapply(result, function(x) x$flowFrame)),
                   gates = lapply(result, function(x) x$gates))
  }
  return(result)
}


