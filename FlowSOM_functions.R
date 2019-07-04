PlotMarkerMFI <- function(files, 
                          markers_to_plot,
                          ff,
                          marker_names,
                          figure_dir) {
  suppressWarnings(dir.create(figure_dir))
  counts <- table(gsub("/.*","",files))
  idc <- counts[1]
  start <- counts[1]
  for (element in counts[-1]){
    start <- start + element
    idc <- c(idc , start)
  }
  marker <- markers_to_plot[1]
  for(marker in markers_to_plot){
    print(paste0("Plotting ",marker_names[marker]," for the aggregated file."))
    name_marker <- gsub("/", "_", marker_names[marker])
    png(file.path(figure_dir,
                  paste0("Aggregate_",name_marker,".png")), width = 1200, height = 1000)
    plot(ff@exprs[,c("File_scattered",marker)],
         pch=".",col="#00000044",
         xlab="Files",
         ylab=marker_names[marker])
    abline(v = idc,col = "red")
    mtext(text =c(13:33), side = 3, cex = 1.3 , at = idc-5)
    dev.off()
  }
}


PlotStars_adapted <- function (fsom, markers = fsom$map$colsUsed, view = "MST", colorPalette = grDevices::colorRampPalette(c("#00007F", 
                                                                                                        "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
                                                                                                        "red", "#7F0000")), starBg = "white", backgroundValues = NULL, 
          backgroundColor = function(n) {
            grDevices::rainbow(n, alpha = 0.3)
          }, backgroundLim = NULL, backgroundBreaks = NULL, backgroundSize = NULL, 
          thresholds = NULL, scaling = NULL, legend = TRUE, query = NULL, main = "") 
{
  add.vertex.shape("star", clip = igraph.shape.noclip, plot = FlowSOM:::mystar, 
                   parameters = list(vertex.data = NULL, vertex.cP = colorPalette, 
                                     vertex.scale = TRUE, vertex.bg = starBg))
  if (is.null(thresholds) & is.null(scaling)) {
    data <- fsom$map$medianValues[, markers, drop = FALSE]
    scale <- TRUE
  } else {
    if (fsom$transform) {
      warning("Thresholds should be given in the transformed space")
    }
    if (scaling == "PosPop"){
      if (!is.null(fsom$scaled.center) & scaling== "Original") {
        thresholds <- scale(t(thresholds), center = fsom$scaled.center[markers], 
                            scale = fsom$scaled.scale[markers])
      }
      data <- t(sapply(seq_len(fsom$map$nNodes), function(i) {
        res = NULL
        for (m in seq_along(markers)) {
          res = c(res, sum(subset(fsom$data, fsom$map$mapping[, 
                                                              1] == i)[, markers[m]] > thresholds[m])/sum(fsom$map$mapping[, 
                                                                                                                           1] == i))
        }
        res
      })) } 
    if (scaling == "MaxValue"){
      if (length(grep("FSC|SSC", colnames(fsom$data)[markers], fixed = FALSE)) != 0){
        
        scatter_idc <- grep("FSC|SSC", colnames(fsom$data)[markers])  
        used_scatter_idc <- markers[scatter_idc]
        other_idc <- markers[-scatter_idc]
        
        
        } else { other_idc <- markers}  
        
        
        max_marker <- which.max(apply(fsom$data[,other_idc],2,max))
        min_marker <- which.min(apply(fsom$data[,other_idc],2,min))
        current_min <- (min(apply(fsom$data[,other_idc],2,min)) * fsom$scaled.scale[other_idc][min_marker]) + fsom$scaled.center[other_idc][min_marker]
        current_max <- (max(apply(fsom$data[,other_idc],2,max)) * fsom$scaled.scale[other_idc][max_marker]) + fsom$scaled.center[other_idc][max_marker]
        
        data <- sapply(other_idc, function(i){
          (fsom$map$medianValues[, i] * fsom$scaled.scale[i]) + fsom$scaled.center[i]
        })
        
  
        data_2 <- sapply(seq_len(ncol(data)), function(i){
          (data[,i] - current_min ) / current_max
        })
        
        if (length(grep("FSC|SSC", colnames(fsom$data)[markers], fixed = FALSE)) != 0){
          
          FSC_idc <- grep("FSC", colnames(fsom$data)[markers])
          data_FSC <- (fsom$map$medianValues[, FSC_idc] *fsom$scaled.scale[FSC_idc]) + fsom$scaled.center[FSC_idc]
          min_FSC <- (min(fsom$data[,FSC_idc]) * fsom$scaled.scale[FSC_idc])+ fsom$scaled.center[FSC_idc] 
          max_FSC <- (max(fsom$data[,FSC_idc]) * fsom$scaled.scale[FSC_idc])+ fsom$scaled.center[FSC_idc]
          data_2_FSC <- (data_FSC - min_FSC)/ max_FSC
          
          SSC_idc <- grep("SSC", colnames(fsom$data)[markers])
          data_SSC <- (fsom$map$medianValues[, SSC_idc] *fsom$scaled.scale[SSC_idc]) + fsom$scaled.center[SSC_idc]
          min_SSC <- (min(fsom$data[,SSC_idc]) * fsom$scaled.scale[SSC_idc])+ fsom$scaled.center[SSC_idc] 
          max_SSC <- (max(fsom$data[,SSC_idc]) * fsom$scaled.scale[SSC_idc])+ fsom$scaled.center[SSC_idc]
          data_2_SSC <- (data_SSC - min_SSC)/ max_SSC
          
          
          data_2 <- cbind(data_2_FSC, data_2_SSC, data_2)
        } 
        
        data <- data_2
      
      
      }
      

    scale <- FALSE
  }
  switch(view, MST = {
    layout <- fsom$MST$l
    lty <- 1
  }, grid = {
    layout <- as.matrix(fsom$map$grid)
    lty <- 0
  }, tSNE = {
    layout <- fsom$MST$l2
    lty <- 0
  }, stop("The view should be MST, grid or tSNE. tSNE will only work if you specified this when building the MST."))
  if (!is.null(backgroundValues)) {
    background <- FlowSOM:::computeBackgroundColor(backgroundValues, 
                                         backgroundColor, backgroundLim, backgroundBreaks)
    if (is.null(backgroundSize)) {
      backgroundSize <- fsom$MST$size
      backgroundSize[backgroundSize == 0] <- 3
    }
  } else {
    background <- NULL
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(1, 1, 1, 1))
  if (legend) {
    if (!is.null(backgroundValues)) {
      graphics::layout(matrix(c(1, 3, 2, 3), 2, 2, byrow = TRUE), 
                       widths = c(1, 2), heights = c(1))
    } else {
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), 
                       widths = c(1, 2), heights = c(1))
    }
    if (is.null(query)) {
      FlowSOM:::plotStarLegend(fsom$prettyColnames[markers], colorPalette(ncol(data)))
    } else {
      FlowSOM:::plotStarQuery(fsom$prettyColnames[markers], values = query == 
                      "high", colorPalette(ncol(data)))
    }
    if (!is.null(backgroundValues)) {
      FlowSOM:::PlotBackgroundLegend(backgroundValues, background)
    }
  }
  igraph::plot.igraph(fsom$MST$g, vertex.shape = "star", vertex.label = NA, 
                      vertex.size = fsom$MST$size, vertex.data = data, vertex.cP = colorPalette(ncol(data)), 
                      vertex.scale = scale, layout = layout, edge.lty = lty, 
                      mark.groups = background$groups, mark.col = background$col[background$values], 
                      mark.border = background$col[background$values], mark.expand = backgroundSize, 
                      main = main)
  graphics::par(oldpar)
  graphics::layout(1)
}
