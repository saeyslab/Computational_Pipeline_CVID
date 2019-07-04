PlotClusters1D <- function (fsom, markerlist, metaclusters, 
                            col = "#FF0000", 
                            maxBgPoints = 10000,
                            main = "",
                            textsize = 11,
                            ff, ...) 
{

  
#  abline(v = fsom$map$medianValues[nodes, marker], col= "blue")
  


#metaclusters <- c(1,14, 15)
#markers <- c("APC-A", "FITC-A")

df_markers <- NULL
df_vline <- NULL

marker <- markerlist[1]
for (marker in markerlist) {
  
  metacluster <- metaclusters[1]
  for (metacluster in metaclusters ){
    forground_data <- data.frame(value = fsom$FlowSOM$data[fsom$FlowSOM$map$mapping[, 1] %in% which(fsom$metaclustering == metacluster), 
                               marker], 
                               marker = fsom$FlowSOM$prettyColnames[marker], 
                               bg = "fg",
                              metacluster = paste0("metacluster ",metacluster),
                               id = seq_len(length(fsom$FlowSOM$data[fsom$FlowSOM$map$mapping[, 1] %in% which(fsom$metaclustering == metacluster), 
                                                   marker])))
    
    if (!is.null(maxBgPoints)) {
      background <- sample(seq_len(nrow(fsom$FlowSOM$data)), min(maxBgPoints, 
                                                         nrow(fsom$FlowSOM$data)))
    }
    else {
      background <- seq_len(nrow(fsom$FlowSOM$data))
    }
    
    background_data <- data.frame(value =fsom$FlowSOM$data[background, marker], 
                                  marker = fsom$FlowSOM$prettyColnames[marker], 
                                  bg = "bg",
                                  metacluster = paste0("metacluster ",metacluster),
                                  id = seq_len(length(background)))
    
    
    #dv_vline_specific <- data.frame(vl = list(fsom$FlowSOM$map$medianValues[fsom$metaclustering == metacluster, marker]), 
    #                                marker = fsom$FlowSOM$prettyColnames[marker], 
    #                                metacluster = paste0("metacluster ", metacluster) )
    
    df_markers <- rbind(df_markers, forground_data, background_data)
    #df_vline <- rbind(df_vline, dv_vline_specific)
  }
}

df_markers$marker <- as.factor(df_markers$marker)
df_markers$metacluster <- as.factor(df_markers$metacluster)
df_markers$bg <- factor(df_markers$bg, levels = c("fg", "bg"))



p <- ggplot(data = df_markers, aes(x=value)) + geom_density(aes(fill=bg), alpha = 0.4)
p <- p + facet_grid( marker ~ metacluster)
p <- p + scale_fill_brewer(palette = "Set1", labels = c("Cluster density", "background density"), name = "") 
p <- p + theme(strip.text.y = element_text(size = textsize))
p <- p + xlim(c(-10,10))

plot(p)

}


















