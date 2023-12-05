#' Resulting clustergram/heatmap.
#' Function to display a categorical data matrix given a user defined number of clusters `nCl`, a categorical distance `distName` and a predefined clustering method `FUNcluster`.
#' The output displays a heatmap separating and color-labelling resulting clusters vertically in the rows and allowing unsupervised clustering on questions in the columns. Each cell is colored according to the categorical values provided or found in the data.
#' The clustergram is based on the `pheatmap` function from the pheatmap R package. Thus, any parameter found in pheatmap can be specified to `clusGapDiscrHeat`.
#' This function can be used to examine number of clusters before running `clusGapDiscrHeat` but also after number of clusters is determined.
#' @param x matrix object or data.frame
#' @param nCl number of clusters to plot
#' @param distName Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This function returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' @param catVals All categorical values to consider for the plot. If `NULL` the function will extract all unique categorical values in `x`.
#' @return pheatmap output.
#' @export

clusGapDiscrHeat <- function(x,
                       FUNcluster,
                       nCl,
                       distName,
                       catVals=NULL,
                       ...){

   if(!is.matrix(x))
      x <- as.matrix(x)
   AssignedCls <- FUNcluster(distancematrix(x, d = distName),
                      k = nCl)$cluster

   x <- data.frame(rowNames = row.names(x),
                   Clust = paste0('Cl', AssignedCls),
                   x,
                   check.names=FALSE) %>%
      dplyr::arrange(Clust, rowNames)
   x$rowNames <- c()
   rGaps <- with(x, table(Clust)) %>% as.vector

   if(length(rGaps) > 1 )
      rGaps <- rGaps[1:(length(rGaps)-1)] %>% cumsum
   else
      rGaps <- NULL

   clusCols = list(Clust = stats::setNames(nm=paste0('Cl', 1:nCl),
                                    object = Polychrome::palette36[1:nCl]))

   myColors = stats::setNames(RColorBrewer::brewer.pal(n = length(catVals),
                                                name = 'Greens'), catVals)
   pheatmap::pheatmap(mat = x[, -1], color = myColors,
                      border.color=NULL, cluster_cols=TRUE,
                      cluster_rows=FALSE,
                      annotation_colors = clusCols,
                      annotation_row = subset(x, select = 'Clust'),
                      gaps_row = rGaps,
                      ...)
}

#' Multidimensional scaling plot for categorical data.
#' Scatter dot dimension reduction plot obtained from user input number of clusters `nCl`, specified categorical distance `distName` and clustering method function `FUNcluster.
#' First two dimensions are used including the estimated percentage explained variability. Dots are color-coded according to the assigned cluster.
#' The function outputs a ggplot object where additional layers can be added for customization purposes.
#' @param x matrix object or data.frame
#' @param nCl number of clusters to plot
#' @param distName Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This function returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' @param title Plot title
#' @return ggplot object.
#' @export
clusGapDiscrMDS <- function(x,
                           distName,
                           nCl,
                           FUNcluster,
                           title=NULL){

   if(is.null(title))
      title <- paste0('Discrete Gap Statistic\nDistance: ', distName,
                      '; nClusters=', nCl)
   if(!is.matrix(x))
      x <- as.matrix(x)
   AssignedCls <- FUNcluster(distancematrix(x, d = distName),
                             k = nCl)$cluster

   distData <- x %>%
      distancematrix(d = distName)
   myMDS <- stats::cmdscale(d = distData, eig=TRUE)

   PerVarExpl <- round(100*c(myMDS$eig^2)/sum(myMDS$eig^2), 2)
   myMDS <- data.frame(myMDS$points)
   colnames(myMDS) <- c("Dim1", "Dim2")
   PerVarExpl <- paste0(c("Dim1", "Dim2"), ' (', PerVarExpl[1:nCl], '%)')

   myClsCols = stats::setNames(nm=paste0('Cl', 1:nCl), object = Polychrome::palette36.colors[1:nCl])
   ggOut <- ggplot2::ggplot(x = myMDS %>%
                               dplyr::mutate(Clust =  paste0('Cl', AssignedCls)),
                            aes(x = Dim1, y = Dim2, fill = Clust)) +
      ggplot2::geom_jitter(size = 2, shape = 21, alpha = 0.8, width =0.01, height=0.01)+
      ggplot2::scale_fill_manual(values = myClsCols)+
      ggplot2::theme_bw()+
      ggplot2::theme(text = ggplot2::element_text(size = 13), legend.position = 'bottom')+
      ggplot2::labs(title = title, x = PerVarExpl[1], y = PerVarExpl[2]) +
      ggplot2::scale_color_manual(values = myClsCols)+
      ggplot2::geom_rug(ggplot2::aes(color=Clust), alpha = 0.5) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 4)),
             color = 'none')
   ggOut
   ## ggsave(filename = paste0(outDir, '/', fileName, '.', fileFormat),
   ##        plot = ggOut, height = 6, width= 7 )
}
