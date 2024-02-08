#' Summary Heatmap for categorical/Likert data
#' Heatmap representation summarizing categorical/likert data.
#' Modified version of `likert.heat.plot` from `likert` package.
#' Does not allow different categorical ranges across questions.
#' The function outputs a ggplot object where additional layers can be added for customization purposes.
#' The output plot preserves the question order given by columns of `x`.
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import dplyr
#' @import likert
#' @param x matrix object or data.frame with categorical data. Columns are questions and rows are observations.
#' @param allLevels vector with all categorical (ordered) levels.
#' @param low.color string; name of color assigned to the first level found in `allLevels`.
#' @param high.color string; name of color assigned to the last level found in `allLevels`.
#' @param text.color string; text color of numbers within cells.
#' @param text.size string; text size for numbers within cells.
#' @param textLen string; maximum length of text-length for question labels (column names)
#' @param ... other valid arguments in pheatmap function
#' @return ggplot object.
#' @export
likert.heat.plot2 <- function(x,
                              allLevels,
                              low.color = "white",
                              high.color = "blue",
                              text.color = "black",
                              text.size = 4,
                              textLen = 50,
                              ...) {

   Item <- variable <- value <- label <- NULL
   xLikertForm <- x %>%
      apply(2, as.character) %>%
      data.frame(check.names=FALSE) %>%
      mutate_if(is.character, function(x) factor(x, levels = allLevels))

   likertOut <- likert(xLikertForm)
   lsum <- summary(object=likertOut, ordered = FALSE)

   results <- melt(likertOut$results, id.vars = "Item")
   results$variable <- as.character(results$variable)
   results$label <- paste(format(results$value,
                                 digits = 2,
                                 drop0trailing = FALSE),
                          "%", sep = "")
   results$Item <- factor(x = results$Item,
                          levels = lsum$Item %>% rev)

   p <- ggplot(results,
               aes(x = Item,
                   y = variable,
                   fill = value,
                   label = label)) +
      scale_y_discrete(limits = names(likertOut$results)[2:ncol(likertOut$results)]) +
      geom_tile() +
      geom_text(size = text.size, colour = text.color) +
      coord_flip() +
      scale_fill_gradient2("Percent", low = "white", mid = low.color,
                           high = high.color, limits = c(0, 100)) + xlab("") +
      ylab("") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) +
      scale_x_discrete(breaks = likertOut$results$Item,
                       labels = base::substring(likertOut$results$Item, first=1, last = textLen))
   # class(p) <- c("likert.heat.plot", class(p))
   return(p)
}

#' sample-to-sample heatmap clustering samples according to a given categorical distance
#' Exploratory tool that helps to visualize/cluster blocks of observations across columns ordered according to given categorical distance. The final output is a clustered distance matrix.
#' This plot is aimed to guide the `DiscreteClusGap` user to give an idea which type of categorical distance would accommodate better to the inputted data.
#' `sample2sampleHeat` is based on the `pheatmap` function from the `pheatmap` R package. Thus, any parameter found in pheatmap can be specified to `sample2sampleHeat`.
#' @importFrom magrittr %>%
#' @import pheatmap
#' @param x matrix object or data.frame
#' @param distName Name of categorical distance to apply.
#' @param ... other valid arguments in pheatmap function
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param clustering_method string; clustering method used by pheatmap
#' @return clustered heatmap
#' @export
distanceHeat <- function(x,
                        distName,
                        clustering_method = 'complete',
                        ...){
   x0 <- x
   if(!is.matrix(x))
      x <- as.matrix(x)

   myDist <- distancematrix(X = x, d = distName)

   distMatr <- as.matrix(myDist)
   colnames(distMatr) <- rownames(x)
   rownames(distMatr) <- rownames(x)

   hmcol <- grDevices::colorRampPalette(colors = rev(RColorBrewer::brewer.pal(8, "Blues")))(87)
   pheatmap(mat = distMatr,
                      color = hmcol,
                      clustering_distance_rows = myDist,
                      clustering_distance_cols = myDist,
                      clustering_method = clustering_method,
                      border_color=NA,
                      treeheight_row=0,
                      ...)

}

#' Heatmap assuming a given a distance function and a known number of clusters.
#' Function to display a categorical data matrix given a user defined number of clusters `nCl`, a categorical distance `distName` and a predefined clustering method `FUNcluster`.
#' The output displays a heatmap separating and color-labelling resulting clusters vertically in the rows and allowing unsupervised clustering on questions in the columns. Each cell is colored according to the categorical values provided or found in the data.
#' The clustergram is based on the `pheatmap` function from the pheatmap R package. Thus, any parameter found in pheatmap can be specified to `clusGapDiscrHeat`.
#' This function can be used to examine number of clusters before running `clusGapDiscrHeat` but also after number of clusters is determined.
#' @importFrom magrittr %>%
#' @import pheatmap
#' @import Polychrome
#' @import RColorBrewer
#' @param x matrix object or data.frame
#' @param nCl number of clusters to plot; if `nCl` is a vector, its length lN is the number of clusters, with values ranging from 1 to lN in desired order to place the clusters.
#' @param distName Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param catVals character string vector with (ordered) categorical values
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This function returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' @param out Specifies the desired output between "heatmap" (default; produce a heatmap), "clusters" (return a `data.frame` with clustering assignments) or "clustersReord" (return a `data.frame` with reorganized clusters)
#' @param clusterNames Either `null` or 'renumber'. When `nCl` the cluster ordering is rearranged. `NULL` leaves cluster names as their original cluster assignment. 'renumber' respects the rearrangements but relabels the cluster numbers from top to bottom in ascending order.
#' @param prefObs character string vector of length 1 with a prefix for the observations, in case they come unlabelled or the user wants to anomymize sample IDs
#' @param filename character string with name of file output
#' @param outDir character string with the directory path to save output file
#' @param height numeric height of output plot in inches
#' @param height numeric width of output plot in inches
#' @return png file or ComplexHeatmap object
#' @export
ResHeatmap <- function(x,
                       nCl,
                       distName,
                       catVals,
                       FUNcluster = cluster::pam,
                       out = 'heatmap',
                       clusterNames = NULL,
                       prefObs = NULL,
                       filename = NULL,
                       outDir = NULL,
                       height = 10, width = 6){
   ## To Do: specify option where the input data is plotted without any cluster rearrangement

   if(length(nCl) > 1){
      clOrd <- nCl
      nCl <- length(nCl)
   }else{
      clOrd <- NULL
   }

   if(!is.matrix(x))
      x <- as.matrix(x)

   if(nCl > 0){
      AssignedCls <- FUNcluster(x = distancematrix(x, d = distName),
                                k = nCl)$cluster
   }else{
      AssignedCls <- rep(0, nrow(x))
   }

   if(!is.null(prefObs))
      row.names(data) <- paste0(prefObs, 1:nrow(x))

   ## Notice that data is reorganized according to clusters!
   data <- data.frame(rowNames = row.names(x),
                      Clust = AssignedCls,
                      x,
                      check.names=FALSE) %>%
      arrange(Clust, rowNames)

   if(out == 'clusters')
      return(data[, 1:2])

   if(!is.null(clOrd)){
      ## Alter cluster ordering
      data$Clust <- clOrd[data$Clust]

      if(!is.null(clusterNames)){
         if(clusterNames == 'renumber'){
            data <- data %>% arrange(Clust, rowNames)
            ClustRun <- rle(data$Clust)$lengths
            data$Clust <- rep(1:length(ClustRun), ClustRun)
         }
      }
   }

   if(out == 'clustersReord')
      return(data[, 1:2])

   ## Tones of green
   myCols = stats::setNames(object = RColorBrewer::brewer.pal(n = length(catVals),
                                                              name = 'Greens'),
                            nm = catVals)
   myColsCH <- structure(myCols, names = names(myCols))

   rowSplit <- subset(data, select = 'Clust')
   rowTitle <- paste0('Cluster ', unique(data$Clust))

   if(!is.null(outDir) & !is.null(filename) ){
      grDevices::png(paste0(outDir, '/', filename,'.png'),
                     width = width, height = height, units = 'in', res=500)

      hmObj <- ComplexHeatmap::Heatmap(mat = data[, -c(1:2)],
                                       col = myCols,
                                       name = ' ',
                                       show_row_names = FALSE,
                                       # left_annotation = rowAnn,
                                       show_heatmap_legend=c(TRUE) ,
                                       # row_split = data[, 'Clust'],
                                       row_split = rowSplit,
                                       row_title = rowTitle)
      ComplexHeatmap::draw(hmObj,
                           column_title = paste0("DiscreteClusGap Assignments\n", distName),
                           column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                           show_annotation_legend=FALSE)
      grDevices::dev.off()

   }else{
      hmObj <-  ComplexHeatmap::Heatmap(mat = data[, -c(1:2)],
                                        col = myCols,
                                        name = ' ',
                                        show_row_names = FALSE,
                                        ## left_annotation = rowAnn,
                                        ## row_split = data[, 'Clust'],
                                        row_split = rowSplit,
                                        show_heatmap_legend=c(TRUE) ,
                                        row_title = ' ')
      ComplexHeatmap::draw(hmObj,
                           column_title = paste0("DiscreteClusGap Assignments\n", distName),
                           column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                           show_annotation_legend=FALSE)
   }
}

#' Multidimensional scaling plot for categorical data.
#' Scatter dot dimension reduction plot obtained from user input number of clusters `nCl`, specified categorical distance `distName` and clustering method function `FUNcluster.
#' First two dimensions are used including the estimated percentage explained variability. Dots are color-coded according to the assigned cluster.
#' The function outputs a ggplot object where additional layers can be added for customization purposes.
#' @import ggplot2
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

   Dim1 <- Dim2 <- Clust <- NULL
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

   myClsCols = stats::setNames(nm=paste0('Cl', 1:nCl), object = Polychrome::palette36.colors()[1:nCl])
   ggOut <- ggplot(data = myMDS %>%
                               dplyr::mutate(Clust =  paste0('Cl', AssignedCls)),
                            aes(x = Dim1, y = Dim2, fill = Clust)) +
      geom_jitter(size = 2, shape = 21, alpha = 0.8, width =0.01, height=0.01)+
      scale_fill_manual(values = myClsCols)+
      theme_bw()+
      theme(text = element_text(size = 13), legend.position = 'bottom')+
      labs(title = title, x = PerVarExpl[1], y = PerVarExpl[2]) +
      scale_color_manual(values = myClsCols)+
      geom_rug(aes(color=Clust), alpha = 0.5) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 4)),
             color = 'none')
   ggOut
   ## ggsave(filename = paste0(outDir, '/', fileName, '.', fileFormat),
   ##        plot = ggOut, height = 6, width= 7 )
}

