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
   tmp <- data.frame(Item = lsum$Item,
                    variable = rep("Mean (SD)",
                                   nrow(lsum)),
                    value = rep(-100, nrow(lsum)),
                    label = paste(format(lsum$mean,
                                         digits = 3, drop0trailing = FALSE),
                                  " (", format(lsum$sd,
                                               digits = 2, drop0trailing = FALSE), ")", sep = ""),
                    stringsAsFactors = FALSE)
   results <- rbind(tmp, results)
   results$Item <- factor(x = results$Item, levels = lsum$Item %>% rev)

   p <- ggplot(results,
              aes(x = Item,
                  y = variable,
                  fill = value,
                  label = label)) +
      scale_y_discrete(limits = c("Mean (SD)",
                                  names(likertOut$results)[2:ncol(likertOut$results)])) +
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

#' Clustergram/heatmap assuming a given number of clusters.
#' Function to display a categorical data matrix given a user defined number of clusters `nCl`, a categorical distance `distName` and a predefined clustering method `FUNcluster`.
#' The output displays a heatmap separating and color-labelling resulting clusters vertically in the rows and allowing unsupervised clustering on questions in the columns. Each cell is colored according to the categorical values provided or found in the data.
#' The clustergram is based on the `pheatmap` function from the pheatmap R package. Thus, any parameter found in pheatmap can be specified to `clusGapDiscrHeat`.
#' This function can be used to examine number of clusters before running `clusGapDiscrHeat` but also after number of clusters is determined.
#' @importFrom magrittr %>%
#' @import pheatmap
#' @import Polychrome
#' @import RColorBrewer
#' @param x matrix object or data.frame
#' @param nCl number of clusters to plot
#' @param distName Name of categorical distance to apply.
#' @param ... other valid arguments in pheatmap function
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This function returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' @param catVals All categorical values to consider for the plot. In the case `x` is a character matrix, the vector should specify category-number mapping. If `NULL` the function will extract all unique categorical values in `x` and assign numbers to them alphabetically.
#' @return pheatmap output.
#' @export
clusGapDiscrHeat <- function(x,
                       FUNcluster,
                       nCl,
                       distName,
                       catVals=NULL,
                       ...){

   Clust <- rowNames <- NULL
   x0 <- x
   if(!is.matrix(x))
      x <- as.matrix(x)

   myUniq <- as.vector(x) %>% unique %>% sort

   if(is.character(myUniq)){
      if(is.null(catVals) ){
         catVals <- stats::setNames(0:(length(myUniq)-1), myUniq)
      }else{
        stopifnot(exprs = length(myUniq) == length(catVals))
      }

      # rownames(x) <- paste0('s', 1:nrow(x))
      xNum <- catVals[x] %>% matrix(byrow=FALSE, nrow = nrow(x), ncol=ncol(x) )
      rownames(xNum) <- rownames(x)
      x <- xNum

      message('Since pheatmap can only generate plots from numerical matrices, ')
      message('categories in the matrix will be transformed to numerical values:')
      message(paste0(names(catVals), ' -> ', catVals, collapse = ', '))
   }

   AssignedCls <- FUNcluster(distancematrix(x, d = distName),
                      k = nCl)$cluster

   x <- data.frame(rowNames = row.names(x),
                   Clust = paste0('Cl', AssignedCls),
                   x,
                   check.names = FALSE) %>%
      dplyr::arrange(Clust, rowNames)
   x$rowNames <- c()
   rGaps <- with(x, table(Clust)) %>% as.vector

   if(length(rGaps) > 1 )
      rGaps <- rGaps[1:(length(rGaps)-1)] %>% cumsum
   else
      rGaps <- NULL

   clusCols = list(Clust = stats::setNames(nm=paste0('Cl', 1:nCl),
                                    object = Polychrome::palette36.colors()[1:nCl]))

   myColors = stats::setNames(brewer.pal(n = length(catVals),
                              name = 'Greens'),
                              catVals)
   pheatmap(mat = x[, -1],
            color = myColors,
            border.color=NULL,
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

