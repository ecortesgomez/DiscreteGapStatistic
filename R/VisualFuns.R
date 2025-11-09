utils::globalVariables(c('Dim1', 'Dim2', 'rowNames'))

#' Summary Heatmap for categorical data
#'
#' Heatmap representation summarizing categorical/likert data.
#' Modified version of `likert.heat.plot` from `likert` package.
#' Does not allow different categorical ranges across questions.
#' The function outputs a ggplot object where additional layers can be added for customization purposes.
#' The output plot preserves the question order given by columns of `x`.
#'
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @param x matrix object or data.frame with categorical data. Columns are questions and rows are observations.
#' @param allLevels vector with all categorical (ordered) levels.
#' @param low.color string; name of color assigned to the first level found in `allLevels`.
#' @param high.color string; name of color assigned to the last level found in `allLevels`.
#' @param text.color string; text color of numbers within cells.
#' @param text.size string; text size for numbers within cells.
#' @param textLen string; maximum length of text-length for question labels (column names)
#' @return ggplot object.
#' @export
likert.heat.plot2 <- function (x,
                              allLevels,
                              low.color = "white",
                              high.color = "blue",
                              text.color = "black",
                              text.size = 4,
                              textLen = 50){

   myPer <- function(x) round(100*x/sum(x), 2)
   Item <- variable <- value <- label <- Freq <- Var1 <- NULL
   levItem <- colnames(x)

   xLik <- x %>% apply(2, as.character) %>% data.frame(check.names = FALSE) %>%
      mutate_if(is.character, function(x) factor(x, levels = allLevels))
   xLik <- lapply(1:ncol(xLik),
                  function(i) table(xLik[[i]]) %>%
                     as.data.frame %>%
                     dplyr::mutate(Freq = myPer(Freq), Item = colnames(xLik)[i])) %>%
      do.call(what=rbind) %>%
      dplyr::mutate(Item = factor(x = Item, levels = levItem)) %>%
      tidyr::spread(key = Var1, value = Freq, fill = 0)

   results <- reshape2::melt(xLik, id.vars = "Item")
   results$variable <- as.character(results$variable)
   results$label <- paste(format(results$value, digits = 2,
                                 drop0trailing = FALSE), "%", sep = "")
   results$Item <- factor(x = results$Item, levels = levItem %>% rev)

   p <- ggplot(results,
               aes(x = Item, y = variable, fill = value,
                   label = label)) +
      scale_y_discrete(limits = names(xLik)[2:ncol(xLik)]) +
      geom_tile() +
      geom_text(size = text.size, colour = text.color) +
      coord_flip() +
      scale_fill_gradient2("Percent",
                           low = "white",
                           mid = low.color,
                           high = high.color,
                           limits = c(0, 100)) +
      xlab("") +
      ylab("") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) +
      scale_x_discrete(breaks = xLik$Item,
                       labels = base::substring(xLik$Item,
                                                first = 1,
                                                last = textLen))
   return(p)
}

#' Sample-to-sample heatmap
#'
#' sample-to-sample heatmap clustering samples according to a given categorical distance
#' Exploratory tool that helps to visualize/cluster blocks of observations across
#' columns ordered according to given categorical distance. The final output is
#' a clustered distance matrix.
#' This plot is aimed to guide the `DiscreteClusGap` user to give an idea which
#' type of categorical distance would accommodate better to the inputted data.
#' `sample2sampleHeat` is based on the `pheatmap` function from the `pheatmap`
#' R package. Thus, any parameter found in pheatmap can be specified to `sample2sampleHeat`.
#'
#' @importFrom magrittr %>%
#' @import pheatmap
#' @param x matrix object or data.frame
#' @param distName Name of categorical distance to apply.
#' @param ... other valid arguments in pheatmap function
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param clustering_method string; clustering method used by pheatmap
#' @param border_color string; color cell borders. By default, border_color = NA, where no border colors are shown.
#' @return clustered heatmap
#' @export
distanceHeat <- function(x,
                        distName,
                        clustering_method = 'complete',
                        border_color = NA,
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
                      border_color = border_color,
                      treeheight_row=0,
                      ...)

}

#' Discrete Data Heatmap
#'
#' Heatmap assuming a given a distance function and a known number of clusters.
#' Function to display a categorical data matrix given a user defined number of
#' clusters `nCl`, a categorical distance `distName` and a predefined clustering
#' method `FUNcluster`.
#' The output displays a heatmap separating and color-labelling resulting
#' clusters vertically in the rows and allowing unsupervised clustering on
#' questions in the columns. Each cell is colored according to the categorical
#' values provided or found in the data.
#' The clustergram is based on the `pheatmap` function from the pheatmap R package.
#' Thus, any parameter found in pheatmap can be specified to `clusGapDiscrHeat`.
#' This function can be used to examine number of clusters before running
#' `clusGapDiscrHeat` but also after the number of clusters is determined.
#'
#' @importFrom magrittr %>%
#' @import pheatmap
#' @import Polychrome
#' @import RColorBrewer
#' @param x matrix object or data.frame
#' @param nCl number of clusters to plot; if `nCl` is a permutation vector of the first lN integers will rearrange clusters according to the original given ordering.
#' @param distName Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param catVals character string vector with (ordered) categorical values
#' @param clusterFUN Character string with one of the available clustering implementations.
#' Available options are: 'pam' (default) from `cluster::pam`, 'diana' from `cluster::diana`, 'fanny' from `cluster::fanny`.
#' 'agnes-\{average, single, complete, ward, weighted\}' from `cluster::agnes`,
#' 'hclust-\{ward.D, ward.D2, single, complete, average, mcquitty, median, centroid\}' from `stats::hclust`,
#' 'kmodes' from `klar::kmodes` (`weighted = FALSE` and `fast= TRUE`).
#' @param out Specifies the desired output between "heatmap" (default; produce a heatmap), "clusters" (return a `data.frame` with clustering assignments) or "clustersReord" (return a `data.frame` with reorganized clusters)
#' @param seed Seed number.
#' @param clusterNames Either `null` or 'renumber'. When `nCl` is a numerical vector, the cluster ordering is rearranged. `NULL` leaves cluster names as their original cluster assignment. 'renumber' respects the rearrangements but relabels the cluster numbers from top to bottom in ascending order.
#' @param prefObs character string vector of length 1 with a prefix for the observations, in case they come unlabelled or the user wants to anomymize sample IDs.
#' @param rowNames character vector with names of rows according to `x`. By default, `rownames(x)` will be printed in the plot. `rowNames=NULL` prevents from showing names. `prefObs` option takes precedence if is different to `NULL`.
#' @param filename character string with name of file output
#' @param outDir character string with the directory path to save output file
#' @param height numeric height of output plot in inches
#' @param width numeric width of output plot in inches
#' @return png file or ComplexHeatmap object
#' @export
ResHeatmap <- function(x,
                       nCl,
                       distName,
                       catVals,
                       clusterFUN,
                       out = 'heatmap',
                       seed=NULL,
                       clusterNames = NULL,
                       prefObs = NULL,
                       rowNames = rownames(x),
                       filename = NULL,
                       outDir = NULL,
                       height = 10, width = 6){
   ## To Do: specify option where the input data is plotted without any cluster rearrangement
   ## To Do: Specify arbitrary cluster labels.

   Clust <- NULL

   if(length(nCl) > 1){
      clOrd <- nCl
      nCl <- length(nCl)
   }else{
      clOrd <- NULL
   }

   if(!is.matrix(x))
      x <- as.matrix(x)

   if(!is.null(seed))
      set.seed(seed)

   if(nCl > 0){
      ## Requesting cluster assignment
      FUNcluster <- clusterFunSel(clustFun = clusterFUN)


      if(!grepl(pattern = '^kmodes.*', x = clusterFUN)){
         AssignedCls <- FUNcluster(x = distancematrix(x, d = distName),
                                   k = nCl)$cluster
      }else{

         PWdistFun <- function(x, y){
            x <- unlist(x)
            y <- unlist(y)
            distancematrix(rbind(x, y), d = distName)[1]
         }

         AssignedCls <- FUNcluster(x = x,
                                   k = nCl,
                                   distFun = PWdistFun)$cluster
      }

   }else{
      AssignedCls <- rep(0, nrow(x))
   }

   if(!is.null(prefObs)){
      row.names(data) <- paste0(prefObs, 1:nrow(x))
      rownames(x) <- paste0(prefObs, 1:nrow(x))
      showRownames <- TRUE
   }else{
      if(is.null(rowNames)){
         showRownames <- FALSE
      }else{
         showRownames <- TRUE
         rownames(x) <- rowNames
      }
   }
   ## Notice that data is reorganized according to clusters!
   if(is.null(rownames(x)))
      row.names(x) <- 1:nrow(x)

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
                                       show_row_names = showRownames,
                                       # row_labels = rowNames,
                                       show_heatmap_legend=c(TRUE) ,
                                       row_split = rowSplit,
                                       row_title = rowTitle)
      ComplexHeatmap::draw(hmObj,
                           column_title = paste0("DiscreteClusGap Assignments\n",
                                                 "Dist: ", distName, " - Alg: ", clusterFUN),
                           column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                           show_annotation_legend=FALSE)
      grDevices::dev.off()

   }else{
      hmObj <-  ComplexHeatmap::Heatmap(mat = data[, -c(1:2)],
                                        col = myCols,
                                        name = ' ',
                                        show_row_names = showRownames,
                                        # row_labels = rowNames,
                                        show_heatmap_legend=c(TRUE) ,
                                        ## row_title = ' ',
                                        row_split = rowSplit,
                                        row_title = rowTitle,
                                        heatmap_width = unit(width, "in"),
                                        heatmap_height = unit(height, "in"))
      ComplexHeatmap::draw(hmObj,
                           column_title = paste0("DiscreteClusGap Assignments\n", distName),
                           column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                           show_annotation_legend=FALSE)
   }
}

#' MDS Plots for Categorical Data
#'
#' Function to visualize the distribution and spread of categorical data in two dimensions
#' starting from the distance matrix.
#' A cluster assignment vector needs to be provided to show color coded points.
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import ggrepel
#' @param x dist object
#' @param cl character Vector with clustering assignments. Important: the vector order should match
#' number (and if possible names) as the observed data matrix used to generate the distance object `x`.
#' @param type character String specifying the class of MDS: either classic `'MDS'` (default) or
#' Non-metrical Dimensional Scaling `'NMDS'`. The core `MDS` function used is `stats::cmdscale(x, k = 2)` and
#' `vegan::metaMDS(comm = x, distance = "none", k = 2, trymax = 100, autotransform = FALSE)` for `NMDS`,
#' which uses `vegan::monoMDS` by default.
#' @param dotSize numerical String specifying the point size. `dotSize = 3` by default.
#' @param cols character Vector of colors to use with clustering labels. Cluster names should match.
#' to the ones provided in `cl`. The by default (`cols = NULL`) de functions produces highly contrasting colors.
#' @param LabTitle character String for the plot's title. By default `LabTitle = type`.
#' @param filename character string with name of file output. `filename = paste0(type, 'plot.png')` by default.
#' @param outDir character string with the directory path to save output file. `outDir = './'` by default.
#' @param addRowNames logical Single value indicating whether to place observation names next to points.
#' The labels used are the names found in `cl`. If `names(cl) == NULL`, the samples
#' will be labelled in number of appearance. ggrepel package is used to locate the labels.
#' @param labSize numeric Single value indicating the size of labels. `labSize = 3` by default.
#' Non-functional if `addRowNames = FALSE`.
#' @param out character String specify output to obtain: `'plot'` (default), `'data.frame'` or `'object'`. Either
#' `plot` return the resulting ggplot object or
#' `data.frame` for a data.frame with coordinates and corresponding clustering assignments or
#' `object` to obtain raw objects from the chosen `type` above.
#' @return png file and either a data.frame with coordinates and labels or either a list with
#' related MDS results or a `metaMDS` object, depending on the `out` option.
#' @export
plotMDS2 <- function(x,
                     cl,
                     type = 'MDS',
                     cols = NULL,
                     dotSize = 3,
                     LabTitle = type,
                     outDir = './',
                     filename = paste0(type, 'plot.png'),
                     addRowNames = FALSE,
                     labSize = 3,
                     out = 'plot'){

   if(type == 'MDS'){
      mdsRaw <- stats::cmdscale(x, k = 2)
      mdsIn <- mdsRaw %>%
         data.frame(row.names = rownames(x))
   }else if(type == 'NMDS'){
      mdsRaw <- vegan::metaMDS(comm = x,
                               distance = "none",
                               k = 2,
                               trymax = 100,
                               autotransform = FALSE)

      mdsIn <- mdsRaw$points %>%
         data.frame(row.names = rownames(x))
   }
   colnames(mdsIn)  <- c('Dim1', 'Dim2')
   mdsIn$cl <- cl

   if(is.null(cols)){
      myCols <- c('orange',  'gray80', 'slateblue', 'firebrick', 'cyan',
                  'yellow', 'linen', 'gray60', 'tomato', 'navy', 'lightgreen', 'gray40',
                  'blue', 'darkolivegreen', 'gray20')

      Labs <- if(is.factor(cl)) levels(cl) else unique(cl)
      cols <- stats::setNames(object = myCols, nm = Labs)
   }

   ggOut <- ggplot(data = mdsIn,
                   mapping = aes(x = Dim1,
                                 y = Dim2,
                                 fill = cl)) +
      geom_point(shape = 21,
                 size = dotSize,
                 alpha = 0.7) +
      scale_fill_manual(values = cols) +
      theme_bw()+
      labs(title = LabTitle)+
      theme(legend.position = 'bottom')

   if(addRowNames){
      if(!is.null(names(cl)))
         mdsIn$rowNames <- names(cl)
      else
         mdsIn$rowNames <- rownames(mdsIn)

      ggOut <- ggOut +
         geom_text_repel(data = mdsIn,
                         mapping = aes(x = Dim1,
                                       y = Dim2,
                                       label = rowNames),
                         segment.alpha = 0.3,
                         segment.size = 0.3,
                         size = labSize,
                         max.overlaps = Inf)
   }

   ggsave(paste0(outDir, '/', filename),
          ggOut,
          width = 7, height = 7.3)

   if(out == 'data.frame')
      return(mdsIn)
   else if(out == 'object')
      return(mdsRaw)
   else if(out == 'plot')
      ggOut

}
