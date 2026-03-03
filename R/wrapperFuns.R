# Internal helper: close any graphics devices opened after `dev_before`
.dgs_snapshot_devices <- function() {
   dl <- grDevices::dev.list()
   if (is.null(dl)) integer(0) else as.integer(dl)
}

.dgs_close_new_devices <- function(dev_before) {
   if (is.null(dev_before)) dev_before <- integer(0)

   repeat {
      dl <- grDevices::dev.list()
      if (is.null(dl)) break

      dev_after <- as.integer(dl)
      new_devs <- setdiff(dev_after, dev_before)
      if (length(new_devs) == 0L) break

      # Close one at a time to be robust to device renumbering
      d <- max(new_devs)
      try(grDevices::dev.off(which = d), silent = TRUE)
   }

   invisible(NULL)
}

#' Retrieve cluster assignments from a DiscreteGapStatistic heatmap run
#'
#' Convenience wrapper around \code{\link{ResHeatmap}} that extracts cluster labels
#' (optionally re-ordered to match the input \code{data}) and returns either a
#' data.frame or a named character vector.
#'
#' @param data A matrix or data.frame of categorical observations (rows are observations,
#'   columns are variables). Row names should be present if \code{ordering = "data"}.
#' @param catVals Vector of possible category values used by \code{\link{ResHeatmap}}.
#' @param clusterFUN Character name of a clustering algorithm supported by the workflow
#'   (e.g., \code{"pam"}).
#' @param distName Character name of the distance to use.
#' @param nCl Integer number of clusters.
#' @param outFormat Output format: \code{"data.frame"} (default) or \code{"vector"}.
#' @param clusterNames Optional cluster names passed to \code{\link{ResHeatmap}}.
#' @param ordering If \code{"data"}, re-order to match \code{rownames(data)}; if
#'   \code{"heatmap"}, keep \code{\link{ResHeatmap}} ordering.
#'
#' @return A data.frame (default) or a named character vector of cluster labels.
#'
#' @export
#'
#' @examples
#' # x <- matrix(sample(letters[1:3], 120, TRUE), nrow = 30)
#' # rownames(x) <- paste0("s", seq_len(nrow(x)))
#' # RetrClustAssign(x, catVals = letters[1:3], distName = "hamming", nCl = 3)
#' # RetrClustAssign(x, catVals = letters[1:3], distName = "hamming", nCl = 3,
#' #                 outFormat = "vector")
RetrClustAssign <- function(data,
                            catVals,
                            clusterFUN = "pam",
                            distName,
                            nCl,
                            outFormat = c("data.frame", "vector"),
                            clusterNames = NULL,
                            ordering = c("data", "heatmap")) {
   outFormat <- match.arg(outFormat)
   ordering <- match.arg(ordering)

   outClust <- ResHeatmap(
      x = data,
      catVals = catVals,
      clusterFUN = clusterFUN,
      distName = distName,
      nCl = nCl,
      clusterNames = clusterNames,
      out = "clustersReord"
   )

   outClust$Clust <- paste0("cl", outClust$Clust)

   if (ordering == "data") {
      if (is.null(rownames(data))) {
         stop("`ordering = \"data\"` requires `rownames(data)` to be non-NULL.", call. = FALSE)
      }
      idx <- match(rownames(data), outClust$rowNames)
      if (anyNA(idx)) {
         stop("Some `rownames(data)` were not found in `outClust$rowNames`.", call. = FALSE)
      }
      outClust <- outClust[idx, , drop = FALSE]
   }

   if (outFormat == "vector") {
      outClust <- stats::setNames(outClust$Clust, nm = outClust$rowNames)
   }

   outClust
}

#' Run the Discrete Gap Statistic workflow and save plots
#'
#' High-level wrapper that (i) saves a distance heatmap, (ii) runs
#' \code{\link{clusGapDiscr}}, (iii) saves a Gap Statistic plot (with the chosen
#' \code{K} marked), (iv) saves the resulting categorical heatmap, and (v) attempts
#' to save an MDS plot colored by the chosen clusters.
#'
#' This function includes a device "safety rail" that closes any graphics devices
#' opened during its execution (useful if downstream plotting code forgets to call
#' \code{dev.off()}).
#'
#' @param x A matrix of categorical observations. Must be a matrix (enforced).
#' @param catVals Vector of possible category values for the variables in \code{x}.
#' @param dataClass Data class indicator passed to \code{\link{clusGapDiscr}}.
#' @param clusterFUN Character name of a clustering algorithm (default \code{"pam"}).
#' @param B Integer number of bootstrap samples for the gap statistic.
#' @param K.max Maximum number of clusters considered.
#' @param value.range Character string passed to \code{\link{clusGapDiscr}}.
#' @param distName Character name of the distance metric (default \code{"hamming"}).
#' @param useLog Logical; passed to \code{\link{clusGapDiscr}}.
#' @param title Optional title prefix used in plot titles and output filenames.
#' @param outDir Output directory where PNG files will be written.
#'
#' @return The object returned by \code{\link{clusGapDiscr}}.
#'
#' @export
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot abline
DGSrun <- function(x,
                   catVals,
                   dataClass,
                   clusterFUN = "pam",
                   B = 100,
                   K.max = 7,
                   value.range = "DS",
                   distName = "hamming",
                   useLog = TRUE,
                   title = NULL,
                   outDir = "./") {
   stopifnot(is.matrix(x))

   if (!dir.exists(outDir)) {
      dir.create(outDir, recursive = TRUE)
   }

   # Safety rail: close any devices opened during this function
   dev_start <- .dgs_snapshot_devices()
   on.exit(.dgs_close_new_devices(dev_start), add = TRUE)

   logDGS <- ifelse(useLog, "DiscreteGapStatistic", "DiscreteGapStatistic*")
   logDGSsh <- ifelse(useLog, "DGS", "DGS-star")

   fileName <- paste0(
      logDGSsh, "_", clusterFUN, "_",
      gsub("_", "-", distName), "_", value.range
   )
   titleOut <- paste0(
      logDGS, "\n", clusterFUN, " - ",
      gsub("_", "-", distName), " - ", value.range
   )

   if (!is.null(title)) {
      fileName <- paste0(title, "_", fileName)
      titleOut <- paste0(title, "\n", titleOut)
   }

   message("Plotting the distance heatmap")
   dev0 <- .dgs_snapshot_devices()
   distanceHeat(
      x = x,
      clusterFUN = clusterFUN,
      distName = distName,
      main = titleOut,
      fontsize = 6.5,
      cluster.cols = TRUE,
      title = titleOut,
      show_rownames = FALSE,
      show_colnames = FALSE,
      filename = file.path(outDir, paste0("distHeat_", fileName, ".png"))
   )
   .dgs_close_new_devices(dev0)

   message("Running DGS")
   curRun <- clusGapDiscr(
      x = x,
      clusterFUN = clusterFUN,
      B = B,
      K.max = K.max,
      value.range = value.range,
      distName = distName,
      dataClass = dataClass,
      useLog = useLog
   )

   k_best <- findK(curRun)

   # Gap Statistic plot (device opened here, so guarantee closure)
   message("Saving Gap Statistic plot")
   dev1 <- .dgs_snapshot_devices()
   grDevices::png(
      filename = file.path(outDir, paste0("GapStatPlot_", fileName, ".png")),
      width = 7, height = 7, res = 500, units = "in"
   )
   on.exit(.dgs_close_new_devices(dev1), add = TRUE)

   graphics::plot(
      curRun,
      main = titleOut,
      cex = 2, cex.lab = 1.2, cex.axis = 1.5, cex.main = 1.5
   )
   graphics::abline(v = k_best, lty = 3, lwd = 2, col = "blue")
   .dgs_close_new_devices(dev1)

   message("Plotting resulting heatmap")
   dev2 <- .dgs_snapshot_devices()
   ResHeatmap(
      x = x,
      distName = distName,
      clusterFUN = clusterFUN,
      catVals = catVals,
      nCl = k_best,
      out = "heatmap",
      rowNames = NULL,
      prefObs = NULL,
      height = 7,
      filename = paste0("ResHeatmap_", fileName),
      outDir = outDir
   )
   .dgs_close_new_devices(dev2)

   ClLabs <- RetrClustAssign(
      data = x,
      catVals = catVals,
      distName = distName,
      nCl = k_best,
      outFormat = "vector"
   )

   message("Plotting MDS")
   dev3 <- .dgs_snapshot_devices()
   try({
      dm <- distancematrix(x, d = distName)
      plotMDS2(
         dm,
         cl = ClLabs,
         type = "MDS",
         filename = paste0("MDSplot_", fileName, ".png"),
         outDir = outDir,
         out = "plot"
      )
   }, silent = TRUE)
   .dgs_close_new_devices(dev3)

   curRun
}
