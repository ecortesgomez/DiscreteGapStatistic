#' Discrete application of clusGap - core function.
#'
#' Based on the implementation of the function found in the `cluster` R package.
#' This function assumes that all attributes have identical categories.
#'
#' @param x A matrix object specifying category attributes in the columns and observations in the rows.
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This function should return a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' Make sure `FUNcluster` and `Input2Alg` agree.
#' @param K.max Integer. Maximum number of clusters `k` to consider
#' @param value.range A length 1 character string, a character string vector or a list of character vector with the length matching the number of columns (nQ) of the array.
#' By DEFAULT value.range = 'DS' (Data Support null model).
#' A vector with all categories (either character for nominal or integer ordinal data) to consider when bootstrapping the null distribution sample (KS: Known Support option).
#' If a list with category vectors is provided, it has to have the same number of columns as the input array. The order of list element corresponds to the array's columns.
#' @param verbose Integer or logical. Determines if “progress” output should be printed. The default prints one bit per bootstrap sample.
#' @param distName String. Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param Input2Alg Specifies the kind of input provided to the algorithm function in `FUNcluster`.
#' For algorithms that only accept a distance matrix use `'distMatr'` option (default).
#' For algorithms that require the dataset and a prespecified distance function (e.g. `stats::dist`) use the `'distFun'` option.
#' This case the distance function is defined internally and determined by parameter `distName`.
#' @param B Number of bootstrap samples. By default `B = nrow(x)`.
#' @param verbose Integer or logical. Determines whether progress output should printed while running. By DEFAULT one bit is printed per bootstrap sample.
#' @param useLog Logical. Use log function after estimating `W.k`. Following the original formulation `useLog=TRUE` by default.
#' @param dataClass character. Either 'nom' for nominal or 'ord' for ordinal.
#' @param offset Numeric. A nonnegative constant added to both observed and bootstrap W.k
#' before taking logarithms when `useLog=TRUE`. `offset = 0` by default.
#' @param ... optionally further arguments for `FUNcluster()`
#'
#' @return A `clusGap` object. Normally its `Tab` component has `K.max` rows. If a
#' non-finite observed or bootstrap dispersion is encountered, the table is truncated to
#' the largest contiguous finite value of k. The object records `K.max.requested`,
#' `K.max.used`, and `truncated`.

clusGapDiscr0 <- function (x,
                          FUNcluster,
                          K.max,
                          B = nrow(x),
                          value.range = "DS",
                          verbose = interactive(),
                          distName = "hamming",
                          useLog = TRUE,
                          Input2Alg = 'distMatr',
                          dataClass = 'nom',
                          offset = 0,
                          ...){

   stopifnot(is.function(FUNcluster),
             length(dim(x)) == 2,
             K.max >= 2,
             (n <- nrow(x)) >= 1, ## Define n !!
             ncol(x) >= 1)

   if(!is.matrix(x))
      x <- as.matrix(x)

   uniLevs <- unique(as.vector(x))
   uniLevs <- sort(uniLevs)

   ## Nominal case
   ## The only way inside of nominal functions is with character vector!
   if (is.numeric(x) & grepl('^nom', dataClass)) {
      message(paste0("x array is numerical and has ", length(uniLevs),
                     " levels."))
      ## stopifnot(length(uniLevs) < 10) ## why?
      message("Array values would be transformed to categorical levels: ")
      myCats <- paste0(paste0(uniLevs, " -> c", uniLevs), collapse = ", ")
      message(myCats)
      x <- apply(x, 1, function(y) paste0('c', y)) %>% t ## Transformed matrix

      message(paste0("Resulting levels: ", paste0(uniLevs, collapse = ", ")))
   }else if(is.numeric(x) & grepl('^ord', dataClass) ){
      ## Ordinal case
      ## Only way to pass ordinal data is with numerical vector.
      message(paste0("Found levels: ", paste0(uniLevs, collapse = ", ")))
   }else if(is.character(x) & grepl('^ord', dataClass) ){
      stopifnot(!is.character(x) & grepl('^ord', dataClass))
   }

   if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
      stop("'B' has to be a positive integer")

   if (length(offset) != 1L || !is.numeric(offset) ||
       !is.finite(offset) || offset < 0) {
      stop("'offset' must be one finite, nonnegative numeric value")
   }

   K.max.requested <- as.integer(K.max)
   cl. <- match.call()

   ii <- seq_len(n)

   if(Input2Alg == 'distMatr'){

      W.k <- function(X, kk) {
         clus <- if (kk > 1) {
            Xdist <- distancematrix(X, d = distName)
            FUNcluster(Xdist, kk)$cluster
         }
         else {
            rep.int(1L, nrow(X))
         }
         mySplit <- split(ii, clus)
         mySplit <- mySplit[sapply(mySplit, length) != 1]
         0.5 * sum(vapply(mySplit, function(I) {
         ## sum(vapply(mySplit, function(I) {
            xs <- X[I, , drop = FALSE]
            sum(distancematrix(xs, d = distName)/nrow(xs), na.rm = TRUE)
         }, 0))
      }

   }else if(Input2Alg == 'distFun'){

      PWdistFun <- function(x, y){
         x <- unlist(x)
         y <- unlist(y)
         distancematrix(rbind(x, y), d = distName)[1]
      }

      myFUNcluster <- function(X, kk){
         ## Condition: FUNcluster is expecting a distFun parameter
         FUNcluster(X,
                    kk,
                    distFun = function(x, y) PWdistFun(x, y) )
      }

      W.k <- function(X, kk){
         clus <- if (kk > 1) {
            ## Xdist <- distancematrix(X, d = distName)
            myFUNcluster(X, kk)$cluster
         }
         else {
            rep.int(1L, nrow(X))
         }
         mySplit <- split(ii, clus) ## Split into a list with potential clusters
         mySplit <- mySplit[sapply(mySplit, length) != 1] ## Removes singletons
         # 0.5 * sum(vapply(mySplit, function(I) {
         sum(vapply(mySplit, function(I) {
            xs <- X[I, , drop = FALSE]
            sum(distancematrix(xs, d = distName)/nrow(xs), na.rm = TRUE)
         }, 0))
      }

   }

   transform_W <- function(w) {
      if (useLog) {
         log(w + offset)
      } else {
         w
      }
   }

   logW <- rep(NA_real_, K.max.requested)
   K.max.used <- K.max.requested
   limiting.event <- NULL

   if (verbose)
      cat("Clustering k = 1,2,..., K.max (= ", K.max.requested,
          "): .. ", sep = "")

   ## Observed-data dispersions are evaluated first. Once a non-finite value
   ## occurs, larger k values are not used because the valid range must be a
   ## contiguous sequence 1, ..., K.max.used.
   for (k in seq_len(K.max.requested)) {
      wk <- W.k(x, k)
      logW[k] <- transform_W(wk)

      if (!is.finite(logW[k])) {
         K.max.used <- k - 1L
         limiting.event <- list(
            source = "observed data",
            k = k,
            bootstrap = NA_integer_,
            W.k = wk
         )
         break
      }
   }

   if (K.max.used < 1L) {
      stop(
         "Non-finite ", if (useLog) "logW" else "W.k",
         " occurred at k = 1 for the observed data. No valid gap-statistic ",
         "range is available. Consider a positive offset or useLog = FALSE."
      )
   }

   if (verbose)
      cat("done\n")

   if(grepl('^nom', dataClass)){

      ## Nominal case
      if(is.numeric(value.range)) {
         stop('value.range cannot be numerical and will be transformed to characters.')
      }
      else if (is.character(value.range) & !is.list(value.range) &
             value.range[1] == "DS" & length(value.range) == 1) {
         ## Data-Support null option

         rng.x1 <- lapply(1:ncol(x), function(i) unique(x[, i]))
      }
      else if (all(is.character(value.range)) & !is.list(value.range) & (length(value.range) > 1)) {
         ## KS: User-defined vector

         rng.x1 <- rep(x = list(value.range), times = ncol(x))
      }
      else if (!is.character(value.range) & is.list(value.range)) {
         ## When user defines the value range for each question

         message('Known Support has variability between attributes.')
         stopifnot(ncol(x) == length(value.range))
         rng.x1 <- value.range

      }
   }else if(grepl('^ord', dataClass)){
      ## Ordinal case

      if(!is.list(value.range) & value.range[1] == 'DS' & length(value.range) == 1){
         rng.x1 <- lapply(1:ncol(x), function(i) unique(x[, i]))
      }
      else if(is.integer(value.range) & is.vector(value.range)){
         ## KS case

         rng.x1 <- rep(x = list(value.range), times = ncol(x))
        }
      else if (!is.integer(value.range) & is.list(value.range)) {
         ## When user defines the value range for each question
         message('Known Support has variability between attributes.')
         stopifnot(ncol(x) == length(value.range))
         rng.x1 <- value.range
      }
   }

   logWks <- matrix(
      NA_real_,
      nrow = B,
      ncol = K.max.used
   )

   if (verbose)
      cat("Bootstrapping, b = 1,2,..., B (= ", B,
          ")  [one \".\" per sample]:\n", sep = "")

   for (b in seq_len(B)) {
      z <- lapply(
         rng.x1,
         function(M, nn) sample(size = nn, x = M, replace = TRUE),
         nn = n
      ) %>%
         do.call(what = cbind)

      ## K.max.used can decrease during bootstrapping. Each bootstrap sample is
      ## evaluated only over the currently valid contiguous range.
      for (k in seq_len(K.max.used)) {
         wk <- W.k(z, k)
         logWks[b, k] <- transform_W(wk)

         if (!is.finite(logWks[b, k])) {
            new.K.max <- k - 1L

            if (new.K.max < 1L) {
               stop(
                  "Non-finite ", if (useLog) "bootstrap logW" else "bootstrap W.k",
                  " occurred at bootstrap sample ", b, ", k = 1. No valid ",
                  "gap-statistic range is available. Consider a positive offset ",
                  "or useLog = FALSE."
               )
            }

            K.max.used <- new.K.max
            limiting.event <- list(
               source = "bootstrap reference data",
               k = k,
               bootstrap = b,
               W.k = wk
            )
            break
         }
      }

      if (verbose)
         cat(".", if (b %% 50 == 0) paste(b, "\n"))
   }

   if (verbose && (B %% 50 != 0))
      cat("", B, "\n")

   ## Retain only the common finite prefix across the observed data and every
   ## bootstrap reference sample.
   keep.k <- seq_len(K.max.used)
   logW <- logW[keep.k]
   logWks <- logWks[, keep.k, drop = FALSE]

   E.logW <- colMeans(logWks)
   SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, stats::var))

   truncated <- K.max.used < K.max.requested

   if (truncated) {
      event.location <- if (is.na(limiting.event$bootstrap)) {
         paste0("observed data at k = ", limiting.event$k)
      } else {
         paste0(
            "bootstrap sample ", limiting.event$bootstrap,
            " at k = ", limiting.event$k
         )
      }

      warning(
         "Non-finite ", if (useLog) "logW" else "W.k",
         " encountered for ", event.location, ". Requested K.max = ",
         K.max.requested, "; results were returned only through K.max = ",
         K.max.used, ".",
         call. = FALSE
      )
   }

   return(structure(
      class = "clusGap",
      list(
         Tab = cbind(
            logW,
            E.logW,
            gap = E.logW - logW,
            SE.sim
         ),
         call = cl.,
         n = n,
         B = B,
         FUNcluster = FUNcluster,
         useLog = useLog,
         distName = distName,
         value.range = value.range,
         offset = offset,
         K.max.requested = K.max.requested,
         K.max.used = K.max.used,
         truncated = truncated,
         limiting.event = limiting.event
      )
   ))
}

#' Discrete application of clusGap
#'
#' Based on the implementation of the function found in the `cluster` R package.
#'
#' @param x A matrix object specifying category attributes in the columns and observations in the rows.
#' @param clusterFUN Character string with one of the available clustering implementations.
#' Available options are: 'pam' (default) from `cluster::pam`, 'diana' from `cluster::diana`, 'fanny' from `cluster::fanny`,
#' 'agnes-\{average, single, complete, ward, weighted\}' from `cluster::fanny`,
#' 'hclust-\{ward.D, ward.D2, single, complete, average, mcquitty, median, centroid\}' from `stats::hclust`,
#' 'kmodes' from `klar::kmodes` (`iter.max = 10`, `weighted = FALSE` and `fast= TRUE`).
#' 'kmodes-N' enables to run the `kmodes` algorithm with a given number N of iterations where `iter.max = N`.
#' @param K.max Integer. Maximum number of clusters `k` to consider
#' @param value.range A length 1 character string, a character string vector or a list of character vector with the length matching the number of columns (nQ) of the array.
#' By DEFAULT value.range = 'DS' (Data Support null model).
#' A vector with all categories (either character for nominal or integer ordinal data) to consider when bootstrapping the null distribution sample (KS: Known Support option).
#' If a list with category vectors is provided, it has to have the same number of columns as the input array. The order of list element corresponds to the array's columns.
#' @param verbose Integer or logical. Determines if “progress” output should be printed. The default prints one bit per bootstrap sample.
#' @param distName String. Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'.
#' @param B Number of bootstrap samples. By default `B = nrow(x)`.
#' @param verbose Integer or logical. Determines whether progress output should printed while running. By DEFAULT one bit is printed per bootstrap sample.
#' @param useLog Logical. Use log function after estimating `W.k`. Following the original formulation `useLog=TRUE` by default.
#' @param dataClass character. Either 'nom' for nominal or 'ord' for ordinal.
#' @param offset Numeric. A nonnegative constant added to both observed and bootstrap W.k
#' before taking logarithms when `useLog=TRUE`. `offset = 0` by default.
#' @param ... optionally further arguments for `FUNcluster()`
#'
#' @return A `clusGap` object. Normally its `Tab` component has `K.max` rows. If a
#' non-finite observed or bootstrap dispersion is encountered, the table is truncated to
#' the largest contiguous finite value of k. The object records `K.max.requested`,
#' `K.max.used`, and `truncated`.
#' @export

clusGapDiscr <- function(x,
                         clusterFUN,
                         K.max,
                         B = nrow(x),
                         value.range = "DS",
                         verbose = interactive(),
                         distName = "hamming",
                         useLog = TRUE,
                         dataClass = 'nom',
                         offset = 0,
                         ...){

   if(clusterFUN == 'pam'){
      clusGapDiscr0(x = x,
                    FUNcluster = cluster::pam,
                    K.max = K.max,
                    B = B,
                    value.range = value.range,
                    verbose = verbose,
                    distName = distName,
                    useLog = useLog,
                    dataClass = dataClass,
                    offset = offset,
                     ...)

   }else if(clusterFUN == 'fanny'){
      clusGapDiscr0(x = x,
                    FUNcluster = cluster::fanny,
                    K.max = K.max,
                    B = B,
                    value.range = value.range,
                    verbose = verbose,
                    distName = distName,
                    useLog = useLog,
                    dataClass = dataClass,
                    offset = offset,
                    ...)

   }else if(clusterFUN == 'diana'){
      dianaK <- function(d, k){
         out <- cluster::diana(x = d) %>%
                 stats::cutree(k = k)
         list(cluster = out)
      }

      clusGapDiscr0(x = x,
                    FUNcluster = dianaK,
                    K.max = K.max,
                    B = B,
                    value.range = value.range,
                    verbose = verbose,
                    distName = distName,
                    useLog = useLog,
                    dataClass = dataClass,
                    offset = offset,
                    ...)

   }else if(grepl(pattern = '^agnes-.+', x = clusterFUN)){

      cMeth <- sub('agnes-', '', clusterFUN)
      agnesK <- function(d = d, k = k){
         out <- cluster::agnes(x = d,
                               diss = TRUE,
                               method = cMeth) %>%
            stats::cutree(k = k)
         list(cluster = out)
      }

      clusGapDiscr0(x = x,
                    FUNcluster = agnesK,
                    K.max = K.max,
                    B = B,
                    value.range = value.range,
                    verbose = verbose,
                    distName = distName,
                    useLog = useLog,
                    dataClass = dataClass,
                    offset = offset,
                    ...)

   }else if(grepl(pattern = '^hclust-.+', x = clusterFUN)){

      cMeth <- sub('hclust-', '', clusterFUN)
      hClustK <- function(d = d, k = k){
            out <- stats::hclust(d = d, method = cMeth) %>%
               stats::cutree(k = k)
            list(cluster = out)
         }

      clusGapDiscr0(x = x,
                    FUNcluster = hClustK,
                    K.max = K.max,
                    B = B,
                    value.range = value.range,
                    verbose = verbose,
                    distName = distName,
                    useLog = useLog,
                    dataClass = dataClass,
                    offset = offset,
                    ...)

   }else if(grepl(pattern = '^kmodes.*', x = clusterFUN)){

      cIter <- ifelse(grepl(pattern = 'kmodes-.+',  x = clusterFUN),
                      as.numeric(sub('kmodes-', '', clusterFUN)),
                      10)

      kmodesD_ <- function(x, k, distFun){
         kmodesD(data = x, modes = k, distFun = distFun, iter.max = cIter)
      }

      clusGapDiscr0(x = x,
                    FUNcluster = kmodesD_,
                    K.max = K.max,
                    B = B,
                    value.range = value.range,
                    verbose = verbose,
                    distName = distName,
                    useLog = useLog,
                    Input2Alg = 'distFun',
                    dataClass = dataClass,
                    offset = offset,
                    ...)
   }else
      message('Clustering algorithm not available.')

}

#' Criteria to determine number of clusters k
#'
#' Similar to `maxSE` function found the `cluster` package.
#'
#' @import dplyr
#' @param cG_obj Output object obtained from `clusGapDiscr`
#' @param meth Method to use to determine optimal k number of clusters.
#' @param SE.fact Standard Error Factor generalizing the 1-SE rule.
#' @return A numerical value from 1 to K.max, contained in the input `cG_obj` object.
#' @export

findK <- function (cG_obj, meth = "Tibs2001SEmax", SE.fact = 1){
   logW <- gap <- SE.sim <- NULL

   outDF <- data.frame(cG_obj['Tab'])
   outDF[, 'nClus'] <- 1:nrow(outDF)
   colnames(outDF) <- c('logW', 'E.logW', 'gap', 'SE.sim', 'nClus')

   myTab <- outDF %>%
      subset(is.finite(logW) & is.finite(gap) & is.finite(SE.sim))

   if (nrow(myTab) < 1L) {
      stop(
         "No finite logW, gap, and SE.sim values are available for cluster ",
         "selection."
      )
   }

   selK <- with(myTab,
                cluster::maxSE(f = gap,
                               SE.f = SE.sim,
                               method = meth,
                               SE.factor = SE.fact) )
   myTab[selK, 'nClus']
}
