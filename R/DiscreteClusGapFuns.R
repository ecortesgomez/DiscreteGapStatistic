#' Discrete application of clusGap
#' Based on the implementation of the function found in the `cluster` R package
#' @param x categorical/number matrix
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This functions returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' @param K.max Maximum number of clusters `k` to consider
#' @param value.range String, character vector or a list of character vector with the length matching the number of columns (nQ) of the array.
#' A vector with all categories to consider when bootstrapping the null distribution sample (FR: Full Range option).
#' By DEFAULT vals=NULL, meaning unique range of categories found in the data will be used when drawing the null (DR: Data Range option).
#' If a character vector of categories is provided, these values would be used for the null distribution drawing across the array.
#' If a list with category character vectors is provided, it has to have the same number of columns as the input array. The order of list element corresponds to the array's columns.
#' @param integer or logical, determining if “progress” output should be printed. The default prints one bit per bootstrap sample.
#' @param distName Name of categorical distance to apply.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming', 'hellinger',
#' @param B Number of bootstrap samples. By default B = nrow(x).
#' @param verbose integer or logical determining whether progress output should printed while running. By DEFAULT one bit is printed per bootstrap sample.
#' @param ... optionally further arguments for `FUNcluster()`
#'
#' @return a matrix with K.max rows and 4 columns, named "logW", "E.logW", "gap", and "SE.sim",
#' where gap = E.logW - logW, and SE.sim correspond to the standard error of `gap`.
#' @export
clusGapUDiscr <- function (x,
                           FUNcluster,
                           K.max,
                           B = nrow(x),
                           value.range = 'DR',
                           verbose = interactive(),
                           distName = 'hamming', ...
                           ## seed = 17,
                           ){

    ## 0. Preliminary checkups
    ## #########################
    ## set.seed(seed)

    stopifnot(is.function(FUNcluster),
              length(dim(x)) == 2,
              K.max >= 2,
              (n <- nrow(x)) >= 1,
              ncol(x) >= 1)

    ## Reformatting if the array is found to be numerical
   uniLevs <- unique(as.vector(x))
   if(is.numeric(x)){
       message(paste0('x array is numerical and has ', length(uniLevs), ' levels.'))
       stopifnot(length(uniLevs) < 10)
       message('Array values would be transformed to categorical levels: ')
       myCats <- paste0(paste0(uniLevs, ' -> c', uniLevs), collapse = ', ')
       message(myCats)
    }else{
       message(paste0('Found levels: ',
                      paste0(uniLevs,
                             collapse = ', ')) )
    }

    if(B != (B. <- as.integer(B)) || (B <- B.) <= 0)
         stop("'B' has to be a positive integer")
    cl. <- match.call()

    ## 1. Estimate log(W(k)) for the data
    ## ###################################

    if(is.data.frame(x))
        x <- as.matrix(x)
    ii <- seq_len(n)
    ## debugonce(W.k)
    ## W.k(X=x, kk=2)
    W.k <- function(X, kk) {
        ## This is where the distance functions comes into place!

        clus <- if(kk > 1){
                    Xdist <- distancematrix(X, d = distName)
                    FUNcluster(Xdist, kk)$cluster
                }else{
                    rep.int(1L, nrow(X))
                }
        ## ---------- =  =       -------- kmeans() has 'cluster'; pam() 'clustering'
        ## clus is a vector with cluster assignments per observation

        ## split(ii, clus) %>% sapply(length) %>% paste0(collapse = '-') %>% message
        ## This exception selects only clusters with more than 1 observation!!!
        mySplit <- split(ii, clus)
        mySplit <- mySplit[sapply(mySplit, length) != 1]

        0.5* sum(vapply(mySplit,
                        function(I) {
                            xs <- X[I,, drop=FALSE]
                            sum(distancematrix(xs, d = distName) / nrow(xs),
                                na.rm=TRUE)
                        }, 0.))
    }

    logW <- E.logW <- SE.sim <- numeric(K.max)

    if(verbose) cat("Clustering k = 1,2,..., K.max (= ", K.max,"): .. ", sep='')

    for(k in 1:K.max)
        logW[k] <- log(W.k(x, k))
    if(verbose) cat("done\n")

    ## Column ranges
    if(is.character(value.range) & !is.list(value.range) & value.range[1] == 'DR' &
       length(value.range) == 1){
         vals <- NULL
         ## rng.x1 <- apply(x, 2L, unique)
         rng.x1 <- lapply(1:ncol(x), function(i) unique(x[, i]))

    }else if(all(is.character(value.range) ) &
             !is.list(value.range) &
             (length(value.range) > 1) ){
         vals <- value.range

    }else if(!is.character(value.range) & is.list(value.range)){
         stopifnot(ncol(x) == length(value.range))
         vals <- value.range
    }

    ## 2. Bootstrap Piece
    logWks <- matrix(0, B, K.max)
    if(verbose)
        cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n",
            sep="")

    for (b in 1:B) {
        ## Generate "H0"-data as "parametric bootstrap sample" :

        if(is.null(vals)){
            ## DR: Data Range option
            z <- lapply(rng.x1,
                       function(M, nn) sample(size = nn,
                                              x = M,
                                              replace = TRUE),
                       nn=n) %>%
               do.call(what = cbind)

        }else if(!is.list(vals) & length(vals) > 2 ){
            ## Single vector with all levels
            z <- matrix(sample(x = vals, size = nrow(x)*ncol(x), replace=TRUE),
                        nrow = nrow(x))

        }else if(is.list(vals)){
           z <- sapply(1:length(vals),
                      function(i, nn) sample(size = nn,
                                             x = vals[[i]], ## Discrete values to use
                                             replace = TRUE),
                      nn=n)
           ##
        }

        for(k in 1:K.max) {
            logWks[b,k] <- log(W.k(z, k))
        }
        if(verbose) cat(".", if(b %% 50 == 0) paste(b,"\n"))
        ## Show progress every 50 observations
    }

    if(verbose && (B %% 50 != 0)) cat("", B,"\n")

    E.logW <- colMeans(logWks)
    SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, stats::var))

    ## 3. Release Output
    ## ####################

    structure(class = "clusGap",
              list(Tab = cbind(logW,
                               E.logW,
                               gap = E.logW - logW,
                               SE.sim),
                   ## K.max == nrow(T)
                   call = cl.,
                   n = n,
                   B = B,
                   FUNcluster = FUNcluster))
}

#' Criteria to determine k
#' @param cG_obj Output object obtained from `clusGapUDiscr`
#' @param meth Method to use to determine optimal k number of clusters.
#'
findK <- function(cG_obj, meth='Tibs2001SEmax'){
    if(!meth %in% c('minSE', 'minGap', 'maxChange')){
        cluster::maxSE(f = cG_obj$Tab[, "gap"],
              SE.f = cG_obj$Tab[, "SE.sim"],
              method = meth)
    }else if(meth == 'minSE'){
        which.min(cG_obj$Tab[, "SE.sim"])
    }else if(meth == 'minGap'){
        which.min(cG_obj$Tab[, "gap"])
    }else if(meth == 'maxChange'){
        ## This is not right, since it fails to select 1 cluster
        which.max(abs(cG_obj$Tab[, "gap"])) + 1
    }
}

