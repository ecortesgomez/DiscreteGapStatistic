#' Discrete application of clusGap
#' Based on the implementation of the function found in the `cluster` R package
#' @param x categorical/number matrix
#' @param FUNcluster a function that accepts as first argument a matrix like `x`; second argument specifies number of `k` (k=>2) clusters
#' This functions returns a list with a component named `cluster`, a vector of length `n=nrow(x)` of integers from `1:k` indicating observation cluster assignment.
#' @param K.max Maximum number of clusters `k` to consider
#' @param vals String or numerical vector. A vector with all categories to consider when bootstrapping the null distribution sample (FR: Full Range option).
#' By DEFAULT vals=NULL, meaning unique range of categories found in the data will be used when drawing the null (DR: Data Range option).
#' @param integer or logical, determining if “progress” output should be printed. The default prints one bit per bootstrap sample.
#' @param distName Name of categorical distance to apply
#' @param B Number of bootstrap samples.
#' @param verbose integer or logical determining whether progress output should printed while running. By DEFAULT one bit is printed per bootstrap sample.
#' @param ... optionally further arguments for `FUNcluster()`
#'
#' @return a matrix with K.max rows and 4 columns, named "logW", "E.logW", "gap", and "SE.sim",
#' where gap = E.logW - logW, and SE.sim correspond to the standard error of `gap`.
#' @export
clusGapUDiscr <- function (x,
                           FUNcluster,
                           K.max,
                           B = 100,
                           vals=NULL,
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
                    FUNcluster(X, kk)$cluster
                }else{
                    rep.int(1L, nrow(X))
                }
        ## ---------- =  =       -------- kmeans() has 'cluster'; pam() 'clustering'
        ## clus is a vector with cluster assignments per observation

        ## split(ii, clus) %>% sapply(length) %>% paste0(collapse = '-') %>% message
        ## This exception selects only clusters with more than 1 observation!!!
        mySplit <- split(ii, clus)
        mySplit <- mySplit[sapply(mySplit, length) != 1]

        ## 0.5* sum(vapply(split(ii, clus),
        0.5* sum(vapply(mySplit,
                        function(I) {
                            xs <- X[I,, drop=FALSE]
                            sum(distancematrix(xs, d=distName) / nrow(xs),
                                na.rm=TRUE)
                        }, 0.))
    }

    logW <- E.logW <- SE.sim <- numeric(K.max)

    if(verbose) cat("Clustering k = 1,2,..., K.max (= ", K.max,"): .. ", sep='')

    for(k in 1:K.max)
        logW[k] <- log(W.k(x, k))
    if(verbose) cat("done\n")

    ## Column ranges
    rng.x1 <- apply(x, 2L, range)

    ## 2. Bootstrap Piece
    logWks <- matrix(0, B, K.max)
    if(verbose)
        cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n",
            sep="")

    for (b in 1:B) {
        ## Generate "H0"-data as "parametric bootstrap sample" :

        if(is.null(vals)){
            ## Here we're doing the discrete uniform resampling across column ranges.
            z <- apply(rng.x1, 2,
                       function(M, nn) sample(size = nn,
                                              x = seq(M[1], M[2]),
                                              replace = TRUE),
                       nn=n)
        }else{
            ## The null samples over all possible discrete values.
            z <- apply(rng.x1, 2,
                       function(M, nn) sample(size = nn,
                                              x = vals, ## Discrete values to use
                                              replace = TRUE),
                       nn=n)
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

