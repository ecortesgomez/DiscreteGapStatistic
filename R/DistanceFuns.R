#' Calculate categorical distance matrix for discrete data
#'
#' Function invoking discrete distance functions.
#' Available distances: 'bhattacharyya', 'chisquare', 'cramerV', 'hamming' and 'hellinger'
#'
#' @param X Matrix where rows are the observations and columns are discrete features
#' @param d Name of distance. Distances available: bhattacharyya, chisquare, cramerV, hamming and hellinger
#'
#' @return R distance object
#' @export
#'
#' @examples
#' X = rbind(matrix(paste0("a", rpois(7*5, 1)), nrow=5),
#'           matrix(paste0("a", rpois(7*5, 3)), nrow=5))
#' distancematrix(X = X, d = "hellinger")

distancematrix <- function (X, d){
   #####################
   # Nominal distances #
   #####################

    if(d == 'bhattacharyya')
        return(dissbhattacharyya(X))
    if(d == 'chisquare')
        return(disschisquare(X))
    if (d == "cramerV")
         return(disscramerv(X))
    if (d == "hamming") ## General version
        return(disshamming(X))
    if (d == "hellinger")
        return(disshellinger(X))

   ## nomclust functions do not like the matrix class!
   ## Transform to data.frame
   if(grepl('nmcl_.+', d)){
      if(grepl('nmcl_good.*', d)){ ## goodall diss need special treatment
         X <- data.frame(X, check.names=FALSE) %>% apply(2, factor)
      }else{
         X <- as.data.frame(X)
      }
   }

   if (d == "nmcl_anderberg")
      return(nomclust::anderberg(X))
   if (d == "nmcl_burnaby")
      return(nomclust::burnaby(X))
   if (d == "nmcl_eskin")
      return(nomclust::eskin(X))
   if (d == "nmcl_gambaryan")
      return(nomclust::gambaryan(X))
   if (d == "nmcl_goodall1")
      return(nomclust::goodall1(X))
   if (d == "nmcl_goodall2")
      return(nomclust::goodall2(X))
   if (d == "nmcl_goodall3")
      return(nomclust::goodall3(X))
   if (d == "nmcl_goodall4")
      return(nomclust::goodall4(X))
   if (d == "nmcl_iof")
      return(nomclust::iof(X))
   if (d == "nmcl_lin")
      return(nomclust::lin(X))
   if (d == "nmcl_lin1")
      return(nomclust::lin1(X))
   if (d == "nmcl_of")
      return(nomclust::of(X))
   if (d == "nmcl_sm")
      return(nomclust::sm(X))
   if (d == "nmcl_smirnov")
      return(nomclust::smirnov(X))
   if (d == "nmcl_ve")
      return(nomclust::ve(X))
   if (d == "nmcl_vm")
      return(nomclust::vm(X))

   #####################
   # Ordinal distances #
   #####################

   if(d == 'absolute')
      return(dissabs(X))

   if(d == 'ks')
      return(dissks(X))

   if(d == 'podani')
      return(disspodani(X))

   if (d == "spearman")
      return(dissspearman(X))

   if (d == 'tau')
      return(disstau(X))

   if (grepl('wasserstein', d)){
      if(grepl('wass.+KS[0-9]+', d)){
         myM <- sub(pattern = 'wass.+KS', replacement = '', x = d) %>%
            as.numeric
         return(dissWass(X, type = 'KS', m = myM))
      }else if(grepl('wass.+DS$', d)){
         return(dissWass(X, type = 'DS'))
      }else if(grepl('wass.+R$', d)){
         return(dissWass(X, type = 'R'))
      }
   }

   stop("Distance metric ", d, " not available")
}

#' Nominal Distances
#' Bhattacharyya distance
#'
#' Bhattacharyya distance core function
#'
#' @param x Matrix
#' @param adj Small quantity added to avoid indefinite log(0) values. DEFAULT=0.001
#'
#' @return Distance R object
#' @export
BhattacharyyaDist <- function(x, adj = 0.01){
    nRow <- nrow(x)
    myP <- ncol(x)
    compInd <- lapply(sort(unique(as.vector(x))),
                      function(k) diag((x == k) %*% t(x == k))/myP)
    names(compInd) <- sort(unique(as.vector(x)))

    myCoords <- utils::combn(1:nRow, m = 2)
    BC <- sapply(compInd,
                 function(x) utils::combn(x, m = 2,
                                   FUN=function(x){
                                       preSq <- sqrt(prod(x))
                                   }))
    if(length(dim(BC)) == 2)
        BC <- -log(apply(BC, 1, sum) + adj)
    else
        BC <- -log(sum(BC) + adj)

    outMat <- matrix(rep(0, nRow^2), nrow = nRow )
    for(k in 1:ncol(myCoords)){
        outMat[myCoords[1, k], myCoords[2, k]] <- BC[k]
    }
    outMat <- outMat + t(outMat)
    outMat[outMat < 0] <- 0 ## remove adjustment from identical distribs
    stats::as.dist(outMat)
}

#' Bhattacharyya's distance (wrapper)
#'
#' Wrapper of `BhattacharyyaDist`
#'
#' @param X Matrix
#'
#' @return Distance R object
#' @export
dissbhattacharyya  <-  function (X) {

    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    BhattacharyyaDist(x=X)
}

#' Chi-square distance
#'
#' Chi-square distance core function
#'
#' @param x Matrix
#'
#' @return Distance R object
#' @export
ChisqDist <- function(x){
        nRow <- nrow(x)
        myP <- ncol(x)
        ## Count incidences of all categories
        ## List of Frequencies of categories per subject
        compInd <- lapply(sort(unique(as.vector(x))),
                          function(k) diag((x == k) %*% t(x == k)))
        names(compInd) <- sort(unique(as.vector(x)))

        myCoords <- utils::combn(1:nRow, m = 2)
        SqDiff <- sapply(compInd,
                         function(x) utils::combn(x, m = 2,
                                           FUN=function(x){
                                               preSq <- (diff(x)^2)/sum(x)
                                               ifelse(is.nan(preSq), 0, preSq)
                                           }))

        if(is.numeric(SqDiff) & !any(class(SqDiff) == 'matrix'))
            SqDiff <- sum(SqDiff)
        else
            SqDiff <- apply(SqDiff, 1, sum)

        outMat <- matrix(rep(0, nRow^2), nrow =nRow )
        for(k in 1:ncol(myCoords)){
            outMat[myCoords[1, k], myCoords[2, k]] <- SqDiff[k]
        }
        outMat <- outMat + t(outMat)
        stats::as.dist(outMat/myP)
}

#' Chi-square distance (wrapper)
#'
#' Wrapper of `ChisqDist`
#'
#' @param X Matrix
#'
#' @return Distance R object
#' @export
disschisquare  <-  function (X) {

    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }
    ChisqDist(X)
}

#' Cramer's V modified pairwise vector function based on the function found in lsr package
#'
#' This is simple wrapper of the usual chisq.test function.
#' This is actually an adjusted version of the pi = sqrt(Chisq2/N)
#' guaranteeing that values are within 0 (no association) and 1 (association)
#'
#' @param x vector of size n
#' @param y vector of size n
#'
#' @return numerical value
#' @export
cramersVmod <- function(x, y){
   if(identical(x, y)){
      return(1)
   }else if(length(unique(x)) == 1 | length(unique(y)) == 1 ){
      return(0)
   }else if(all(is.na(x)) | all(is.na(y))){
      return(NA)
   }else{
      test <- stats::chisq.test(x=x, y=y, correct=FALSE)
      chi2 <- test$statistic
      N <- sum(test$observed)
      k <- min(dim(test$observed))
      V <- sqrt(chi2/(N * (k - 1)))
      names(V) <- NULL
      return(V)
   }
}

#' Cramer's V distance
#'
#' Cramer's V core function
#'
#' @param X matrix
#'
#' @return Distance matrix
#' @export
CramerV <- function(X){
    myR = nrow(X)
    matOut <- matrix(NA, nrow=myR, ncol=myR)

    utils::combn(myR, 2,
          function(x){
              ## aux <- suppressWarnings(lsr::cramersV(X[x[1], ], X[x[2], ]))
              aux <- suppressWarnings(cramersVmod(X[x[1], ], X[x[2], ]))
              matOut[x[1], x[2]] <<- aux
          })
    matOut <- t(matOut)
    stats::as.dist(matOut)
}

#' Cramer's V distance (wrapper)
#'
#' Wrapper of `CramerV`
#'
#' @param X Matrix
#'
#' @return Distance R object
#' @export
disscramerv <-  function (X) {

    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    1 - CramerV(X)
}

#' Hamming distance wrapper function
#'
#' Function based on cultevo's package implementation
#'
#' @param X matrix
#'
#' @return Distance matrix
#' @export
disshamming <- function (X){

   out <- cultevo::hammingdists(X)/ncol(X)
   return(out)
}

#' Hellinger distance
#'
#' Hellinger distance core function
#' @param x matrix
#'
#' @return Distance matrix
#' @export
HellingerDist <- function(x){
    nRow <- nrow(x)
    myP <- ncol(x)
    compInd <- lapply(sort(unique(as.vector(x))),
                      function(k) rowSums(x == k) )
   names(compInd) <- sort(unique(as.vector(x)))

    myCoords <- utils::combn(1:nRow, m = 2)
    He <- sapply(compInd,
                 function(x) utils::combn(x, m = 2,
                                   FUN=function(x){
                                       preSq <- (diff(sqrt(x)))^2
                                   }))
    if(length(dim(He)) == 2)
        He <- (1/sqrt(2*myP))*sqrt(apply(He, 1, sum))
    else
        He <- (1/sqrt(2*myP))*sqrt(sum(He))

    outMat <- matrix(rep(0, nRow^2), nrow = nRow )
    for(k in 1:ncol(myCoords)){
        outMat[myCoords[1, k], myCoords[2, k]] <- He[k]
    }
    outMat <- outMat + t(outMat)
    outMat[outMat < 0] <- 0 ## remove adjustment from identical distribs
    stats::as.dist(outMat)
}

#' Hellinger distance (wrapper)
#'
#' Wrapper of `HellingerDist`
#'
#' @param X Matrix
#'
#' @return Distance R object
#' @export
disshellinger  <-  function (X) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    HellingerDist(X)
}

#' ####################
#' Ordinal distances ##
#' ####################

#' Absolute distance
#'
#' Absolute distance for ordinal data
#' @param X Numerical matrix with integer/ordinal values.
#'
#' @return Distance matrix
#' @export

dissabs <- function(X){

   stopifnot(is.integer(X))
   n <- nrow(X)

   ecdfMat <- apply(X,
                    2,
                    function(w){
                       stats::ecdf(w)(w)
                       }
                    )
   stats::dist(ecdfMat, method = 'manhattan')

}

#' KS ordinal pairwise distance
#'
#' KS-ordinal pairwise function
#' @param X Numerical matrix with integer/ordinal values.
#'
#' @return Distance matrix
#' @export
dissks <- function(X){

   stopifnot(is.integer(X))

   n <- nrow(X)
   p <- ncol(X)
   dist_mat <- matrix(NA, nrow = n, ncol = n)
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))
   for (i in 1:(n-1)) {
      for (j in (i+1):n) {
         lvls <- sort(unique(c(X[i, ], X[j, ]))) ## This step is unavoidable!!
         # Fx <- stats::ecdf(X[i, ])(lvls)
         # Fy <- stats::ecdf(X[j, ])(lvls)
         Fx <- ecdfL[[i]](lvls)
         Fy <- ecdfL[[j]](lvls)
         dist_mat[i, j] <- max(abs(Fx - Fy))
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   return(stats::as.dist(dist_mat))
}

#' Podani distance wrapper function
#'
#' Function based on FD's implementation
#'
#' @param X Numerical matrix with integer/ordinal values.
#'
#' @return Distance matrix
#' @export
disspodani <- function (X){
   stopifnot(is.integer(X))

   out <- FD::gowdis(X, ord = 'podani')
   return(out)
}

#' Spearman distance wrapper function
#'
#' Function based on bioDist's package implementation.
#'
#' @param X Numerical matrix with integer/ordinal values.
#'
#' @return Distance matrix
#' @export
dissspearman <- function (X){
   ## out <- bioDist::spearman.dist(X, abs = FALSE)
   out <- bioDist::spearman.dist(X, abs = TRUE)
   if(any(is.na(out))){
      out <- as.matrix(out)
      out[is.na(out)] <- 0 ## Constant vectors
      out <- stats::as.dist(out)
   }
   return(out)
}

#' Kendall's Tau distance wrapper function
#'
#' Function based on bioDist's package implementation.
#'
#' @param X Numerical matrix with integer/ordinal values.
#'
#' @return Distance matrix
#' @export
disstau <- function (X){
   ## out <- bioDist::tau.dist(X, abs = FALSE)
   out <- bioDist::tau.dist(X, abs = TRUE)

   if(any(is.na(out))){
      out <- as.matrix(out)
      out[is.na(out)] <- 0 ## Constant vectors
      out <- stats::as.dist(out)
   }
   return(out)
}

#' Wasserstein's DS distance function
#'
#' Estimates eCDF based on the union of observed values
#'
#' @param X numeric Matrix with integer values.
#' @param p numeric Order of the distance. p=1 by default.
#' @return Distance matrix
#' @export
#'
dissWassDS <- function(X, p=1){
   stopifnot(is.integer(X))

   n <- nrow(X)
   dist_mat <- matrix(NA, nrow = n, ncol = n)
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))
   for (i in 1:(n-1)) {
      for (j in (i+1):n) {
         lvls <- sort(unique(c(X[i, ], X[j, ])))
         Fx <- ecdfL[[i]](lvls)
         Fy <- ecdfL[[j]](lvls)
         dist_mat[i, j] <- (mean((abs(Fx - Fy))^p))^(1/p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   return(stats::as.dist(dist_mat))
}

#' Wasserstein's Ranked-based distance function
#'
#' Estimates eCDF based on the union of observed values
#'
#' @param X numeric Matrix with integer values.
#' @param p numeric Order of the distance. p=1 by default.
#' @return Distance matrix
#' @export
#'
dissWassR <- function(X, p=1){
   stopifnot(is.integer(X))

   n <- nrow(X)
   dist_mat <- matrix(NA, nrow = n, ncol = n)
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))
   for (i in 1:(n-1)) {
      for (j in (i+1):n) {
         dist_mat[i, j] <- mean(abs(sort(rank(X[i, ])) - sort(rank(X[j, ]))^p))^(1/p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   return(stats::as.dist(dist_mat))
}

#' Wasserstein's KS distance function
#'
#' Assumes known and identical number of ordinal categories `maxM`
#'
#' @param X numeric Matrix with integer values.
#' @param maxM integer Number of ordinal values which will range from 1 to `maxM`
#' @param p numeric Order of the distance. p=1 by default
#' @return Distance matrix
#' @export
#'
dissWassKS <- function(X, maxM, p=1){
   ## Based on all available values
   ## eCDF based
   stopifnot(is.integer(X))

   n <- nrow(X)
   dist_mat <- matrix(NA, nrow = n, ncol = n)
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))
   for (i in 1:(n-1)) {
      for (j in (i+1):n) {
         lvls <- 1:maxM
         Fx <- ecdfL[[i]](lvls)
         Fy <- ecdfL[[j]](lvls)
         dist_mat[i, j] <- (mean((abs(Fx - Fy))^p))^(1/p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   return(stats::as.dist(dist_mat))
}

#' Wasserstein's distance function
#'
#' Wrapper for the proposed Wasserstein's distance types
#'
#' @param X numerical matrix containing integer/ordinal values.
#' @param type character Select between three WD implementations: `KS` (Known Support),
#' `DS` (Data Support), `R` (Rank).
#' @param m integer Specifies the number of nominal categories for `type='KS'` option.
#' @param p numeric Order of the distance. p=1 by default
#' @return Distance matrix
#' @export
#'
dissWass <- function(X, type='DS', m=NULL, p=1){
   n <- nrow(X)
   dist_mat <- matrix(0, n, n)

   if(type == 'DS'){
      dist_mat <- dissWassDS(X = X, p = p)

   }else if(type == 'KS'){
      m <- as.integer(m)
      dist_mat <- dissWassKS(X = X, maxM = m, p = p)

   }else if(type == 'R'){
      dist_mat <- dissWassR(X = X, p = p)

   }
   return(dist_mat)

}


