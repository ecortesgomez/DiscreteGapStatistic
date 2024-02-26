#' Function invoking discrete distance functions
#'
#' @param X Matrix where rows are the observations and columns are discrete features
#' @param d Name of distance. Distances available: bhattacharyya, chisquare, cramerV, hamming and hellinger
#' @param na.rm Remove NAs default=TRUE
#'
#' @return R distance object
#' @export
#'
#' @examples
#' # X = rbind(matrix(paste0("a", rpois(7*5, 1)), nrow=5),
#' #           matrix(paste0("a", rpois(7*5, 3)), nrow=5))
#' # distancematrix(X = X, d = "hellinger")
distancematrix <- function (X, d, na.rm=TRUE){

    if(d == 'bhattacharyya')
        return(dissbhattacharyya(X, na.rm))
    if(d == 'chisquare')
        return(disschisquare(X, na.rm))
    if (d == "cramerV")
         return(disscramerv(X, na.rm))
    if (d == "hamming") ## General version
        return(disshamming(X, na.rm))
    if (d == "hellinger")
        return(disshellinger(X, na.rm))

    stop("Distance metric ", d, " not available")
}

#' Bhattacharyya distance core function
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

#' Bhattacharyya's wrapper Function
#' @param X Matrix
#' @param na.rm Remove NAs default=TRUE
#'
#' @return Distance R object
#' @export
dissbhattacharyya  <-  function (X, na.rm = TRUE) {

    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    BhattacharyyaDist(x=X)
}

#' Chi-square distance core function
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

#' Chi-square distance wrapper function
#' @param X Matrix
#' @param na.rm logical
#'
#' @return Distance R object
#' @export
disschisquare  <-  function (X, na.rm = TRUE) {
    ## chisquare2 distance for discrete distributions

    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    ChisqDist(X)
}

#' Cramer's V modified pairwise vector function based on the function found in lsr package
#' This is simple wrapper of the usual chisq.test fun
#' This is actually an adjusted version of the pi = sqrt(Chisq2/N)
#' guaranteeing that values are within 0 (no association) and 1 (association)
#' @param x vector of size n
#' @param y vector of size n
#'
#' @return numerical value
#' @export

cramersVmod <- function(x, y){
   if(length(unique(x)) == 1 | length(unique(y)) == 1 )
      return(0)
   else if(all(is.na(x)) | all(is.na(y)))
      return(NA)
   else
      lsr::cramersV(x, y)
}


#' Cramer's V core function
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

#' Cramer's V distance wrapper function
#' @param X Matrix
#' @param na.rm logical
#'
#' @return Distance R object
#' @export
disscramerv <-  function (X, na.rm = TRUE) {

    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    1 - CramerV(X)
}

#' Hellinger distance core function
#' @param x matrix
#'
#' @return Distance matrix
#' @export
HellingerDist <- function(x){
    nRow <- nrow(x)
    myP <- ncol(x)
    compInd <- lapply(sort(unique(as.vector(x))),
                      function(k) diag((x == k) %*% t(x == k))/myP)
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

#' Hellinger's distance wrapper Function
#' @param X Matrix
#' @param na.rm logical
#'
#' @return Distance R object
#' @export
disshellinger  <-  function (X, na.rm = TRUE) {
    if (!is.matrix(X)) {
        stop(paste(sQuote("X"), "not a matrix"))
    }

    HellingerDist(X)
}

#' Hamming distance wrapper function
#' Function based on cultevo's package implementation
#' @param X matrix
#' @param na.rm logical
#'
#' @return Distance matrix
#' @export
disshamming <- function (X, na.rm = TRUE){
    ## Euclidean Distance from the rows of a matrix

    out <- cultevo::hammingdists(X)
    return(out)
}
