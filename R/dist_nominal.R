# Nominal distances -----------------------------------------------------------

# Bhattacharyya (prefer Rcpp, keep R fallback)
.bhattacharyya_dist <- function(X, offset = 1e-8) {
   D <- BhattacharyyaDist_rcpp(x = X, offset = offset)
   stats::as.dist(D)
}

# Optional pure-R fallback (kept internal)
.BhattacharyyaDist_R <- function(x, adj = 0.01) {
   nRow <- nrow(x)
   myP <- ncol(x)
   compInd <- lapply(sort(unique(as.vector(x))),
                     function(k) diag((x == k) %*% t(x == k)) / myP)
   names(compInd) <- sort(unique(as.vector(x)))

   myCoords <- utils::combn(seq_len(nRow), m = 2)
   BC <- sapply(compInd,
                function(v) utils::combn(v, m = 2, FUN = function(xx) sqrt(prod(xx))))

   if (length(dim(BC)) == 2) {
      BC <- -log(apply(BC, 1, sum) + adj)
   } else {
      BC <- -log(sum(BC) + adj)
   }

   outMat <- matrix(0, nRow, nRow)
   for (k in seq_len(ncol(myCoords))) {
      outMat[myCoords[1, k], myCoords[2, k]] <- BC[k]
   }
   outMat <- outMat + t(outMat)
   outMat[outMat < 0] <- 0
   stats::as.dist(outMat)
}

# Chi-square (prefer Rcpp, keep R fallback)
.chisq_dist <- function(X) {
   D <- ChisqDist_rcpp(X)
   stats::as.dist(D)
}

# Optional pure-R fallback (kept internal)
.ChisqDist_R <- function(x) {
   nRow <- nrow(x)
   myP <- ncol(x)

   compInd <- lapply(sort(unique(as.vector(x))),
                     function(k) diag((x == k) %*% t(x == k)))
   names(compInd) <- sort(unique(as.vector(x)))

   myCoords <- utils::combn(seq_len(nRow), m = 2)
   SqDiff <- sapply(compInd,
                    function(v) utils::combn(v, m = 2, FUN = function(xx) {
                       preSq <- (diff(xx)^2) / sum(xx)
                       ifelse(is.nan(preSq), 0, preSq)
                    }))

   if (is.numeric(SqDiff) && !inherits(SqDiff, "matrix")) SqDiff <- sum(SqDiff) else SqDiff <- apply(SqDiff, 1, sum)

   outMat <- matrix(0, nRow, nRow)
   for (k in seq_len(ncol(myCoords))) {
      outMat[myCoords[1, k], myCoords[2, k]] <- SqDiff[k]
   }
   outMat <- outMat + t(outMat)
   stats::as.dist(outMat / myP)
}

# Cramer's V distance = 1 - association
.cramerv_dist <- function(X) {
   myR <- nrow(X)
   matOut <- matrix(NA_real_, nrow = myR, ncol = myR)

   utils::combn(myR, 2, function(idx) {
      aux <- suppressWarnings(.cramersVmod(X[idx[1], ], X[idx[2], ]))
      matOut[idx[1], idx[2]] <<- aux
   })

   matOut <- t(matOut)
   stats::as.dist(1 - matOut)
}

.cramersVmod <- function(x, y) {
   if (identical(x, y)) return(1)
   if (length(unique(x)) == 1L || length(unique(y)) == 1L) return(0)
   if (all(is.na(x)) || all(is.na(y))) return(NA_real_)

   test <- stats::chisq.test(x = x, y = y, correct = FALSE)
   chi2 <- unname(test$statistic)
   N <- sum(test$observed)
   k <- min(dim(test$observed))
   sqrt(chi2 / (N * (k - 1)))
}

.hamming_dist <- function(X) {
   cultevo::hammingdists(X) / ncol(X)
}

.hellinger_dist <- function(X) {
   nRow <- nrow(X)
   myP <- ncol(X)

   compInd <- lapply(sort(unique(as.vector(X))),
                     function(k) rowSums(X == k))
   names(compInd) <- sort(unique(as.vector(X)))

   myCoords <- utils::combn(seq_len(nRow), m = 2)
   He <- sapply(compInd,
                function(v) utils::combn(v, m = 2, FUN = function(xx) (diff(sqrt(xx)))^2))

   if (length(dim(He)) == 2) {
      He <- (1 / sqrt(2 * myP)) * sqrt(apply(He, 1, sum))
   } else {
      He <- (1 / sqrt(2 * myP)) * sqrt(sum(He))
   }

   outMat <- matrix(0, nRow, nRow)
   for (k in seq_len(ncol(myCoords))) {
      outMat[myCoords[1, k], myCoords[2, k]] <- He[k]
   }
   outMat <- outMat + t(outMat)
   outMat[outMat < 0] <- 0
   stats::as.dist(outMat)
}

# nomclust dispatch -----------------------------------------------------------

.nomclust_prepare <- function(X, d) {
   if (grepl("^nmcl_good.*", d)) {
      df <- data.frame(X, check.names = FALSE)
      df[] <- lapply(df, factor)
      return(df)
   }
   as.data.frame(X)
}

.nomclust_distance <- function(X, d) {
   df <- .nomclust_prepare(X, d)
   metric <- sub("^nmcl_", "", d)
   fn <- get(metric, envir = asNamespace("nomclust"), inherits = FALSE)
   fn(df)
}
