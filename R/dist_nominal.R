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
   if (!requireNamespace("nomclust", quietly = TRUE)) {
      stop("Distance metric '", d, "' requires the 'nomclust' package. Install it with install.packages('nomclust').",
           call. = FALSE)
   }

   df <- .nomclust_prepare(X, d)
   metric <- sub("^nmcl_", "", d)

   # Prefer exported-value lookup (works even if nomclust is not attached)
   fn <- tryCatch(getExportedValue("nomclust", metric),
                  error = function(e) NULL)

   if (is.null(fn) || !is.function(fn)) {
      stop("Unknown nomclust distance '", d, "'. ",
           "Expected 'nmcl_' followed by an exported function name from the 'nomclust' package (e.g., 'nmcl_goodall1').",
           call. = FALSE)
   }

   out <- fn(df)

   # Ensure output matches the rest of DGSA distances (a 'dist' object)
   if (!inherits(out, "dist")) out <- stats::as.dist(out)
   out
}

# Nominal distances -----------------------------------------------------------

.hamming_dist <- function(X) {
   .use_rcpp_or_fallback(
      .as_dist_with_names(HammingDist_rcpp(x = X), X),
      function() .HammingDist_R(X),
      "hamming"
   )
}

.bhattacharyya_dist <- function(X, offset = 0.01) {
   .use_rcpp_or_fallback(
      .as_dist_with_names(BhattacharyyaDist_rcpp(x = X, offset = offset), X),
      function() .rename_dist(.BhattacharyyaDist_R(X, adj = offset), X),
      "bhattacharyya"
   )
}

.chisq_dist <- function(X) {
   .use_rcpp_or_fallback(
      .as_dist_with_names(ChisqDist_rcpp(x = X), X),
      function() .rename_dist(.ChisqDist_R(X), X),
      "chisquare"
   )
}

.hellinger_dist <- function(X) {
   .use_rcpp_or_fallback(
      .as_dist_with_names(HellingerDist_rcpp(x = X), X),
      function() .HellingerDist_R(X),
      "hellinger"
   )
}

.cramerv_dist <- function(X) {
   .use_rcpp_or_fallback(
      .as_dist_with_names(CramerVDist_rcpp(x = X), X),
      function() .CramerVDist_R(X),
      "cramerV"
   )
}

# Keep current fallbacks, but give them consistent internal names.
# Your existing .BhattacharyyaDist_R and .ChisqDist_R can be kept as-is,
# except that .BhattacharyyaDist_R should accept adj = offset.

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

.HammingDist_R <- function(X) {
   n <- nrow(X)
   p <- ncol(X)
   D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
   for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
         D[i, j] <- sum(X[i, ] != X[j, ]) / p
         D[j, i] <- D[i, j]
      }
   }
   stats::as.dist(D)
}

.HellingerDist_R <- function(X) {
   nRow <- nrow(X)
   myP <- ncol(X)

   compInd <- lapply(sort(unique(as.vector(X))), function(k) rowSums(X == k))
   names(compInd) <- sort(unique(as.vector(X)))

   myCoords <- utils::combn(seq_len(nRow), m = 2)
   He <- sapply(compInd, function(v) {
      utils::combn(v, m = 2, FUN = function(xx) (diff(sqrt(xx)))^2)
   })

   if (length(dim(He)) == 2) {
      He <- (1 / sqrt(2 * myP)) * sqrt(apply(He, 1, sum))
   } else {
      He <- (1 / sqrt(2 * myP)) * sqrt(sum(He))
   }

   outMat <- matrix(0, nRow, nRow, dimnames = list(rownames(X), rownames(X)))
   for (k in seq_len(ncol(myCoords))) outMat[myCoords[1, k], myCoords[2, k]] <- He[k]
   outMat <- outMat + t(outMat)
   outMat[outMat < 0] <- 0
   stats::as.dist(outMat)
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

.CramerVDist_R <- function(X) {
   myR <- nrow(X)
   matOut <- matrix(NA_real_, nrow = myR, ncol = myR, dimnames = list(rownames(X), rownames(X)))

   utils::combn(myR, 2, function(idx) {
      aux <- suppressWarnings(.cramersVmod(X[idx[1], ], X[idx[2], ]))
      matOut[idx[1], idx[2]] <<- aux
   })

   matOut <- t(matOut)
   stats::as.dist(1 - matOut)
}
