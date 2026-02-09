# Ordinal distances -----------------------------------------------------------

.abs_ordinal_dist <- function(X) {
   .assert_integer_matrix(X)
   nCol <- ncol(X)

   ecdfMat <- apply(X, 2, function(w) stats::ecdf(w)(w))
   stats::dist(ecdfMat, method = "manhattan") / nCol
}

.ks_ordinal_dist <- function(X) {
   .assert_integer_matrix(X)

   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n)
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))

   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         lvls <- sort(unique(c(X[i, ], X[j, ])))
         Fx <- ecdfL[[i]](lvls)
         Fy <- ecdfL[[j]](lvls)
         dist_mat[i, j] <- max(abs(Fx - Fy))
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   stats::as.dist(dist_mat)
}

.podani_dist <- function(X) {
   .assert_integer_matrix(X)
   FD::gowdis(X, ord = "podani")
}

.spearman_dist <- function(X) {
   out <- bioDist::spearman.dist(X, abs = TRUE)
   if (any(is.na(out))) {
      m <- as.matrix(out)
      m[is.na(m)] <- 0
      out <- stats::as.dist(m)
   }
   out
}

.tau_dist <- function(X) {
   out <- bioDist::tau.dist(X, abs = TRUE)
   if (any(is.na(out))) {
      m <- as.matrix(out)
      m[is.na(m)] <- 0
      out <- stats::as.dist(m)
   }
   out
}
