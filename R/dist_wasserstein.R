# Wasserstein family ----------------------------------------------------------

.wasserstein_dispatch <- function(X, d) {
   .assert_integer_matrix(X)

   # Existing patterns from your original dispatcher:
   #  - "wasserstein ... KS<m>"
   #  - "wasserstein ... DS"
   #  - "wasserstein ... R"
   if (grepl("wass.+KS[0-9]+", d)) {
      myM <- sub("wass.+KS", "", d)
      myM <- as.integer(myM)
      return(.wass_KS(X, maxM = myM))
   }

   if (grepl("wass.+DS$", d)) return(.wass_DS(X))
   if (grepl("wass.+R$", d))  return(.wass_R(X))

   stop("Unknown Wasserstein spec: '", d, "'.", call. = FALSE)
}

.wass_DS <- function(X, p = 1) {
   .assert_integer_matrix(X)
   .use_rcpp_or_fallback(
      .as_dist_with_names(WassersteinDS_rcpp(x = X, power = p), X),
      function() .WassDS_R(X, p = p),
      "wassersteinDS"
   )
}

.wass_KS <- function(X, maxM, p = 1) {
   .assert_integer_matrix(X)
   .use_rcpp_or_fallback(
      .as_dist_with_names(WassersteinKS_rcpp(x = X, maxM = maxM, power = p), X),
      function() .WassKS_R(X, maxM = maxM, p = p),
      "wassersteinKS"
   )
}

.wass_R <- function(X, p = 1) {
   .assert_integer_matrix(X)
   .use_rcpp_or_fallback(
      .as_dist_with_names(WassersteinRank_rcpp(x = X, power = p), X),
      function() .WassRank_R(X, p = p),
      "wassersteinR"
   )
}

.WassDS_R <- function(X, p = 1) {
   ## .wass_DS previously
   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(rownames(X), rownames(X)))
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))

   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         lvls <- sort(unique(c(X[i, ], X[j, ])))
         Fx <- ecdfL[[i]](lvls)
         Fy <- ecdfL[[j]](lvls)
         dist_mat[i, j] <- (mean((abs(Fx - Fy))^p))^(1 / p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   stats::as.dist(dist_mat)
}

.WassKS_R <- function(X, maxM, p = 1) {
   ## .wass_KS
   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(rownames(X), rownames(X)))
   ecdfL <- apply(X, 1, function(w) stats::ecdf(w))

   lvls <- seq_len(maxM)
   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         Fx <- ecdfL[[i]](lvls)
         Fy <- ecdfL[[j]](lvls)
         dist_mat[i, j] <- (mean((abs(Fx - Fy))^p))^(1 / p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   stats::as.dist(dist_mat)
}

.WassRank_R <- function(X, p = 1) {
   ## .wass_R previously
   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(rownames(X), rownames(X)))

   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         dist_mat[i, j] <- mean(abs(sort(rank(X[i, ])) - sort(rank(X[j, ])))^p)^(1 / p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   stats::as.dist(dist_mat)
}


