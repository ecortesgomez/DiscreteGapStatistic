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
   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n)
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

.wass_R <- function(X, p = 1) {
   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n)

   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         dist_mat[i, j] <- mean(abs(sort(rank(X[i, ])) - sort(rank(X[j, ])))^p)^(1 / p)
         dist_mat[j, i] <- dist_mat[i, j]
      }
   }
   stats::as.dist(dist_mat)
}

.wass_KS <- function(X, maxM, p = 1) {
   n <- nrow(X)
   dist_mat <- matrix(NA_real_, nrow = n, ncol = n)
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
