test_that("Rcpp nominal distances match R fallbacks", {
  set.seed(20260702)
  X <- matrix(sample(LETTERS[1:5], 15 * 8, replace = TRUE), nrow = 15)
  rownames(X) <- paste0("s", seq_len(nrow(X)))

  old_opt <- getOption("DiscreteGapStatistic.use_rcpp", TRUE)
  on.exit(options(DiscreteGapStatistic.use_rcpp = old_opt), add = TRUE)

  check_one <- function(fun, tol = 1e-12) {
    options(DiscreteGapStatistic.use_rcpp = TRUE)
    D_cpp <- fun(X)
    options(DiscreteGapStatistic.use_rcpp = FALSE)
    D_R <- fun(X)

    expect_s3_class(D_cpp, "dist")
    expect_s3_class(D_R, "dist")
    expect_equal(as.matrix(D_cpp), as.matrix(D_R), tolerance = tol)
  }

  check_one(.hamming_dist)
  check_one(.chisq_dist)
  check_one(.bhattacharyya_dist)
  check_one(.hellinger_dist)
  check_one(.cramerv_dist, tol = 1e-10)
})

test_that("Rcpp ordinal and Wasserstein distances match R fallbacks", {
  set.seed(20260702)
  X <- matrix(sample(1:7, 12 * 10, replace = TRUE), nrow = 12)
  storage.mode(X) <- "integer"
  rownames(X) <- paste0("s", seq_len(nrow(X)))

  old_opt <- getOption("DiscreteGapStatistic.use_rcpp", TRUE)
  on.exit(options(DiscreteGapStatistic.use_rcpp = old_opt), add = TRUE)

  check_one <- function(fun, tol = 1e-12) {
    options(DiscreteGapStatistic.use_rcpp = TRUE)
    D_cpp <- fun(X)
    options(DiscreteGapStatistic.use_rcpp = FALSE)
    D_R <- fun(X)

    expect_s3_class(D_cpp, "dist")
    expect_s3_class(D_R, "dist")
    expect_equal(as.matrix(D_cpp), as.matrix(D_R), tolerance = tol)
  }

  check_one(.ks_ordinal_dist)
  check_one(function(z) .wass_DS(z, p = 1))
  check_one(function(z) .wass_KS(z, maxM = 7, p = 1))
  check_one(function(z) .wass_R(z, p = 1))
})

test_that("basic distance invariants hold", {
  Xn <- matrix(c("A", "B", "A", "B", "A", "A"), nrow = 3, byrow = TRUE)
  D <- as.matrix(.hamming_dist(Xn))
  expect_true(all(diag(D) == 0))
  expect_equal(D, t(D))
  expect_true(all(D >= 0))

  Xo <- matrix(c(1L, 2L, 3L, 1L, 2L, 3L, 3L, 3L, 3L), nrow = 3, byrow = TRUE)
  Dks <- as.matrix(.ks_ordinal_dist(Xo))
  expect_true(all(diag(Dks) == 0))
  expect_equal(Dks, t(Dks))
  expect_true(all(Dks >= 0))
})
