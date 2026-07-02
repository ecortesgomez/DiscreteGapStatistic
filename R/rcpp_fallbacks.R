# Helpers ---------------------------------------------------------------------

.dgsa_use_rcpp <- function() {
   isTRUE(getOption("DiscreteGapStatistic.use_rcpp", TRUE))
}

.as_dist_with_names <- function(D, X) {
   D <- as.matrix(D)
   rn <- rownames(X)
   if (!is.null(rn)) dimnames(D) <- list(rn, rn)
   stats::as.dist(D)
}

.use_rcpp_or_fallback <- function(expr, fallback, distance_name) {
   if (!.dgsa_use_rcpp()) return(fallback())

   tryCatch(
      expr,
      error = function(e) {
         warning(
            "Rcpp implementation for '", distance_name,
            "' failed; using R fallback. Original error: ",
            conditionMessage(e),
            call. = FALSE
         )
         fallback()
      }
   )
}

.rename_dist <- function(D, X) {
   D <- as.matrix(D)

   rn <- rownames(X)
   if (!is.null(rn)) {
      dimnames(D) <- list(rn, rn)
   }

   stats::as.dist(D)
}
