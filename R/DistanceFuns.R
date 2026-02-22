#' Calculate categorical distance matrix for discrete data
#'
#' Dispatcher for discrete distance functions (nominal + ordinal).
#'
#' @param X Matrix where rows are observations and columns are discrete features.
#' @param d Character scalar naming the distance.
#'
#' @return An object of class \code{dist}.
#' @export
#'
#' @useDynLib DiscreteGapStatistic, .registration = TRUE
#'
#' @examples
#' X <- rbind(matrix(paste0("a", rpois(7*5, 1)), nrow=5),
#'            matrix(paste0("a", rpois(7*5, 3)), nrow=5))
#' distancematrix(X = X, d = "hellinger")
distancematrix <- function(X, d) {
   if (!is.matrix(X)) stop("'X' must be a matrix.", call. = FALSE)
   if (!is.character(d) || length(d) != 1L) stop("'d' must be a single character string.", call. = FALSE)

   # Registry-based dispatch for "simple" names
   fun <- .distance_registry()[[d]]
   if (!is.null(fun)) return(fun(X))

   # nomclust family
   # if (grepl("^nmcl_.+", d)) return(row.names = rownames(X, d))
   if (grepl("^nmcl_.+", d)) return(.nomclust_distance(X, d))

   # wasserstein family (legacy strings include 'wasserstein' + suffix)
   if (grepl("wasserstein", d)) return(.wasserstein_dispatch(X, d))

   stop("Distance metric '", d, "' not available.", call. = FALSE)
}

.distance_registry <- function() {
   list(
      # Nominal
      bhattacharyya = .bhattacharyya_dist,
      chisquare     = .chisq_dist,
      cramerV       = .cramerv_dist,
      hamming       = .hamming_dist,
      hellinger     = .hellinger_dist,

      # Ordinal
      absolute      = .abs_ordinal_dist,
      ks            = .ks_ordinal_dist,
      podani        = .podani_dist,
      spearman      = .spearman_dist,
      tau           = .tau_dist
   )
}

.assert_integer_matrix <- function(X) {
   if (!is.matrix(X) || !is.integer(X)) {
      stop("This distance requires an *integer* matrix. ",
           "Convert with storage.mode(X) <- 'integer' if appropriate.",
           call. = FALSE)
   }
   invisible(TRUE)
}
