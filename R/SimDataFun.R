#' Simulate Data
#'
#' @param N Integer. Number of observations.
#' @param nQ Integer. Number of questions.
#' @param pi Numeric vector. Vector of probabilities adding up to 1. Names of elements can be characters.
#'
#' @returns N x nQ matrix with simulated categories distributed according to vector pi
#' @export
#'
#' @examples
#' Pix <- setNames(c(0.1, 0.2, 0.3, 0.4, 0), paste0('a', 1:5))
#' X <- SimData(N=10, nQ=5, Pix)
#' head(X)
#'
#' Piy <- setNames(c(0.6, 0.1, 0.1, 0.1, 0.1), letters[1:5])
#' Y <- SimData(N=10, nQ=3, Piy)
#' head(Y)
SimData <- function(N, nQ, pi){
   out <- matrix(sample(names(pi), N*nQ, replace=TRUE, prob = pi),
                 nrow = N,
                 dimnames = list(paste0('s', 1:N),
                                 paste0('q', 1:nQ)))
    out
}
