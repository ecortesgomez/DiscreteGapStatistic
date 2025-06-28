#' Simulate Data
#'
#' A function to simulate data based on a multinomial vector parameter vector or
#' a list of parameter vectors.
#'
#' @param N Integer. Number of observations.
#' @param nQ Integer. Number of questions.
#' @param pi Numeric vector. Vector of probabilities adding up to 1.
#' Alternatively, pi can be list of vectors as previously described with length equal to `nQ`.
#' This case, notice that that the vectors within the list can be different.
#' The order of the pi vectors in the list will be reflected in the resulting column names.
#' If `dataClass = 'ordinal'` it is required that vector names of `pi` be integers and decimals or other numeric values should be avoided.
#' If `dataClass = nominal` and pi vector names are numerical, these will remain characters.
#' @param dataClass Character. Either 'nominal' or 'ordinal'.
#' @param seed Integer. Numerical seed for the RNG.
#'
#' @returns N x nQ matrix with simulated categories distributed according to vector pi
#' @export
#'
#' @examples
#' Pix <- setNames(c(0.1, 0.2, 0.3, 0.4, 0), paste0('a', 1:5))
#' X <- SimData(N=10, nQ=5, Pix, dataClass = 'nominal')
#' head(X)
#'
#' Piy <- setNames(c(0.3, 0.2, 0.5), paste0('b', 1:3))
#' Y <- SimData(N=10, nQ=3, Piy, dataClass = 'nominal')
#' head(Y)
#'
#' PiZ <- list(x1 = Pix, y1 = Piy, y2 = Piy)
#' Z <- SimData(N = 10, nQ = length(PiZ), PiZ)
#'
#' Piw <- setNames(Piy, 1:3)
#' W <- SimData(N=10, nQ=3, Piw, dataClass = 'ord')
#' head(W)

SimData <- function(N, nQ, pi,
                    dataClass = 'nominal',
                    seed = NULL){

   if(!is.null(seed))
      set.seed(seed)

   if(dataClass == 'nominal' | grepl('^nom', dataClass, ignore.case=TRUE)){

      if(is.numeric(pi)){ ## Just a numerical vector
         ## stopifnot(sum(pi) == 1 )
         stopifnot( (sum(pi) %>% round(digits = 5)) == 1 )

         out <- matrix(sample(names(pi), N*nQ, replace=TRUE, prob = pi),
                       nrow = N,
                       dimnames = list(paste0('s', 1:N),
                                       paste0('q', 1:nQ)))
         return(out)

      }else if(is.list(pi) & length(pi) == nQ){

         ## stopifnot(all(lapply(pi, sum) == 1))
         stopifnot(all(lapply(pi, function(x) sum(x) %>% round(digits = 5)) == 1 ))

         out <-lapply(pi,
                      function(y) sample(x = names(y),
                                         size = N,
                                         replace = TRUE,
                                         prob = y)) %>%
            do.call(what = cbind)
         rownames(out) <- paste0('s', 1:N)

         if(all( class(names(pi)) == 'character' ) )
            colnames(out) <- names(pi)
         else
            colnames(out) <- paste0('Q', 1:nQ)

         return(out)
      }else{
         message(paste0('Either pi is numeric vector with no names ',
                        'or pi is a list with numeric vectors of length no equal to nQ.',
                        'If pi is a list, make sure that length(pi) == nQ'))
      }
   }else if(dataClass == 'ordinal' | grepl('^ord', dataClass, ignore.case=TRUE)){

      if(is.numeric(pi) ){ ## Just a vector
         ## stopifnot(sum(pi) == 1)
         stopifnot( (sum(pi) %>% round(digits = 5)) == 1 )

         out <- matrix(sample(names(pi) %>% as.integer,
                              N*nQ,
                              replace=TRUE,
                              prob = pi),
                       nrow = N,
                       dimnames = list(paste0('s', 1:N),
                                       paste0('q', 1:nQ)))
         return(out)

      }else if(is.list(pi) & length(pi) == nQ){
         stopifnot(all(lapply(pi, sum) == 1))
         out <-lapply(pi,
                      function(y) sample(x = names(y) %>% as.integer,
                                         size = N,
                                         replace = TRUE,
                                         prob = y)) %>%
            do.call(what = cbind)
         rownames(out) <- paste0('s', 1:N)

         if(all( class(names(pi)) == 'character' ) )
            colnames(out) <- names(pi)
         else
            colnames(out) <- paste0('Q', 1:nQ)

         return(out)
      }else{
         message(paste0('Either pi is numeric vector with no names ',
                        'or pi is a list with numeric vectors of length no equal to nQ.',
                        'If pi is a list, make sure that length(pi) == nQ'))
      }
   }
}
