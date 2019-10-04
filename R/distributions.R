# distributions.R - functions related to random number generation.

##' Bounded and unbounded power-law distributions
##'
##' For the unbounded and bounded power-law distributions (PL and PLB
##' respectively), probability density function (`dPL` and `dPLB`),
##' cumulative distribution function [P(X <= x)] (`pPL` and `pPLB`),
##' and random generation of values (`rPL` and `rPLB`), with exponent `b`, minimum `x_min` and maximum
##' (for bounded distribution) `x_max` as described in Edwards et al. (2017, Methods in Ecology and
##' Evolution, 8:57-67). Random generation uses the inverse method (e.g. p1215 of Edwards
##' 2008, Journal of Animal Ecology, 77:1212-1222). Unbounded distribution
##' included for completeness [TODO: may not get used in remaining code].
##' @param x vector of values to compute the density and distribution functions.
##' @param n number of random numbers to be generated (if `length(n) > 1` then
##' generate `length(n)` values)
##' @param b exponent of the distribution (must be <-1 for unbounded)
##' @param x_min minimum bound of the distribution, `x_min > 0`
##' @param x_max maximum bound for bounded distribution, `x_max > x_min`
##' @return `dPL` and `dPLB` return vector of probability density values
##' corresponding to `x`. `pPL` and `pPLB` return vector of cumulative
##' distribution values P(X <= x) corresponding to `x`. `rPL` and `rPLB` return
##' a vector (of length `n`) of independent random draws from the distribution.
##' @name Distributions
NULL

##' @rdname Distributions
##' @export
dPL <- function(x = 1, b = -2, x_min = 1)
  {
    if(b >= -1 | x_min <= 0) stop("Parameters out of bounds in dPL")
    C <- - (b+1) / x_min^(b+1)
    y <- 0 * x     # so have zeros where x < x_min
    y[x >= x_min] <- C * x[x >= x_min]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPL <- function(x = 10, b = -2, x_min = 1)
  {
    if(b >= -1 | x_min <= 0) stop("Parameters out of bounds in qPL")
    y <- 0 * x     # so have zeros where x < x_min
    y[x >= x_min] <- 1 - (x[x >= x_min]/x_min)^(b+1)
    return(y)
  }

##' @rdname Distributions
##' @export
rPL <- function(n = 1, b = -2, x_min = 1)
  {
    if(b >= -1 | x_min <= 0) stop("Parameters out of bounds in rPL")
    u <- runif(n)
    y <- x_min * ( 1 - u ) ^ (1/(b+1))
    return(y)
  }

##' @rdname Distributions
##' @export
dPLB <- function(x = 1, b = -2, x_min = 1, x_max = 100)
  {
    if(x_min <= 0 | x_min >= x_max) stop("Parameters out of bounds in dPLB")
    if(b != -1)
        { C <- (b+1) / ( x_max^(b+1) - x_min^(b+1) )
        } else
        { C <- 1/ ( log(x_max) - log(x_min) )
        }
    y <- 0 * x     # so have zeros where x < x_min or x > x_max
    y[x >= x_min & x <= x_max] <- C * x[x >= x_min & x <= x_max]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPLB <- function(x = 10, b = -2, x_min = 1, x_max = 100)
  {
    if(x_min <= 0 | x_min >= x_max) stop("Parameters out of bounds in pPLB")
    y <- 0 * x        # so have zeros where x < x_min
    y[x > x_max] <- 1
    if(b != -1)
        {  x_mintobplus1 <- x_min^(b+1)
           denom <- x_max^(b+1) - x_mintobplus1
           y[x >= x_min & x <= x_max] <-
                                        ( x[x >= x_min & x <= x_max]^(b + 1) -
                                          x_mintobplus1 ) / denom
        } else
        {  logx_min <- log(x_min)
           denom <- log(x_max) - logx_min
           y[x >= x_min & x <= x_max] =
               ( log( x[x >= x_min & x <= x_max] ) - logx_min ) / denom
        }
    return(y)
  }

##' @rdname Distributions
##' @export
rPLB <- function(n = 1, b = -2, x_min = 1, x_max = 100)
  {
    if(x_min <= 0 | x_min >= x_max) stop("Parameters out of bounds in rPLB")
    u <- runif(n)
    if(b != -1)
        { y <- ( u*x_max^(b+1) +  (1-u) * x_min^(b+1) ) ^ (1/(b+1))
        } else
        { y <- x_max^u * x_min^(1-u)
        }
    return(y)
  }
