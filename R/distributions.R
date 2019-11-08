# distributions.R - functions related to random number generation.

##' Bounded and unbounded power-law distributions
##'
##' For the unbounded and bounded power-law distributions (PL and PLB
##' respectively), probability density function (`dPL` and `dPLB`),
##' cumulative distribution function P(X <= x) (`pPL` and `pPLB`),
##' and random generation of values (`rPL` and `rPLB`), with exponent `b`, minimum `xmin` and maximum
##' (for bounded distribution) `xmax` as described in Edwards et al. (2017, Methods in Ecology and
##' Evolution, 8:57-67). Random generation uses the inverse method (e.g. p1215 of Edwards
##' 2008, Journal of Animal Ecology, 77:1212-1222). Unbounded distribution
##' included for completeness TODO: may not get used in remaining code.
##' @param x vector of values to compute the density and distribution functions.
##' @param n number of random numbers to be generated (if `length(n) > 1` then
##' generate `length(n)` values)
##' @param b exponent of the distribution (must be <-1 for unbounded)
##' @param xmin minimum bound of the distribution, `xmin > 0`
##' @param xmax maximum bound for bounded distribution, `xmax > xmin`
##' @return `dPL` and `dPLB` return vector of probability density values
##' corresponding to `x`. `pPL` and `pPLB` return vector of cumulative
##' distribution values P(X <= x) corresponding to `x`. `rPL` and `rPLB` return
##' a vector (of length `n`) of independent random draws from the distribution.
##' @name Distributions
NULL

##' @rdname Distributions
##' @export
dPL <- function(x = 1, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in dPL")
    C <- - (b+1) / xmin^(b+1)
    y <- 0 * x     # so have zeros where x < xmin
    y[x >= xmin] <- C * x[x >= xmin]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPL <- function(x = 10, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in qPL")
    y <- 0 * x     # so have zeros where x < xmin
    y[x >= xmin] <- 1 - (x[x >= xmin]/xmin)^(b+1)
    return(y)
  }

##' @rdname Distributions
##' @export
rPL <- function(n = 1, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in rPL")
    u <- runif(n)
    y <- xmin * ( 1 - u ) ^ (1/(b+1))
    return(y)
  }

##' @rdname Distributions
##' @export
dPLB <- function(x = 1, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in dPLB")
    if(b != -1)
        { C <- (b+1) / ( xmax^(b+1) - xmin^(b+1) )
        } else
        { C <- 1/ ( log(xmax) - log(xmin) )
        }
    y <- 0 * x     # so have zeros where x < xmin or x > xmax
    y[x >= xmin & x <= xmax] <- C * x[x >= xmin & x <= xmax]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPLB <- function(x = 10, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y <- 0 * x        # so have zeros where x < xmin
    y[x > xmax] <- 1
    if(b != -1)
        {  xmintobplus1 <- xmin^(b+1)
           denom <- xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] <-
                                        ( x[x >= xmin & x <= xmax]^(b + 1) -
                                          xmintobplus1 ) / denom
        } else
        {  logxmin <- log(xmin)
           denom <- log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)
  }

##' @rdname Distributions
##' @export
rPLB <- function(n = 1, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in rPLB")
    u <- runif(n)
    if(b != -1)
        { y <- ( u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
        } else
        { y <- xmax^u * xmin^(1-u)
        }
    return(y)
  }
