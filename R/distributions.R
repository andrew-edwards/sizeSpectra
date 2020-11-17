# distributions.R - functions related to random number generation for bounded
# and unbounded power-law distributions. Also the biomass distribution function
# (MEE equations A.4 and A.8).

##' Bounded and unbounded power-law distributions
##'
##' For the unbounded and bounded power-law distributions (PL and PLB
##' respectively), probability density function (`dPL` and `dPLB`),
##' cumulative distribution function P(X <= x) (`pPL` and `pPLB`),
##' and random generation of values (`rPL` and `rPLB`), with exponent `b`, minimum `xmin` and maximum
##' (for bounded distribution) `xmax` as described in Edwards et al. (2017, Methods in Ecology and
##' Evolution, 8:57-67). Random generation uses the inverse method (e.g. p1215 of Edwards
##' 2008, Journal of Animal Ecology, 77:1212-1222). Unbounded distribution
##' included for completeness but is not used in remaining code.
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

##' Biomass distribution function from MEE equations A.4 and A.8
##'
##' Total biomass between two body masses WILL MAKE MORE GENERAL `x1` and `x2`, assuming a bounded power-law
##' distribution of body masses between `xmin` and `xmax` and a given value of
##' exponent `b`, with a total of `n` individuals.
##' Given by MEE equations A.4 and A.8 (with the integration calculated between
##' general `x1` and `x2` rather than `xmin` and `x`).
##' @return
##' @export
##' @author Andrew Edwards
##' @examples
##' @donttest{
##' @
##' @}
pBiomass <- function(b = -2,
                     xmin = 1,
                     xmax = 100,
                     n = 1000,
                     x1 = 1,
                     x2 = 100){
  if(xmin <= 0 | xmin >= xmax | x1 < xmin | x2 > xmax | x1 >= xmax | n <= 0){
    stop("Parameters out of bounds in pBiomass")
  }

  if(b != -1){
    C <- (b+1) / ( xmax^(b+1) - xmin^(b+1) )
  } else {
    C <- 1/ ( log(xmax) - log(xmin) )
  }

  if(b != -2){
    biomass <- n * C * (x2^(b+2) - x1^(b+2)) / (b + 2)
  } else {
    biomass <- n * C * (log(x2) - log(x1))
  }

  return(biomass)
}
