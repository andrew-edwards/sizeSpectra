# distributions.R - functions related to random number generation and statistical fitting.

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

# dPL - probability density function for unbounded power-law distribution
# pPL - probability distribution function for unbounded power-law distribution
# rPL - random numbers from an unbounded power-law distribution
# dPLB - probability density function for bounded power-law distribution
# pPLB - probability distribution function for bounded power-law distribution
# rPLB - random numbers from a bounded power-law distribution

##' @rdname Distributions
##' @export
dPL = function(x = 1, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in dPL")
    C = - (b+1) / xmin^(b+1)
    y = 0 * x     # so have zeros where x < xmin
    y[x >= xmin] = C * x[x >= xmin]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPL = function(x = 10, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in qPL")
    y = 0 * x     # so have zeros where x < xmin
    y[x >= xmin] = 1 - (x[x >= xmin]/xmin)^(b+1)
    return(y)
  }

##' @rdname Distributions
##' @export
rPL = function(n = 1, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in rPL")
    u = runif(n)
    y = xmin * ( 1 - u ) ^ (1/(b+1))
    return(y)
  }

##' @rdname Distributions
##' @export
dPLB = function(x = 1, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in dPLB")
    if(b != -1)
        { C = (b+1) / ( xmax^(b+1) - xmin^(b+1) )
        } else
        { C = 1/ ( log(xmax) - log(xmin) )
        }
    y = 0 * x     # so have zeros where x < xmin or x > xmax
    y[x >= xmin & x <= xmax] = C * x[x >= xmin & x <= xmax]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPLB = function(x = 10, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y = 0 * x        # so have zeros where x < xmin
    y[x > xmax] = 1
    if(b != -1)
        {  xmintobplus1 = xmin^(b+1)
           denom = xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] =
               ( x[x >= xmin & x <= xmax]^(b + 1) - xmintobplus1 ) / denom
        } else
        {  logxmin = log(xmin)
           denom = log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)
  }

##' @rdname Distributions
##' @export
rPLB = function(n = 1, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in rPLB")
    u = runif(n)
    if(b != -1)
        { y = ( u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
        } else
        { y = xmax^u * xmin^(1-u)
        }
    return(y)
  }

# PLBfunctions.r - functions to be sourced, including density, distribution
#
# CONTENTS
#
# Statistical functions:
# negLL.PLB - negative log-likelihood function for PLB model
# sum.bins - calculate the total sum of values within a bin
# Llin.method - fitting data using the Llin method
# LBNbiom.method - fitting data using the LBbiom and LBNbiom methods
# LBmizbinsFuns - calculate the bin breaks for the LBmiz method from mizer
# log2bins - construct bins that double in size and encompass the data
# eightMethods - computes exponents for a dataset using all eight methods, for
#  a second manuscript
# negLL.PLB.binned - negative log-likelihood function for PLB when data are
#  only available in binned form
#
# Plotting functions:
# lm.line - plot straight line of lm fit but restricted to the x values
# gap.barplot.cust - customised version of Jim Lemon's gap.barplot for
#  histograms with a break in an axis
# qqtab - constructs automated LaTeX code for tables of quantiles
# confPlot - plotting of the confidence intervals for Figure 4
# histAxes - histogram axes for histogram plots of estimated b values
#  (Figure 3 and others)
# histAxes2 - histAxes adapted for fitting3rep-n10000.r, for n=10,000 sample size
# logTicks - add axes and tick marks to a log-log plot to represent
#  unlogged values (e.g. Figures 2(h) and 6(b))
# legJust - add legend to a plot
#
#  2nd Sept 2014 onwards.
