# stats.R - functions related to random number generation and statistical fitting.


# dPL - probability density function for unbounded power-law distribution
# pPL - probability distribution function for unbounded power-law distribution
# rPL - random numbers from an unbounded power-law distribution
# dPLB - probability density function for bounded power-law distribution
# pPLB - probability distribution function for bounded power-law distribution
# rPLB - random numbers from a bounded power-law distribution

# Statistical functions:
##' Bounded and unbounded power-law distributions
##'
##' .. content for \details{} ..
##' @title
##' @param x
##' @param b
##' @param xmin
##' @return
##' @author
dPL = function(x = 1, b = -2, xmin = 1)
  {
  # Computes probability density function for an unbounded power-law
  #  (Pareto) distribution
  #
  # Args:
  #   x: vector of values to compute the density function
  #   b: exponent of probability density function, b < -1
  #   xmin: minimum bound of the distribution, xmin > 0
  #
  # Returns:
  #   vector of probability density corresponding to vector x
  #
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in dPL")
    C = - (b+1) / xmin^(b+1)
    y = 0 * x     # so have zeros where x < xmin
    y[x >= xmin] = C * x[x >= xmin]^b
    return(y)
  }

pPL = function(x = 10, b = -2, xmin = 1)
  {
  # Computes probability distribution function, P(X <= x),  for an
  #   unbounded power-law (Pareto) distribution
  #
  # Args:
  #   x: vector of values at which to compute the distribution function
  #   b: exponent of probability density function, b < -1
  #   xmin: minimum bound of the distribution, xmin > 0
  #
  # Returns:
  #   vector of probability distribution values P(X <= x) corresponding to x
  #
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in qPL")
    y = 0 * x     # so have zeros where x < xmin
    y[x >= xmin] = 1 - (x[x >= xmin]/xmin)^(b+1)
    return(y)
  }

rPL = function(n = 1, b = -2, xmin = 1)
  {
  # Computes random numbers from an unbounded power-law (Pareto) distribution
  #
  # Args:
  #   n: number of random numbers in sample. If 'length(n) > 1', the length is
  #       taken to be the number required.
  #   b: exponent of probability density function, b < -1
  #   xmin: minimum bound of the distribution, xmin > 0
  #
  # Returns:
  #   vector of length n of independent random draws from the distribution
  #
  # Uses the inverse method (e.g. p1215 of Edwards 2008, Journal
  #   of Animal Ecology, 77:1212-1222).
  #
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in rPL")
    u = runif(n)
    y = xmin * ( 1 - u ) ^ (1/(b+1))
    return(y)
  }

# Bounded power-law distribution:

dPLB = function(x = 1, b = -2, xmin = 1, xmax = 100)
  {
  # Computes probability density function for a bounded power-law
  #  (Pareto) distribution
  #
  # Args:
  #   x: vector of values to compute the density function
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  # Returns:
  #   vector of probability density corresponding to vector x
  #
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

pPLB = function(x = 10, b = -2, xmin = 1, xmax = 100)
  {
  # Computes probability distribution function, P(X <= x),  for a
  #   bounded power-law (Pareto) distribution
  #
  # Args:
  #   x: vector of values at which to compute the distribution function
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  # Returns:
  #   vector of probability distribution values P(X <= x) corresponding to x
  #
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y = 0 * x     # so have zeros where x < xmin
    y[x > xmax] = 1  # 1 for x > xmax
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

rPLB = function(n = 1, b = -2, xmin = 1, xmax = 100)
  {
  # Computes random numbers from a bounded power-law distribution
  #
  # Args:
  #   n: number of random numbers in sample. If 'length(n) > 1', the length is
  #       taken to be the number required.
  #   b: exponent of probability density function
  #   xmin: minimum bound of the distribution, xmin > 0
  #   xmax: maximum bound of the distribution, xmax > xmin
  #
  # Returns:
  #   vector of length n of independent random draws from the distribution
  #
  # Uses the inverse method (e.g. p1215 of Edwards 2008, Journal
  #   of Animal Ecology, 77:1212-1222).
  #
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
#  function and random number generator for unbounded and bounded power-law
#  distributions. Associated with paper in Methods in Ecology and Evolution:
#  http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12641/full
#
#  andrew.edwards@dfo-mpo.gc.ca   or  andrew.edwards.dfo@gmail.com
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
