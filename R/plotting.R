# Functions for plotting

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

##' Plot bounded straight line of results of `lm()` fit
##'
##' Plot a straight line between lowest and highest values of `x.vector`
##'  based on `lm()` results, to use instead of `abline(lm.output)`
##'  which extends beyond the fitted data.
##'
##' @param x.vector values of x (that have been fitted to), fitted line will be
##'   bound by the range of `x.vector`
##' @param lm.results output from an `lm()` fit, with intercept `lm.results$coeff[1]`
  #     and gradient `lm.results$coeff[2]`
##' @param ... arguments to lines
##' @return Adds a line to an existing plot
##' @export
##' @author Andrew Edwards
lm.line = function(x.vector, lm.results, ...)
  {
     intercept = lm.results$coeff[1]
     grad = lm.results$coeff[2]
     lines(c( min(x.vector), max(x.vector)), c(grad * min(x.vector) + intercept,
      grad * max(x.vector) + intercept), ...)
}
