# Functions for plotting and customised functions for output for latex tables

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

##' Customising `plotrix::gap.barplot` for a histogram with a gap in y-axis
##'
##'  For Figure 1 of MEE paper, to make a histogram (barplot) with a gap in the y-axis.
##'   Customising `gap.barplot()` from the package plotrix by Jim Lemon and others
##'   Several options here are customised for the particular plot (and to
##'   change a few of the defaults in gap.barplot) so the code would
##'   require some modifiying to use more generally.
##'
##' This function modifies `plotrix::gap.barplot()`, between 2nd Sept 2014 and
##' finalised here in October 2019. `plotrix` was written by Jim Lemon and
##' others and is available from CRAN at [https://cran.r-project.org/web/packages/plotrix/index.html].
##'
##' @param y vector of data values
##' @param gap range of values to be left out
##' @param xaxlab labels for the x axis ticks
##' @param xtics position of the x axis ticks
##' @param yaxlab labels for the y axis ticks
##' @param ytics position of the y axis ticks
##' @param midpoints TODO
##' @param breakpoints TODO
##' @param xlim optional x limits for the plot
##' @param ylim optional y limits for the plot
##' @param xlab label for the x axis
##' @param ylab label for the y axis
##' @param horiz whether to have vertical or horizontal bars
##' @param col color(s) in which to plot the values
##' @param N TODO: to do with ticks I think
##' @param ... arguments passed to 'barplot'.
##' @return Barplot with a gap in the y-axis
##' @export
##' @author Andrew Edwards
gap.barplot.cust = function (y,
                             gap,
                             xaxlab,
                             xtics,
                             yaxlab,
                             ytics,
                             midpoints,
                             breakpoints,
                             xlim,
                             ylim = NA,
                             xlab = NULL,
                             ylab = NULL,
                             horiz = FALSE,
                             col = NULL,
                             N = n, ...)
    {
    if (missing(y))
        stop("y values required")
    # x <- 1:length(y)        # Original
    x <- midpoints            # AE adding; midpoints is defined in raw1infexamp.r
    if (missing(gap))
        stop("gap must be specified")
    if (is.null(ylab))
        ylab <- deparse(substitute(y))
    if (is.null(col))
        col <- color.gradient(c(0, 1), c(0, 1, 0), c(1, 0), length(y))
    else if (length(col) < length(y))
        rep(col, length.out = length(y))
    littleones <- which(y <= gap[1])
    bigones <- which(y >= gap[2])
    if (any(y > gap[1] & y < gap[2]))
        warning("gap includes some values of y")
    gapsize <- gap[2] - gap[1]
    if (missing(xaxlab))
        xaxlab <- as.character(x)
    # xlim <- c(min(x) - 0.4, max(x) + 0.4)    # Original
    xlim <- xlim
    if (is.na(ylim[1]))
        ylim <- c(min(y), max(y) - gapsize)
    if (missing(ytics))
        ytics <- pretty(y)
    if (missing(yaxlab))
        yaxlab <- ytics
    littletics <- which(ytics < gap[1])
    bigtics <- which(ytics >= gap[2])
    halfwidth <- min(diff(x))/2
    if (horiz) {                   # AE not editing anything for horizontal
        if (!is.null(xlab)) {
            tmplab <- xlab
            xlab <- ylab
            ylab <- tmplab
        }
        plot(0, xlim = ylim, ylim = xlim, ylab = ylab, axes = FALSE,
            type = "n", ...)
        plot.lim <- par("usr")
        botgap <- ifelse(gap[1] < 0, gap[1], plot.lim[1])
        box()
        axis(2, at = x, labels = xaxlab)
        axis(1, at = c(ytics[littletics], ytics[bigtics] - gapsize),
            labels = c(ytics[littletics], ytics[bigtics]))
        rect(botgap, x[y < gap[1]] - halfwidth, y[y < gap[1]],
            x[y < gap[1]] + halfwidth, col = col[y < gap[1]])
        rect(botgap, x[bigones] - halfwidth, y[bigones] - gapsize,
            x[bigones] + halfwidth, col = col[bigones])
        plotrix::axis.break(1, gap[1], style = "gap")
    }
    else {                        # AE editing here for vertical bars
        plot(0, xlim = xlim, ylim = ylim, ylab = ylab, axes = FALSE, xlab=xlab,
            type = "n", ...)               # AE added xlab=xlab
        plot.lim <- par("usr")
        botgap <- ifelse(gap[1] < 0, gap[1], plot.lim[3])
        box()
        # axis(1, at = x, labels = xaxlab)  # - original was x=1,2,3,...,58, is
                                            #  now the midpoints, but I want
                                            #  ticks at breakpoints,
                                            #  hmadeup$breaks
        axis(1, breakpoints, tcl=-0.2, labels=rep("", length(breakpoints)))
                                                # short ones at breaks
        # axis(1, at=seq(0, 600, 50), labels = rep("", 13), mgp=c(1.7,0.6,0))
                                             # long ticks where I want
        axis(1, at=seq(0, 1000, 100),
             labels = seq(0, 1000, 100), mgp=c(1.7,0.6,0))
                  # long ticks labelled where I want, not automated though
        axis(2, at = c(ytics[littletics], ytics[bigtics] - gapsize),
            labels = c(ytics[littletics], ytics[bigtics]), mgp=c(1.7,0.6,0))
                       # AE adding mgp, these are the ones with numbers on
        # Short tics:
        bottomShortTics = seq(0, gap[1], 2)        # below gap
        axis(2, bottomShortTics, tcl=-0.2,
             labels=rep("", length(bottomShortTics)))
        topShortTics = seq(gap[2], N, 2)-(gap[2]-gap[1])  # above gap
        axis(2, topShortTics, tcl=-0.2,
             labels=rep("", length(topShortTics)))
        rect(x[y < gap[1]] - halfwidth, botgap, x[y < gap[1]] +
            halfwidth, y[y < gap[1]], col = col[y < gap[1]])
        rect(x[bigones] - halfwidth, botgap, x[bigones] + halfwidth,
            y[bigones] - gapsize, col = col[bigones])
        plotrix::axis.break(2, gap[1], style = "gap")
    }
    invisible(x)
}

##' Produce LaTeX table code for quantiles of results
##'
##'  Quantile table function, to construct lines of quantiles to copy
##'   into LaTeX. Does `quants[1]`, 50\% (median), mean, `quants[2]` quantiles,
##'   and the \%age of values < the true value.
##'
##' @param xx vector of values to give quantiles for
##' @param dig number of decimal places to give
##' @param true the true value of the quantity being estimated, will
##'     depend on method
##' @param quants quantiles to use, with defaults of 0.25 and 0.75 (25\% and 75\%).
##' @return LaTeX output for table to copy into a .tex file
##' @export
##' @author Andrew Edwards
qqtab = function(xx,
                 dig=2,
                 true=b.known,
                 quants = c(0.25, 0.75))
  {
      paste( c( prettyNum(round(quantile(xx, quants[1]),
                                  digits=dig), big.mark=","),
         " & ", prettyNum(round(quantile(xx, 0.50), digits=dig),
                                  big.mark=","),
         " & ", prettyNum(round(mean(xx), digits=dig),
                                  big.mark=","),
         " & ", prettyNum(round(quantile(xx, quants[2]), digits=dig),
                                  big.mark=","),
         " & ", prettyNum(round(sum(xx < true)/length(xx)*100, digits=0),
                                  big.mark=",")),
         sep="", collapse="")             # , "\\%"),
  }
