# Functions for plotting and customised functions for output for latex and
#  Rmarkdown tables

# lm.line - plot straight line of lm fit but restricted to the x values
# gap.barplot.cust - customised version of Jim Lemon's gap.barplot for
#  histograms with a break in an axis
# qqtab - constructs automated row of LaTeX or Rmarkdown code for tables of
#  quantiles
# confPlot - plotting of the confidence intervals for Figure 4
# histAxes - histogram axes for histogram plots of estimated b values
#  (Figure 3 and others)
# histAxes2 - histAxes adapted for fitting3rep-n10000.r, for n=10,000 sample size
# logTicks - add axes and tick marks to a log-log plot to represent
#  unlogged values (e.g. Figures 2(h) and 6(b))
# legJust - add legend to a plot
# MLEmid.MLEbin.hist - one figure with eight histograms for each of MLEmid and
#  MLEbin methods and four binning types
# MLEmid.MLEbin.conf - one figure with eight confidence interval plots for each
#  of MLEmid and MLEbin methods and four binning types
# MLEmid.MLEbin.table - make dataframe of results from MLEmid and MLEbin methods
#  and four binning types
# timeSerPlot - plot time series of estimated *b* with confidence intervals
# ISD_bin_plot - recommended plotting for binned data (MEPS Figures 7 and
#  S.5-S.34).

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
##'   Several default options here are customised for the particular plot (and to
##'   change a few of the defaults in gap.barplot) so the code would
##'   require some modifiying to use more generally.
##'
##' This function modifies `plotrix::gap.barplot()`, between 2nd Sept 2014 and
##' finalised here in October 2019. `plotrix` was written by Jim Lemon and
##' others and is available from CRAN at https://cran.r-project.org/web/packages/plotrix/index.html.
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
                             gap = c(9,980),
                             xaxlab,
                             xtics,
                             yaxlab,
                             ytics = c(seq(0, 8, by=4), seq(980, 988, by=4)),
                             midpoints,
                             breakpoints,
                             xlim,
                             ylim = c(0, 17),  # max = max(y) - gap[2] + gap[1]
                                               # + some space
                             xlab = expression(paste("Values, ", italic(x))),
                             ylab = "Count in each bin",
                             horiz = FALSE,
                             col = NULL,
                             N = n, ...)
  {
    if (missing(y))
        stop("y values required")
    x <- midpoints
    # if (missing(gap))
    #    stop("gap must be specified")
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
    # if (missing(ytics))
    #    ytics <- pretty(y)
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

##' Produce row for a dataframe, LaTeX or Rmarkdown table code for quantiles of results.
##'
##'  Quantile table function, to construct row for a dataframe or a line of quantiles to copy
##'   into either LaTeX or render directly in an Rmarkdown document. Does `quants[1]`,
##'   50\% (median), mean, `quants[2]` quantiles,
##'   and the \%age of values < the true value.
##'
##' @param xx vector of values to give quantiles for
##' @param dig number of decimal places to give
##' @param true the true value of the quantity being estimated, will
##'     depend on method
##' @param quants quantiles to use, with defaults of 0.25 and 0.75 (25\% and
##'   75\%).
##' @param type `data.frame` produces a row to put into a dataframe,
##'   `latex` gives output to copy into a `.tex` file, `Rmd` gives Rmarkdown
##'    output (see vignette TODO).
##' @return row of a data.frame, line of LaTeX output for table to copy into a
##'   .tex file, or line of Rmarkdown to copy into an Rmarkdown document
##'
##' @export
##' @author Andrew Edwards
qqtab = function(xx,
                 dig = 2,
                 true = NA,
                 quants = c(0.25, 0.75),
                 type = "Rmd")
{
  if (!(type %in% c("data.frame", "latex", "Rmd"))) stop("Invalid type in qqtab().")

  if(type == "latex"){
     noquote(
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
             sep="", collapse=""))             # , "\\%"),
  } else if(type == "Rmd") {
    noquote(
      paste( c( prettyNum(round(quantile(xx, quants[1]),
                                  digits=dig), big.mark=","),
             " | ", prettyNum(round(quantile(xx, 0.50), digits=dig),
                                  big.mark=","),
             " | ", prettyNum(round(mean(xx), digits=dig),
                                  big.mark=","),
             " | ", prettyNum(round(quantile(xx, quants[2]), digits=dig),
                                  big.mark=","),
             " | ", prettyNum(round(sum(xx < true)/length(xx)*100, digits=0),
                                  big.mark=",")),
               sep=" ", collapse=" "))
  } else if(type == "data.frame")
  {
    return(c( prettyNum(round(quantile(xx, quants[1]),
                              digits=dig),
                      big.mark=","),
           prettyNum(round(quantile(xx, 0.50),
                           digits=dig),
                     big.mark=","),
           prettyNum(round(mean(xx), digits=dig),
                     big.mark=","),
           prettyNum(round(quantile(xx, quants[2]),
                           digits=dig),
                     big.mark=","),
           prettyNum(round(sum(xx < true)/length(xx)*100,
                           digits=0),
                     big.mark=",")))
               # sep=" ", collapse=" "))
  }
}


##' Plot confidence intervals of repeated estimates for one method
##'
##' Plot confidence intervals of the repeated estimates of `b` for one
##' method. Gets called eight times to produce Figure 4 of MEE paper and TODO of
##' MEPS paper. Plots
##' horizontal lines for the intervals, colour coded as to whether or not they
##' include the true value of b.
##'
##' @param repConf data frame with columns `confMin` and `confMax`, the
##'   minimum and maximum of the 95\% confidence interval for b, and rows
##'   corresponding to each (TODO a, yes?) simulated data set. When `confPlot` is called, for
##'   some methods `b = slope - 1` or `b = slope - 2`, which gets done in the
##'   call.
##' @param legName legend name for that panel
##' @param b.true true known value of b (as used for simulating the data)
##' @param inCol colour for intervals that `b.true` is inside
##' @param outCol colour for intervals that `b.true` is outside
##' @param vertCol colour of vertical line for `b.true`
##' @param pchVal `pch` for points at endpoints of intervals
##' @param cexVal size of points (default of 0 does not plot them)
##' @param xLim x-axis range
##' @param colourCode colour code the figure (may not completely work TODO)
##' @param vertThick thickness of vertical line for `b.true`
##' @param horizThick thickness of horizontal lines
##' @param thin number of values to thin the values for plotting. Only works for
##'      33 or 99 for now (since they are factors of 9999, the value used in simulations).
##' @param horizLines whether to plot horizontal grey lines or not as example
##'        (not needed now since doing lines for intervals)
##' @param horizLinesOut whether to plot horizontal lines for
##'      intervals for which true value is outside the interval
##' @param horizLinesIn whether to plot horizontal lines for
##'      intervals for which true value is inside the interval
##' @param yLab label for y-axis
##' @param yTicks where to have tickmarks on y-axis
##' @param yLabels whether or not to label tickmarks on y-axis
##' @param vertFirst whether or not to plot vertical line first in order to see
##'      the horizontal lines better (for LCD plot)
##' @param insetVal inset shift for naming the panel
##' @param insetVal2 inset shift for printing observed coverage percentage
##' @param xsmallticks where to put unlabelled small tick marks on x-axis
##' @param legLoc where to put the legend, as the first argument in `legend()`
##' @return plots a panel and returns what is plotted as a sorted data frame TODO
##' @export
##' @author Andrew Edwards
confPlot = function(repConf,
                    legName,
                    b.true = b.known,
                    inCol="darkgrey",
                    outCol="blue",
                    vertCol="red",
                    pchVal = 20,
                    cexVal = 0.0,
                    xLim = NULL,
                    colourCode = TRUE,
                    vertThick = 1,
                    horizThick = 0.3,
                    thin = 33,
                    horizLines = FALSE,
                    horizLinesOut = TRUE,
                    horizLinesIn = TRUE,
                    yLab = "Sample number",
                    yTicks = seq(0, 300, 50),
                    yLabels = TRUE,
                    vertFirst = FALSE,
                    insetVal = c(-0.08, -0.06),
                    insetVal2 = c(-0.08, 0.07),
                    xsmallticks=NULL,
                    legLoc = "topleft")
    {
    if(!colourCode) outCol = inCol
    if(!(thin %in% c(33,99))) stop("Need to edit confPlot if thin not 33 or 99")
    if(is.null(xLim))
        {
            rangeVal = range(repConf)
            xLim = c(floor(rangeVal[1]), ceiling(rangeVal[2]))
                  # min(repConf[,1])), ceiling(max(repConf[,2])))
        }
    repConf = dplyr::mutate(repConf,
      inConf = (b.true > confMin & b.true < confMax))  # does CI cover b.true?

    repConf = dplyr::mutate(repConf, confCol = NA)
    # Couldn't figure this out in dplyr:
    repConf[which(repConf$inConf), "confCol"] = inCol
    repConf[which(!repConf$inConf), "confCol"] = outCol
    repConf[, "num.rep"] = as.numeric(row.names(repConf))
              # To preserve the iteration number in case needed later

    repConf.sort.full = dplyr::arrange(repConf, confMin)   # arranged by min of CI
    sum.inConf = sum(repConf.sort.full$inConf, na.rm=TRUE) /
        dim(repConf.sort.full)[1]      # get NA's for Llin, LT with xmax=10000
    # subset for better plotting:
    if(thin == 33)
        { repConf.sort = repConf.sort.full[c(1, 1+thin*(1:303)), ] }else
        { repConf.sort = repConf.sort.full[c(1, 1+thin*(1:101)), ] }
    # print(head(repConf.sort))
    # repConf.sort[, "num.sorted"] = as.numeric(row.names(repConf.sort))
    #  That doesn't work since it preserves the repConf.sort.full row names
    repConf.sort[, "num.sorted"] = 1:(dim(repConf.sort)[1])
              # To preserve the new row number since needed later
    plot(repConf.sort$confMin, repConf.sort$num.sorted, col=repConf.sort$confCol,
     xlim = xLim, ylim = c(0, 1.02*dim(repConf.sort)[1]),
     pch=pchVal, cex=cexVal,
     xlab = "",
     ylab = yLab, yaxt="n")
    if(vertFirst) {abline(v=b.true, lwd=vertThick, col=vertCol)}
    axis(2, at = yTicks, labels = yLabels, tck=-0.04)
    if(!is.null(xsmallticks))
        { axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)), tcl=-0.2)}
    # xlab = expression(paste("Estimate of slope (or ", italic(b), ")")),
    points(repConf.sort$confMax, repConf.sort$num.sorted,
           col=repConf.sort$confCol, pch=pchVal, cex=cexVal)
    legend(legLoc, legName, bty="n", inset=insetVal)
    legend(legLoc, paste(round(sum.inConf*100, digits = 0), "%", sep=""),
           bty="n", inset=insetVal2)
    # legend("topleft", "hello", bty="n", inset=c(-0.08, -0.2))

    # Plot just the rows for which true value lie outside CI
    if(horizLinesOut)
        {
        outConf = dplyr::filter(repConf.sort, (!inConf))
        for(jj in 1:(dim(outConf)[1]))
          {

            lines(c(outConf$confMin[jj], outConf$confMax[jj]),
              c(outConf$num.sorted[jj], outConf$num.sorted[jj]),
              col=outCol, lwd=horizThick)
          }
        }

    # Plot just the rows for which true value lie inside CI
    if(horizLinesIn)
        {
        in.Conf = dplyr::filter(repConf.sort, (inConf))
        for(jj in 1:(dim(in.Conf)[1]))
          {

            lines(c(in.Conf$confMin[jj], in.Conf$confMax[jj]),
              c(in.Conf$num.sorted[jj], in.Conf$num.sorted[jj]),
              col=inCol, lwd=horizThick)
          }
        }

    # Plot equally spaced horizontal lines between 2.5 and 97.5 values
    #  (though even from the likelihood ratio test for MLE a 95% CI
    #  isn't actually the 2.5-97.5 range, it's a 95% interval).
    if(horizLines)
        {
        # row numbers to use to plot equally spaced horizontal lines
        horLineInd = round(c(2.5, 21.5, 40.5, 59.5, 78.5, 97.5) * 0.01 *
            dim(repConf.sort)[1])
        horiz.lines = repConf.sort[horLineInd,]
        for(jj in 1:length(horLineInd))
            {
              lines(c(horiz.lines$confMin[jj], horiz.lines$confMax[jj]),
              c(horiz.lines$num.sorted[jj], horiz.lines$num.sorted[jj] ),
              col="grey", lwd=horizThick)
            }
        }
    if(!vertFirst) {abline(v=b.true, lwd=vertThick, col=vertCol)}
    return(repConf.sort.full)     # repConf.sort is what's plotted though
}

##' Histogram axes for figures TODO
##'
##' Do the histogram axes in `fitting3rep.r` TODO and later, since almost all
##' panels will have same axes. Not very flexible, use `histAxes2()`
##' to plot up to 10,000. TODO Note that `xbsmallticks` and `xbigticks` need to be
##' prespecified, but
##' they haven't yet been made arguments here.
##'
##' @param yBigTickLab y-axis big ticks to label
##' @param yBigTickNoLab y-axis big ticks to not label
##' @param ySmallTick y-axis small ticks (unlabelled)
##' @param cexAxis font size for axis labels TODO defaults are for
##'   MEE_reproduce_2.Rmd. May have to check other figs
##' @param xsmallticks where to put unlabelled small tick marks on x-axis
##' @param xbigticks where to put big tick marks on x-axis
##' @param vertCol vertCol colour of vertical line for `b.known`
##' @param vertThick thickness of vertical line for `b.known`
##' @param b.known known value of $b$
##' @return Adds axes to existing histogram TODO: check
##' @export
##' @author Andrew Edwards
histAxes = function(yBigTickLab = c(0, 2000, 4000),
                    yBigTickNoLab = seq(0, 6500, 1000),
                    ySmallTick = seq(0, 6500, 500),
                    cexAxis = 0.9,
                    xsmallticks = seq(-3.5, 0.5, 0.1),
                    xbigticks = seq(-3.5, 0.5, 0.5),
                    vertCol = "red",
                    vertThick = 1,
                    b.known = -2
                    )
    {
    axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks
    axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)), tcl=-0.2)
    axis(2, at=yBigTickLab,
         labels = yBigTickLab,
         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled
      axis(2, at=yBigTickNoLab,
         labels = rep("", length(yBigTickNoLab)),
         mgp=c(1.7,0.7,0))  # big ticks unlabelled
      axis(2, at=ySmallTick,
         labels = rep("", length(ySmallTick)), mgp=c(1.7,0.7,0), tcl=-0.2)
                           # small ticks
      abline(v=b.known, col=vertCol, lwd=vertThick)
}

##' Do histogram axes for Figure TODO
##'
##' Do the histogram axes for, in particular, TODO `fitting3rep-n10000.r`,
##' since the larger n sample size means affects the resulting histograms
##' of estimated b. Not very flexible, just doing for the one figure.
##' Copying what is used for Llin method (so could go back and use this function
##'   throughout earlier code).
##'
##' @return Adds axes to existing histogram TODO: check
##' @export
##' @author Andrew Edwards
histAxes2 = function()
    {
      axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0),
           cex.axis=cexAxis)  # big ticks
      axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)),
           tcl=-0.2)
      axis(2, at=c(0, 5000, 10000),
         labels = c(0, 5000, 10000),
         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled
      axis(2, at=seq(0, 10000, 1000),
         labels = rep("", 11),
         mgp=c(1.7,0.7,0))  # big ticks unlabelled
      abline(v=b.known, col=vertCol, lwd=vertThick)
    }


##' Add axes and tick marks to a log-log plot to represent unlogged values
##'
##' Useful because you can then interpret the unlogged, e.g. Figure TODO
##'
##' @param xLim the x limits for the plot (unlogged scale); if NULL then do not
##'   add anything to x-axis
##' @param yLim  the y limits for the plot (unlogged scale); if NULL then
##'   do not add anything to y-axis
##' @param tclSmall size of small tick marks
##' @param xLabelSmall which small tick marks on x-axis to label
##' @param yLabelSmall which small tick marks on y-axis to label
##' @param xLabelBig which big tick marks on the x-axis to label
##'   (when automated they can overlap, so may need to specify)
##' @param mgpVal `mgp` values for axes. See `?par`
##' @return Adds axes and big and small tick marks to the plot. Returns NULL
##' @examples
##' \dontrun{
##' # Adapt the following:   TODO make explicit
##'   plot(..., log="xy", xlab=..., ylab=..., xlim=..., ylim=..., axes=FALSE)
##'   xLim = 10^par("usr")[1:2]
##'   yLim = 10^par("usr")[3:4]
##'   logTicks(xLim, yLim, xLabelSmall = c(5, 50, 500))
##' }
##' @export
##' @author Andrew Edwards
logTicks = function(xLim,
                    yLim = NULL,
                    tclSmall = -0.2,
                    xLabelSmall = NULL,
                    yLabelSmall = NULL,
                    xLabelBig = NULL,
                    mgpVal=c(1.6,0.5,0))
    {
    ll = 1:9
    log10ll = log10(ll)
    box()
    # x axis
    if(!is.null(xLim))               # if NULL then ignore
      {
      # Do enough tick marks to encompass axes:
      xEncompassLog = c(floor(log10(xLim[1])), ceiling(log10(xLim[2])))
      xBig = 10^c(xEncompassLog[1]:xEncompassLog[2])
      # Big unlabelled, always want these:
      axis(1, at= xBig, labels = rep("", length(xBig)), mgp = mgpVal)
      # Big labelled:
      if(is.null(xLabelBig)) { xLabelBig = xBig }
      axis(1, at= xLabelBig, labels = xLabelBig, mgp = mgpVal)
      # axis(1, at=c(1, 10, 100), labels = c(1, 10, 100), mgp=c(1.7,0.7,0))
      # Small unlabelled:
      axis(1, xBig %x% ll, labels=rep("", length(xBig %x% ll)), tcl=tclSmall)
      # Small labelled:
      if(!is.null(xLabelSmall))
          {
          axis(1, at=xLabelSmall, labels=xLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }
    # Repeat for y axis:
    if(!is.null(yLim))               # if NULL then ignore
      {
      # Do enough tick marks to encompass axes:
      yEncompassLog = c(floor(log10(yLim[1])), ceiling(log10(yLim[2])))
      yBig = 10^c(yEncompassLog[1]:yEncompassLog[2])
      # Big labelled:
      axis(2, at= yBig, labels = yBig, mgp = mgpVal)
      # Small unlabelled:
      axis(2, yBig %x% ll, labels=rep("", length(yBig %x% ll)), tcl=tclSmall)
      # Small labelled:
      if(!is.null(yLabelSmall))
          {
          axis(2, at=yLabelSmall, labels=yLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }
}
##' Add legend with right-justification
##'
##' Add legend with right-justification, functionalising Uwe Ligges'
##'   example in `?legend`. Really a way of adding text automatically in
##'   the corner (which `legend()` works out the positioning for).
##'
##' @param textVec text for the legend, one element for each row
##' @param pos position of the legend (not tested for all positions)
##' @param textWidth width to make the text
##' @param inset inset distance
##' @param logxy TRUE if axes are logarithmic
##' @return Adds legend to exiting plot. Returns NULL.
##' @export
##' @examples
##' \dontrun{ TODO check
##'   plot(10:1)
##'  legJust(c("Method", paste("value=", mean(1:10), sep=""), "x=7.0"))
##' }
##' @author Andrew Edwards
legJust = function(textVec,
                   pos="topright",
                   textWidth = "slope=-*.**",
                   inset=0,
                   logxy=FALSE)
    {
      leg = legend(pos,
                   legend = rep(" ", length(textVec)),
                   text.width = strwidth(textWidth),
                   bty="n",
                   inset=inset)
    if(logxy)
    { text(10^(leg$rect$left + leg$rect$w),
           10^(leg$text$y),
           textVec,
           pos = 2) } else
    { text(leg$rect$left + leg$rect$w,
           leg$text$y,
           textVec,
           pos = 2) }
}

##' One figure with eight histograms for each of  MLEmid and MLEbin methods and four binning types
##'
##' The default plot here reproduces Figure 4 of MEPS paper, showing the
##' estimated exponent $b$ for 10,000 simulated data sets, binned using four
##' methods and fitted using the MLEmid and MLEbin methods.
##'
##' @param results.list output list from `MLEbin.simulate()`; see
##'   `?MLEbin.simulate()` for details
##' @param vertCol colour for vertical line at true value of `b`
##' @param vertThick thickness of vertical line at true value of `b`
##' @param xrange range of x values TODO can change to xlimA for consistency
##' @param xbigticks.by increment between big tick marks on x-axis (all labelled)
##' @param xsmallticks.by increment between small tick marks on x-axis
##' @param ylimA range of y values
##' @param yBigTickLab.by increment between labelled big tick marks on y-axis
##' @param yBigTickNoLab.by increment between all big tick marks on y-axis
##' @param ySmallTick.by increment between small tick marks on y-axis
##' @param binwidth bin width for histograms; if NA then gets calculated such
##'   that `b.known` is a midpoint of a bin
##' @param cexAxis font size for axis labels
##' @param legLabMid character vector of length 4 for labelling each panel in
##'   the MLEmid column
##' @param legLabBin character vector of length 4 for labelling each panel in
##'   the MLEBin column
##' @param omi,mfrow,mai,xaxs,yaxs,mgp,cex,inset standard options for `par()`,
##'   defaults are for Figure 4
##' @return figure with eight histograms, one for each combination of binning
##'   type and fitting method.
##' @export
##' @author Andrew Edwards
MLEmid.MLEbin.hist = function(results.list,
                              vertCol = "red",
                              vertThick = 1,
                              xrange = c(-2.4, -1.1),
                              xbigticks.by = 0.4,
                              xsmallticks.by = 0.1,
                              ylimA = c(0, 7000),
                              yBigTickLab.by = 3000,
                              yBigTickNoLab.by = 1000,
                              ySmallTick.by = 500,
                              binwidth = NA,
                              omi = c(0.12, 0.05, 0.2, 0.0),
                              mfrow = c(4, 2),
                              mai = c(0.5, 0.5, 0.0, 0.3),
                              xaxs="i",
                              yaxs="i",
                              mgp = c(2.0, 0.5, 0),
                              cex = 0.8,
                              cexAxis = 0.9,
                              inset = c(0, -0.04),
                              legLabMid = c("(a)", "(c)", "(e)", "(g)"),
                              legLabBin = c("(b)", "(d)", "(f)", "(h)"))
{
  # Extract required components
  MLE.array = results.list$MLE.array
  num.reps = results.list$MLE.array.parameters$num.reps
  b.known = results.list$MLE.array.parameters$b.known
  binTypes = results.list$MLE.array.parameters$binTypes
  binType.name = results.list$MLE.array.parameters$binType.name

  xbigticks = seq(xrange[1], xrange[2], by = xbigticks.by)
  xsmallticks = seq(xrange[1], xrange[2], by = xsmallticks.by)

  yBigTickLab = seq(0, num.reps, yBigTickLab.by)
  yBigTickNoLab = seq(0, num.reps, yBigTickNoLab.by)
  ySmallTick = seq(0, num.reps, ySmallTick.by)

  # Want b.known (-2 for default) to be a midpoint of the Nth bin, which is yy above the minimum.
  #  Say the 20th bin contains -2, so solve ((N+1)w + Nw)/2 = yy
  #  gives  w = 2 * yy / (2*N + 1). Not flexible yet.
  binwidth = 2 * (b.known - xrange[1]) / ( 2 * 10 + 1)
  breakshist =  seq(xrange[1],
                    length = ceiling( (xrange[2] - xrange[1])/binwidth ) + 1,
                    by = binwidth)

  par(omi = omi,
      mfrow = mfrow,
      mai = mai,
      xaxs = xaxs,
      yaxs = yaxs,
      mgp = mgp,
      cex = cex)

  for(binTypeInd in 1:binTypes)
    {
      hist(MLE.array[ ,binTypeInd, "MLEmid"],
           xlim=xrange,
           breaks=breakshist,
           xlab="",
           ylab="",
           main="",
           axes=FALSE,
           ylim = ylimA)
      histAxes(yBigTickLab = yBigTickLab,
               yBigTickNoLab = yBigTickNoLab,
               ySmallTick = ySmallTick,
               cexAxis = cexAxis,
               xbigticks = xbigticks,
               xsmallticks = xsmallticks,
               vertCol = vertCol,
               vertThick = vertThick)
      legend("topright",
             paste(legLabMid[binTypeInd], binType.name[binTypeInd]),
             bty="n",
             inset=inset)

      hist(MLE.array[ , binTypeInd, "MLEbin"],
           xlim=xrange,
           breaks=breakshist,
           xlab="",
           ylab="",
           main="",
           axes=FALSE,
           ylim = ylimA)
      histAxes(yBigTickLab = yBigTickLab,
               yBigTickNoLab = yBigTickNoLab,
               ySmallTick = ySmallTick,
               cexAxis = cexAxis,
               xbigticks = xbigticks,
               xsmallticks = xsmallticks,
               vertCol = vertCol,
               vertThick = vertThick)
      legend("topright",
             paste(legLabBin[binTypeInd], binType.name[binTypeInd]),
             bty="n",
             inset=inset)
  }

  mtext(expression(paste("Estimate of ", italic(b))),
        side=1,
        outer=TRUE,
        line=-1)
  mtext("Frequency",
        side=2,
        outer=TRUE,
        line=-1)

  mtext("    MLEmid                                                MLEbin",
        side=3,
        outer=TRUE,
        line=0)
}

##' One figure with eight confidence interval plots for each of  MLEmid and MLEbin methods and four binning types
##'
##' The default plot here reproduces Figure 5 of MEPS paper, showing the
##' confidence intervals of the exponent $b$ for 10,000 simulated data sets, binned using four
##' methods and fitted using the MLEmid and MLEbin methods.
##'
##' @param results.list output list from `MLEbin.simulate()`; see
##'   `?MLEbin.simulate()` for details
##' @param vertCol colour for vertical line at true value of `b`
##' @param vertThick thickness of vertical line at true value of `b`
##' @param xrange range of x values TODO can change to xlimA for consistency
##' @param xbigticks.by increment between big tick marks on x-axis (all labelled)
##' @param xsmallticks.by increment between small tick marks on x-axis
##' @param legLabMid character vector of length 4 for labelling each panel in
##'   the MLEmid column
##' @param legLabBin character vector of length 4 for labelling each panel in
##'   the MLEBin column
##' @param insetMat matrix of inset values for MLEmid column (Figure 5(e) legend
##'   has to be shifted), each row is the inset for that MLEmid panel
##' @param omi,mfrow,mai,xaxs,yaxs,mgp,cex,inset standard options for `par()`,
##'   defaults are for Figure 4
##' @return figure with eight panels of confidence intervals, one for each combination of binning
##'   type and fitting method.
##' @export
##' @author Andrew Edwards
MLEmid.MLEbin.conf = function(results.list,
                              vertCol = "red",
                              vertThick = 1,
                              xrange = c(-2.4, -1.1),
                              xbigticks.by = 0.4,
                              xsmallticks.by = 0.1,
                              omi = c(0.2, 0.05, 0.22, 0.1),
                              mfrow = c(4, 2),
                              mai = c(0.3, 0.5, 0.08, 0),
                              xaxs = "i",
                              yaxs = "i",
                              mgp = c(2.0, 0.5, 0),
                              cex = 0.8,
                              inset = c(0, -0.04),
                              insetMat = matrix(rep(c(-0.01, -0.04),
                                                    4),
                                                ncol=2,
                                                byrow=TRUE),
                              legLabMid = c("(a)", "(c)", "(e)", "(g)"),
                              legLabBin = c("(b)", "(d)", "(f)", "(h)"))
{
  # Extract required components
  MLEconf.array = results.list$MLEconf.array
  num.reps = results.list$MLE.array.parameters$num.reps
  b.known = results.list$MLE.array.parameters$b.known
  binTypes = results.list$MLE.array.parameters$binTypes
  binType.name = results.list$MLE.array.parameters$binType.name

  #  xbigticks = seq(xrange[1], xrange[2], by = xbigticks.by)
  xsmallticks = seq(xrange[1], xrange[2], by = xsmallticks.by)

  #  yBigTickLab = seq(0, num.reps, yBigTickLab.by)
  #  yBigTickNoLab = seq(0, num.reps, yBigTickNoLab.by)
  #  ySmallTick = seq(0, num.reps, ySmallTick.by)

  par(omi = omi,
      mfrow = mfrow,
      mai = mai,
      xaxs = xaxs,
      yaxs = yaxs,
      mgp = mgp,
      cex = cex)

for(binTypeInd in 1:binTypes)
  {
  # confPlot returns a data.frame of intervals, but no need to save them
  #  for each binType.
    res = confPlot(as.data.frame(MLEconf.array[ ,binTypeInd, "MLEmid", ]),
                   legName = paste(legLabMid[binTypeInd],
                                   binType.name[binTypeInd]),
                   b.true = b.known,
                   xLim = xrange,
                   xsmallticks = xsmallticks,
                   insetVal = insetMat[binTypeInd,],
                   insetVal2 = insetMat[binTypeInd,] + c(0, 0.12),
                   legLoc="topright",
                   yLab="",
                   vertCol = vertCol,
                   vertThick = vertThick)

    res = confPlot(as.data.frame(MLEconf.array[ ,binTypeInd, "MLEbin", ]),
                   legName=paste(legLabBin[binTypeInd],
                                 binType.name[binTypeInd]),
                   b.true = b.known,
                   xLim = xrange,
                   xsmallticks = xsmallticks,
                   insetVal=inset,  # TODO could be wrong
                   insetVal2=inset + c(0, 0.12),
                   yLabels=FALSE,
                   legLoc="topright",
                   yLab="",
                   vertCol = vertCol,
                   vertThick = vertThick)
  }

  mtext(expression(paste("Estimate of ", italic(b))),
        side = 1,
        outer = TRUE,
        line = 0)
  mtext("Sample number",
        side = 2,
        outer = TRUE,
        line = -1)

  mtext("         MLEmid                                              MLEbin",
        side = 3,
        outer = TRUE,
        line = 0)
}


##' Make dataframe of results from MLEmid and MLEbin methods and four binning types
##'
##' Makes dataframe to produce Tables S.3, S.4 and S.5 of MEPS paper, showing the
##' statistics from estimating exponent $b$ for 10,000 simulated data sets, binned using four
##' methods and fitted using the MLEmid and MLEbin methods.
##'
##' @param results.list output list from `MLEbin.simulate()`; see
##'   `?MLEbin.simulate()` for details
##' @return dataframe of table with columns `Binning.type`, `Method`,
##'   `Quantile.5`, `Median, `Mean`, `Quantile.95`, `Percent.below.true`, one
##'   row for each combination of binning type and fitting method.
##' @export
##' @author Andrew Edwards
MLEmid.MLEbin.table = function(results.list)
{
  # Extract required components
  MLE.array = results.list$MLE.array
  num.reps = results.list$MLE.array.parameters$num.reps
  b.known = results.list$MLE.array.parameters$b.known
  binTypes = results.list$MLE.array.parameters$binTypes
  binType.name = results.list$MLE.array.parameters$binType.name

  res.table <- data.frame("Binning type" = NA,
                          "Method" = NA,
                          "Quantile 5" = NA,
                          "Median" = NA,
                          "Mean" = NA,
                          "Quantile 95" = NA,
                          "Percent below true" = NA)
  for(i in 1:binTypes)
  {
    for(j in 1:2)
    {
      res.table[2*i + j - 2, ] = c(binType.name[[i]],
                                   names(MLE.array[1,1,])[j],
                                   qqtab(MLE.array[ ,i,j],
                                         quants=c(0.05, 0.95),
                                         type="data.frame",
                                         true = b.known)
                                   )
    }
  }
  return(res.table)
}

##' Plot time series of estimated *b* with confidence intervals
##'
##' Plot time series of estimated *b* with confidence intervals, for
##' data that are analysed year-by-year by a single method. And then
##' fit a linear regression with its confidence intervals.
##' This will get called eight times to produce a comparison figure of
##' the methods.
##'
##' @param bForYears dataframe with columns `Year`, `Method` (optional), `b`
##'   (the estimate of *b*), `confMin` and `confMax` (the 95\% lower and upper
##'   confidence limits) and `stdError` (the standard error of the estimate of
##'   *b*).
##' @param legName legend name for that panel
##' @param method method used to obtain the inputted estimates of `b`
##' @param weightReg  TRUE if doing weighted regression (using standard errors)
##'   or FALSE to not do weighted regression.
##' @param bCol colour for points for *b*
##' @param pchVal pch for points for *b*
##' @param cexVal size of points for *b*
##' @param confCol colour for confidence intervals for *b*
##' @param confThick thickness of vertical line for confidence intervals
##' @param xLim x-axis limits
##' @param yLim y-axis limits
##' @param xLab label for x-axis
##' @param yLab label for y-axis
##' @param sep TODO: not needed, make the value if yLab is NULL - expect will
##'   fail sometimes
##' @param xTicksSmallInc increments for where to have small (unlabelled)
##'   tickmarks on x-axis
##' @param xTicksSmallTck tick length for small (unlabelled) tickmarks on x-axis
##' @param yLabels whether or not to label main tickmarks on y-axis
##' @param yTicksSmallInc increments for where to have small (unlabelled)
##' tickmarks on y-axis
##' @param yTicksSmallTck tick length for small (unlabelled) tickmarks on y-axis
##' @param legPos legend position
##' @param newPlot TRUE to create a new plot, FALSE to add to existing
##' @param regPlot TRUE to plot the regression line and conf intervals
##' @param regColNotSig colour for regression line (and its confidence intervals)
##' if the trend is not significant
##' @param regColSig colour for regression line (and its confidence intervals)
##' if the trend is significant
##' @param legExtra extra manually-specified legend (e.g. to distinguish two
##' sets of results)
##' @param legExtraPos position for extra manually-specified legend
##' @param legExtraCol colours (vector) for extra manually-specified legend
##' @param insetVal inset shift for naming the panel
##' @param xJitter value to jitter the x-values by (for comparison plot the
##'   confidence intervals overlap)
##' @return dataframe of just one row (TODO: check) with columns:
##'   * Method: method used
##'   * Low: lower bound of 95\% confidence interval
##'   * Trend: gradient of regression fit
##'   * High: upper bound of 95\% confidence interval
##'   * p: p-value of regression fit
##'   * Rsquared: r-squared of regression fit
##'   * adjRsquared: adjusted r-squared of regression fit
##' @export
##' @author Andrew Edwards
timeSerPlot = function(bForYears, legName, method, weightReg = FALSE,
    bCol="black",
    pchVal = 20, cexVal = 1, confCol="black", confThick = 1,
    xLim = NULL, yLim = NULL, xLab="",
    yLab = expression(paste("Estimate of ", italic(b)), sep=""),
    xTicksSmallInc = NULL, xTicksSmallTck = 0.01,
    yLabels = TRUE, yTicksSmallInc = NULL, yTicksSmallTck = 0.01,
    legPos = "topleft", newPlot = TRUE,
    regPlot = TRUE,
    regColNotSig = "darkgrey", regColSig = "red",
    legExtra = NULL, legExtraPos = "topleft", legExtraCol = "",
    insetVal = c(-0.08, -0.06),
    xJitter = 0)
    {
    if(is.null(xLim))
        {
            xLim = range(bForYears$Year)
        }
    if(is.null(yLim))        # just do the yLim for this set of results
        {
        yLim = range(c(bForYears$confMin, bForYears$confMax), na.rm=TRUE)
        }        # For Llin and LT can get only two bins and so NaN for
                 #  conf intervals, so need na.rm.
    if(newPlot)
       {
       plot(bForYears$Year+xJitter, bForYears$b, xlim=xLim, ylim=yLim, col=bCol,
         pch=pchVal, cex=cexVal, xlab=xLab, ylab=yLab)    # yaxt="n")

       legend(legPos, legName, bty="n", inset=insetVal)
       if(!is.null(yTicksSmallInc))
           { yTicksSmall = seq(yLim[1], yLim[-1], by=yTicksSmallInc)
             axis(2, at = yTicksSmall, labels = rep("", length(yTicksSmall)),
                tck=-yTicksSmallTck)
           }
       if(!is.null(xTicksSmallInc))
           { xTicksSmall = seq(xLim[1], xLim[-1], by=xTicksSmallInc)
             axis(1, at = xTicksSmall, labels = rep("", length(xTicksSmall)),
                tck=-xTicksSmallTck)
           }
       # Confidence intervals (instead of plotCI from regress2.Snw):
       segments(x0=bForYears$Year+xJitter, y0=bForYears$confMin,
                x1=bForYears$Year+xJitter, y1=bForYears$confMax)
       if(!is.null(legExtra)) legend(legExtraPos, legExtra, bty="n",
                                     col=legExtraCol, pch=pchVal, cex=cexVal)
       } else    # Add to existing plot
       {
       points(bForYears$Year+xJitter, bForYears$b, col=bCol,
         pch=pchVal, cex=cexVal)
       segments(x0=bForYears$Year+xJitter, y0=bForYears$confMin,
                x1=bForYears$Year+xJitter, y1=bForYears$confMax)
       }


    # Now just fit a linear regression through the points,
    #  and colour code it red if significant trend and grey if not. This
    #  is taken and adapted from regress2.Snw from RBR14 assessment.
    if(weightReg == TRUE)
        { lm = lm(b ~ Year, data = bForYears, weights = 1/(stdErr^2))   } else
        { lm = lm(b ~ Year, data = bForYears) }

    yearInc = seq(xLim[1], xLim[2], 0.1)
    p.conf = predict(lm, newdata=data.frame(Year=yearInc), interval="confidence")
    pVal = summary(lm)$coeff["Year",4]
    if(regPlot)
      {
        if (pVal > 0.05) regCol = regColNotSig else regCol= regColSig
        lm.line(xLim, lm, col=regCol)
        matlines(yearInc,
                 p.conf[ ,c("lwr", "upr")],
                 col=regCol,
                 lty=2)
      }
    confVals = confint(lm, "Year", level=0.95)

    res = data.frame(Method = method,
                     Low = confVals[1],
                     Trend = lm$coeff[2],
                     High = confVals[2],
                     p = pVal,
                     Rsquared = summary(lm)$r.squared,
                     adjRsquared = summary(lm)$adj.r.squared,
                     row.names=NULL)
    return(res)
}

##' Recommended plots of individual size distribution and fit for binned data
##'
##' Plots the individual size distribution and the fit from the MLEbins method,
##' with linear and then logarithmic y-axes, in the recommended way that takes
##' into account the bin structure of the data, as in Figures 7 and S.5-S.34 of
##' the MEPS paper.
##'
##' @param data.year tibble containing columns Year, wmin, wmax, Number,
##'   countGTEwmin, lowCount, highCount
##' @param b.MLE maximum likelihood estimate of *b* (ideally from the MLEbins method)
##' @param b.confMin lower 95\% confidence limits of *b*
##' @param b.confMax upper 95\% confidence limits of *b*
##' @param year year of data to go into legend (use NA if not applicable)
##' @param xRange limits of x-axis
##' @param MLE.round number of decimal places to show MLE of b on the top plot
##' @param xLabel.small which small tickmarks to label on the x-axis
##' @param yBig.inc increment for big tickmarks on the y-axis
##' @param ySmall.inc increment for small unlabelled tickmarks on the y-axis
##' @param ySmall.tcl length of small y-axis tick marks - only for (a)??TODO
##' @param xLab label for x-axis
##' @param yLab label for y-axis
##' @param inset.a how far to inset (a) and (b)
##' @param inset.year how far to inset the year
##' @param seg.col colour for horizontal green lines showing range of body sizes
##'   for each bin
##' @param rect.col colour to fill in the rectangles for each bin
##' @param fit.col colour to plot the MLE curve, and those for the confidence intervals
##' @param fit.lwd line weight for MLE curve
##' @param conf.lty line type for two curves for the MLE confidence intervals
##' @param par.mfrow vector giving the layout of the figures (number of rows,
##'   number of columns)
##' @param par.mai margin size in inches
##' @param par.cex magnification of plotting text and symbols relative to default
##' @param mgp.vals margin line for axis title, axis labels and axis line
##' @return two-panel figure of the recommended plot of binned data and the
##'   fitted individual size distribution, like Figures 7 and S.5-S.34 of the
##'   MEPS paper. See the vignette TODO for explicit example.
##' @export
##' @author Andrew Edwards
ISD_bin_plot <- function(data.year,
                         b.MLE,
                         b.confMin,
                         b.confMax,
                         year = NA,
                         xRange = NA,
                         MLE.round = 2,
                         xLabel.small = c(5, 50, 500, 5000),
                         yBig.inc = 1000,
                         ySmall.inc = 250,
                         ySmall.tcl = -0.2,
                         xLab = expression(paste("Body mass (", italic(x), "), g")),
                         yLab = expression(paste("Number of ", values >= italic(x))),
                         inset.a = c(0, 0),
                         inset.year = c(0, 0.04),
                         seg.col = "green",
                         rect.col = "grey",
                         fit.col = "red",
                         fit.lwd = 2,
                         conf.lty = 2,
                         par.mfrow = c(2, 1),
                         par.mai = c(0.4, 0.5, 0.05, 0.3),
                         par.cex = 0.7,
                         mgp.vals = c(1.6,0.5,0)
                         )
  {
  sumNumber = sum(data.year$Number)

  par(mfrow = par.mfrow)
  par(mai = par.mai, cex = par.cex)  # Affects all figures   TODO put into one

  if(is.na(xRange[1])){
    xRange = c(min(data.year$wmin),
               max(data.year$wmax))  # For PLB line
    }

  x.PLB = seq(xRange[1],
              xRange[2],
              length=10000)
          # x values to plot PLB, need high resolution for both plots, but
          #  insert value close to xmax to make log-log curve go down further
  x.PLB.length = length(x.PLB)
  x.PLB = c(x.PLB[-x.PLB.length],
            0.99999 * x.PLB[x.PLB.length],
            x.PLB[x.PLB.length])
  y.PLB = (1 - pPLB(x = x.PLB,
                    b = b.MLE,
                    xmin = min(x.PLB),
                    xmax = max(x.PLB))) * sumNumber
  # To add curves for the limits of the 95% confidence interval of b:
  y.PLB.confMin = (1 - pPLB(x = x.PLB,
                   b = b.confMin,
                   xmin = min(x.PLB),
                   xmax = max(x.PLB))) * sumNumber
  y.PLB.confMax = (1 - pPLB(x = x.PLB,
                   b = b.confMax,
                   xmin = min(x.PLB),
                   xmax = max(x.PLB))) * sumNumber

  # yRange = c(min(data.year$lowCount), max(data.year$highCount))
  # The above does not work because first val is 0 which is not permissable on
  #  log axis. Which also means that the rectangle that goes to 0 has to be
  #  added manually (below). Picking the y-axis to go down to 0.75 of the
  #  minimum value of CountGTEwmin.
  yRange = c(0.75*min(data.year$countGTEwmin), max(data.year$highCount))   # TODO
                                        # still don't quite get this

  # y-axis not logged
  plot(data.year$wmin,
       data.year$countGTEwmin,
       log="x",
       xlab=xLab,
       ylab=yLab,
       xlim = xRange,
       ylim = yRange,
       type = "n",
       axes = FALSE,
       mgp = mgp.vals)

  xLim = 10^par("usr")[1:2]
  yLim = 10^par("usr")[3:4]

  logTicks(xLim,
           yLim=NULL,
           xLabelSmall = xLabel.small)
  yBig = seq(0, yRange[2], yBig.inc)  # May have to tweak for some years
  # Big labelled:
  axis(2, at = yBig, labels = yBig, mgp=mgp.vals)
  # Small unlabelled:
  axis(2,
       seq(yBig[1],
           yRange[2]*1.1,
           by = ySmall.inc),
       labels = rep("",
                    length(seq(yBig[1],
                               yRange[2]*1.1,
                               by=ySmall.inc))),
       tcl = ySmall.tcl,
       mgp = mgp.vals)

  rect(xleft = data.year$wmin,
       ybottom = data.year$lowCount,
       xright = data.year$wmax,
       ytop = data.year$highCount,
       col = rect.col)
  segments(x0 = data.year$wmin,
           y0 = data.year$countGTEwmin,
           x1 = data.year$wmax,
           y1 = data.year$countGTEwmin,
           col = seg.col)
  lines(x.PLB, y.PLB, col = fit.col, lwd = fit.lwd)   # Plot line last so can see it
  lines(x.PLB, y.PLB.confMin, col = fit.col, lty = conf.lty)
  lines(x.PLB, y.PLB.confMax, col = fit.col, lty = conf.lty)

  legend("topright", "(a)",
         bty = "n",
         inset = inset.a)
  if(!is.na(year)){
    legend("topright",
           legend = year,
           bty = "n",
           inset = inset.year)
  }

  legend("topright",
         legend = paste0("b=", round(b.MLE, MLE.round)),
         bty = "n",
         inset = 2 * inset.year)

  legend("topright",
         legend = paste0("n=", round(yRange[2], 2)),
         bty = "n",
         inset = 3 * inset.year)

  # y-axis logged
  # empty plot:
  plot(data.year$wmin,
       data.year$countGTEwmin,
       log = "xy",
       xlab = xLab,
       ylab = yLab,
       xlim = xRange,
       ylim = yRange,
       type = "n",
       axes = FALSE,
       mgp = mgp.vals)

  xLim = 10^par("usr")[1:2]
  yLim = 10^par("usr")[3:4]

  logTicks(xLim,
           yLim,
           xLabelSmall = xLabel.small)

  rect(xleft = data.year$wmin,
       ybottom = data.year$lowCount,
       xright = data.year$wmax,
       ytop = data.year$highCount,
       col = rect.col)

  # Need to manually draw the rectangle with lowCount = 0 since it doesn't
  #  get plotted on log-log plot
  extra.rect = dplyr::filter(data.year,
                             lowCount == 0)
  # if(nrow(extra.rect) > 1) stop("Check rows of extra rect.")
  rect(xleft = extra.rect$wmin,
       ybottom = rep(0.01 * yRange[1], nrow(extra.rect)),
       xright = extra.rect$wmax,
       ytop = extra.rect$highCount,
       col = rect.col)

  segments(x0 = data.year$wmin,
           y0 = data.year$countGTEwmin,
           x1 = data.year$wmax,
           y1 = data.year$countGTEwmin,
           col = seg.col)

  lines(x.PLB,
        y.PLB,
        col = fit.col,
        lwd = fit.lwd)
  lines(x.PLB,
        y.PLB.confMin,
        col = fit.col,
        lty = conf.lty)
  lines(x.PLB,
        y.PLB.confMax,
        col = fit.col,
        lty = conf.lty)
  legend("topright",
         "(b)",
         bty="n",
         inset = inset.a)
}
