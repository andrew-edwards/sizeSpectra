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

##' Produce LaTeX or Rmarkdown table code for quantiles of results
##'
##'  Quantile table function, to construct lines of quantiles to copy
##'   into LaTeX, or render directly in an Rmarkdown document. Does `quants[1]`,
##'   50\% (median), mean, `quants[2]` quantiles,
##'   and the \%age of values < the true value.
##'
##' @param xx vector of values to give quantiles for
##' @param dig number of decimal places to give
##' @param true the true value of the quantity being estimated, will
##'     depend on method
##' @param quants quantiles to use, with defaults of 0.25 and 0.75 (25\% and
##'   75\%).
##' @param latex if TRUE gives `latex` output to copy into a `.tex` file, else gives `Rmd` gives Rmarkdown output
##'   (see vignette TODO).
##' @return LaTeX output for table to copy into a .tex file
##'
##' @export
##' @author Andrew Edwards
qqtab = function(xx,
                 dig=2,
                 true=b.known,
                 quants = c(0.25, 0.75),
                 latex = TRUE)
   { if(latex){
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
  } else {
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
    legend(legLoc, paste(round(sum.inConf*100, dig=0), "%", sep=""),
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
##' @return Adds axes to existing histogram TODO: check
##' @export
##' @author Andrew Edwards
histAxes = function(yBigTickLab = c(0, 2000, 4000),
                    yBigTickNoLab = seq(0, 6500, 1000),
                    ySmallTick = seq(0, 6500, 500))
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
