##' Example length-weight relationships for two species, demonstrating
##'  consequences of length bins
##'
##' Reproduces Figure 2 of MEPS paper (for default values), showing how the same
##' length bins translate to different body-mass bins for different species. Has
##' not been fully tested with other input values (and some axis values etc are
##' still hardwired).
##'
##' @return plots Figure, default is Figure 2 of MEPS
##' @export
##' @author Andrew Edwards
##' @param sp1 Species 1 name (for legend)
##' @param sp2 Species 2 name (for legend)
##' @param LWa Vector of two alpha values (see MEPS equation 1), for species 1
##'   then species 2
##' @param LWb Vector of two beta values (see MEPS equation 1), for species 1
##'   then species 2
##' @param col1 Vector of two colours to alternate for species 1
##' @param col2 Vector of two colours to alternate for species 2
##' @param thick Thickness of bin segments
##' @param curve_thick Thickness for curves
##' @param length Sequence of length values to plot the length-weight
##'   relationships
##' @param lenBins Example length bin breaks to show as examples
##' @param xlim,ylim Limits for x and y axes
##' @param xlab,ylab Labes for x and y axes
##' @param inset Inset for the legend
##' @param xaxs,yaxs,mgp,lend standard options for `par()`,
##'   defaults are for Figure 2 of MEPS
length_weight_plot <- function(sp1 = "Common Ling",
                               sp2 = "Lemon Sole",
                               LWa = c(0.001, 0.0255),
                               LWb = c(3.4362, 2.7643),
                               col1 = c("red", "pink"),
                               col2 = c("blue", "lightblue"),
                               thick = 7,
                               curve_thick = 3,
                               length = 10:50,
                               lenBins = seq(10, 40, by=5),
                               xlim=c(-6, 50),
                               ylim=c(-60, 800),
                               inset = c(0.1,-0.02),
                               xlab="Length, cm",
                               ylab="Body mass, g",
                               xaxs = "i",
                               yaxs = "i",
                               mgp = c(2.0, 0.5, 0),
                               lend = "butt")
{
  par(xaxs = xaxs,
      yaxs = yaxs,
      mgp = mgp,
      lend = lend)
  # Species specific body masses corresponding to each value of length
  mass = rbind(lengthToMass(length, LWa[1], LWb[1]),
               lengthToMass(length, LWa[2], LWb[2]))
  plot(length,
       mass[1,],
       col = col1[1],
       type = "l",
       xlab = xlab,
       ylab = ylab,
       xlim = xlim,
       ylim = ylim,
       lwd=curve_thick)
  lines(length,
        mass[2,],
        col=col2[1],
        lwd=curve_thick)

  # Switching sp1 and sp2 around in the legend (and thus the text), since
  #  looks more consistent:
  legend("topleft",
         legend = c(sp2, sp1),
         lty = 1,
         lwd = curve_thick,
         col = c(col2[1], col1[1]),
         bty = "n",
         inset = inset)
  axis(1,
       at = seq(0, 50, by=5),
       labels = rep("", 11),
       tck=-0.02)
  axis(2,
       at = seq(0, 800, by=100),
       labels=rep("", 9),
       tck=-0.02)
  axis(2,
       at = seq(0, 800, by=50),
       labels=rep("", 17),
       tck=-0.01)

  # Example bins:
  numBinBreaks = length(lenBins)
  lenBinsStart = lenBins[-numBinBreaks]     # starting point for length bins
  lenBinsEnd = lenBins[-1]
  # resulting species-specific bin breaks for body mass:
  massBins = rbind(lengthToMass(lenBins,
                                LWa[1],
                                LWb[1]),
                   lengthToMass(lenBins,
                                LWa[2],
                                LWb[2]))
  massBinsStart = massBins[, -numBinBreaks]
  massBinsEnd = massBins[, -1]

  # colours for bins for each species:
  colBins = rbind(rep(col1,
                      length=round(length(lenBinsStart))),
                  rep(col2,
                      length=round(length(lenBinsStart))))
  yLengths = c(-40, -20)  # y values to plot length bins, for each species
  xMasses  = c(-3, -1)    # x values to plot mass bins

  # length bins:
  segments(x0 = lenBinsStart,
           y0 = yLengths[1],
           x1 = lenBinsEnd,
           col = colBins[1,],
           lwd = thick)  # y1=y0
  segments(x0 = lenBinsStart,
           y0 = yLengths[2],
           x1 = lenBinsEnd,
           col = colBins[2,],
           lwd = thick)

  # resulting mass bins:
  segments(x0 = xMasses[1],
           y0 = lengthToMass(lenBinsStart, LWa[1], LWb[1]),
           y1 = lengthToMass(lenBinsEnd, LWa[1], LWb[1]),
           col = colBins[1,],
           lwd = thick)
  segments(x0 = xMasses[2],
           y0 = lengthToMass(lenBinsStart, LWa[2], LWb[2]),
           y1 = lengthToMass(lenBinsEnd, LWa[2], LWb[2]),
           col = colBins[2,],
           lwd = thick)
  egBin = 5         # example bin number to highlight
  offset = 0.1      # offset to shift species vertical lines so they show up
  midOfBin = mean(lenBinsStart[c(egBin, egBin+1)])    # midpoint of example bin
  massMidOfBin = rbind(lengthToMass(midOfBin, LWa[1], LWb[1]),
                       lengthToMass(midOfBin, LWa[2], LWb[2]))

  # example bin for species 1:
  lines(c(rep(lenBinsStart[egBin],
              2),
          xMasses[1]) - offset,
        c(yLengths[1], rep(lengthToMass(lenBinsStart[egBin], LWa[1], LWb[1]),
                           2)),
        lty = 3,
        lwd = 1,
        col = col1[1])
  lines(c(rep(lenBinsStart[egBin+1],
              2),
          xMasses[1]) - offset,
        c(yLengths[1],
          rep(lengthToMass(lenBinsStart[egBin+1], LWa[1], LWb[1]),
              2)),
        lty = 3,
        lwd = 1,
        col = col1[1])
  lines(c(rep(midOfBin,
              2),
          xMasses[1]) - offset,
        c(yLengths[1],
          rep(massMidOfBin[1,],
              2)),
        lty = 1,
        lwd = 1,
        col = col1[1])

  # example bin for species 2:
  lines(c(rep(lenBinsStart[egBin],
              2),
          xMasses[2]) + offset,
        c(yLengths[2],
          rep(lengthToMass(lenBinsStart[egBin],
                           LWa[2],
                           LWb[2]),
              2)),
        lty = 3,
        lwd = 1,
        col = col2[1])
  lines(c(rep(lenBinsStart[egBin+1],
              2),
          xMasses[2]) + offset,
        c(yLengths[2],
          rep(lengthToMass(lenBinsStart[egBin+1], LWa[2], LWb[2]),
              2)),
        lty = 3,
        lwd = 1,
        col = col2[1])
  lines(c(rep(midOfBin, 2),
          xMasses[2]) + offset,
        c(yLengths[2],
          rep(massMidOfBin[2,],
              2)),
        lty = 1,
        lwd = 1,
        col = col2[1])
}

##' Demonstration of how binned body-mass values get assigned to logarithmic
##'  size-class bins (Figure 3 of MEPS)
##'
##' @param num2 Number of log2 bin breaks
##' @param lenBins Example length bin breaks to use for the example body-mass bins
##' @param LWa Vector of two alpha values (see MEPS equation 1), for species 1
##'   then species 2, though only species 2 is used, just keeping consistency
##'   with `length_weight_plot()`
##' @param LWb Vector of two beta values (see MEPS equation 1), for species 1
##'   then species 2, though only species 2 is used, just keeping consistency
##'   with `length_weight_plot()`
##' @param thick Thickness of bin segments
##' @param col1 Vector of two colours to alternate for species 1
##' @param col2 Vector of two colours to alternate for species 2
##' @param xaxs,yaxs,mgp,lend standard options for `par()`, defaults are for
##'   Figure 3 of MEPS
##' @return Plots figure, default is Figure 3 of MEPS
##' @export
##' @author Andrew Edwards
bins_assignment_plot <- function(num2 = 11,
                                 lenBins = seq(10, 40, by=5),
                                 LWa = c(0.001, 0.0255),
                                 LWb = c(3.4362, 2.7643),
                                 col1 = c("red", "pink"),
                                 col2 = c("blue", "lightblue"),
                                 thick = 7,
                                 xaxs = "i",
                                 yaxs = "i",
                                 mgp = c(2.0, 0.5, 0),
                                 lend = "butt")
{
  par(xaxs = xaxs,
      yaxs = yaxs,
      mgp = mgp,
      lend = lend)

  log2bins = 2^(0:(num2-1))    # on unlogged axes
  log2binsStart = log2bins[-num2]
  log2binsEnd = log2bins[-1]   # see segments in next figure
  log2binsMid = colMeans(rbind(log2binsStart,
                               log2binsEnd))  # midpoints (unlogged)

  # Axes ranges not automated - should do if want to start changing them.
  plot(100,
       100,
       xlab="",
       ylab="Body mass, g",
       xlim=c(0, 12),
       ylim=c(-60, 800),
       xaxt="n")      # dummy points to set up axes.

  axis(2,
       at = seq(0, 800, by=100),
       labels=rep("", 9),
       tck=-0.02)
  axis(2,
       at = seq(0, 800, by=50),
       labels=rep("", 17),
       tck=-0.01)
  axis(1,
       at = seq(1, 11, by=2),
       labels=c("Bin 1", "Bin 2", "Bin 3", "Bin 4",
                "Bin 5", "Bin 6"),
       tck=0,
       col.axis="blue")

  numBinBreaks = length(lenBins)
  lenBinsStart = lenBins[-numBinBreaks]     # starting point for length bins
  lenBinsEnd = lenBins[-1]

  # resulting species-specific bin breaks for body mass:
  massBins = rbind(lengthToMass(lenBins,
                                LWa[1],
                                LWb[1]),
                   lengthToMass(lenBins,
                                LWa[2],
                                LWb[2]))
  massBinsStart = massBins[, -numBinBreaks]
  massBinsEnd = massBins[, -1]

  colBins = rbind(rep(col1,
                      length = round(length(lenBinsStart))),
                  rep(col2,
                      length = round(length(lenBinsStart))))

  # Could go back and use this throughout, as clearer (? - think for
  #  length_weight_plot() also?). Do for one species here
  #  then expand to do for two species (by having columns spName, LWa, LWb).
  binsSp2 = data.frame(lenStart = lenBinsStart,
                       lenEnd = lenBinsEnd)
  binsSp2 = dplyr::mutate(binsSp2,
                          lenMid = (lenStart +lenEnd)/2)
  binsSp2 = dplyr::mutate(binsSp2,
                          massStart = lengthToMass(lenStart,
                                                   LWa[2],
                                                   LWb[2]),
                          massEnd = lengthToMass(lenEnd,
                                                 LWa[2],
                                                 LWb[2]),
                          massMid = (massStart + massEnd) / 2,
                          massOfLenMid = lengthToMass(lenMid, LWa[2], LWb[2]))
        # massOfLenMid is the mass corresponding to the converted midpoint
        #  of length bins, which is actually what gets used in MLEmid

  # Now assign each original length bin a log2 mass bin:
  ind=vector()             # index of which log2 bin each mass bin ends up in
  for(ii in 1:dim(binsSp2)[1])
    {
      ind[ii] = which(binsSp2[ii, "massOfLenMid"] >= log2binsStart &
                      binsSp2[ii,"massOfLenMid"] < log2binsEnd)
  }
  binsSp2 = cbind(binsSp2,
                  log2binStart = log2binsStart[ind],
                  log2binEnd = log2binsEnd[ind],
                  log2binMid = log2binsMid[ind])

  # Show where example bin ends up

  # 6 bins, therefore contain each pair in a span of 2 wide.
  xVals = 1.5 + seq(0, 10, 2)     # xvals for vertical mass bins
  for(egBin2 in 1:6)              # example bin number to highlight here
    {
      xVal = xVals[egBin2]        # where to have vertical bars
      turn= xVal - 0.5            # where to turn the line
      end = xVal - 0.8            # where to end the line
      # log2 bins:
      xLog2 = end - 0.2             # where to place log2 bins
      abline(v = seq(0, 12, 2),
             col = "grey")
      segments(x0 = rep(xLog2, num2),
               y0 = log2binsStart,
               y1 = log2binsEnd,
               col = c("black", "grey"),
               lwd = thick)
      points(rep(xLog2, num2-1),
             log2binsMid,
             pch = "-",
             cex = 2,
             col = "purple")

    segments(x0 = xVal,
             y0 = massBinsStart[2,],
             y1 = massBinsEnd[2,],
             col = colBins[2,],
             lwd = thick)
    lines(c(xVal, turn, end),
          c(rep(binsSp2[egBin2, "massOfLenMid"],
                2),
            binsSp2[egBin2, "log2binMid"]),
          lty = 1,
          lwd = 1,
          col = col2[1])
    lines(c(xVal, turn, end),
          c(rep(binsSp2[egBin2, "massStart"],
                2),
            binsSp2[egBin2, "log2binMid"]),
          lty = 2,
          lwd = 1,
          col = col2[1])
    lines(c(xVal, turn, end),
          c(rep(binsSp2[egBin2, "massEnd"],
                2),
            binsSp2[egBin2, "log2binMid"]),
          lty = 2,
          lwd = 1,
          col = col2[1])
    # Arrow for the final part:
    shape::Arrows(x0 = turn,
                  y0 = binsSp2[egBin2, "massOfLenMid"],
                  x1  =  end,
                  y1  =  binsSp2[egBin2, "log2binMid"],
                  arr.adj = 1,
                  col  =  col2[1],
                  lwd = 1,
                  arr.lwd = 0.1,
                  arr.type = "triangle")
   }
}
