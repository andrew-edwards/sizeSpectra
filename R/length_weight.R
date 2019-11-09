
##' Example length-weight relationships for two species, demonstrating
##'  consequences of length bins
##'
##' Reproduces Figure 2 of MEPS paper (for default values), showing how the same
##' length bins translate to different body-mass bins for different species. Has
##' not been fully tested with other input values (and some axis values etc are
##' still hardwired).
##'
##' @return plots figure
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

length_weight <- function(sp1 = "Common Ling",
                          sp2 = "Lemon Sole",
                          LWa = c(0.001, 0.0255),
                          LWb = c(3.4362, 2.7643),
                          col1 = c("red", "pink"),
                          col2 = c("blue", "lightblue"),
                          thick=7,
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
         legend=c(sp2, sp1),
         lty=1,
         lwd=curve_thick,
         col=c(col2[1], col1[1]),
         bty="n",
         inset = inset)
  axis(1, at = seq(0, 50, by=5), labels=rep("", 11), tck=-0.02)
  axis(2, at = seq(0, 800, by=100), labels=rep("", 9), tck=-0.02)
  axis(2, at = seq(0, 800, by=50), labels=rep("", 17), tck=-0.01)

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
  segments(x0 = xMasses[1], y0 = lengthToMass(lenBinsStart, LWa[1], LWb[1]),
         y1 = lengthToMass(lenBinsEnd, LWa[1], LWb[1]),
         col = colBins[1,], lwd = thick)
segments(x0 = xMasses[2], y0 = lengthToMass(lenBinsStart, LWa[2], LWb[2]),
         y1 = lengthToMass(lenBinsEnd, LWa[2], LWb[2]),
         col = colBins[2,], lwd = thick)

egBin = 5         # example bin number to highlight
offset = 0.1      # offset to shift species vertical lines so they show up
midOfBin = mean(lenBinsStart[c(egBin, egBin+1)])    # midpoint of example bin
massMidOfBin = rbind(lengthToMass(midOfBin, LWa[1], LWb[1]),
    lengthToMass(midOfBin, LWa[2], LWb[2]))
# example bin for species 1:
lines(c(rep(lenBinsStart[egBin], 2), xMasses[1])-offset,
       c(yLengths[1], rep(lengthToMass(lenBinsStart[egBin], LWa[1], LWb[1]), 2)),
       lty = 3, lwd = 1, col = col1[1])
lines(c(rep(lenBinsStart[egBin+1], 2), xMasses[1])-offset,
       c(yLengths[1], rep(lengthToMass(lenBinsStart[egBin+1], LWa[1], LWb[1]),
         2)), lty = 3, lwd = 1, col = col1[1])
lines(c(rep(midOfBin, 2), xMasses[1]) - offset,
       c(yLengths[1], rep(massMidOfBin[1,], 2)),
       lty = 1, lwd = 1, col = col1[1])

# example bin for species 2:
lines(c(rep(lenBinsStart[egBin], 2), xMasses[2]) + offset,
       c(yLengths[2], rep(lengthToMass(lenBinsStart[egBin], LWa[2], LWb[2]), 2)),
       lty = 3, lwd = 1, col = col2[1])
lines(c(rep(lenBinsStart[egBin+1], 2), xMasses[2]) + offset,
       c(yLengths[2], rep(lengthToMass(lenBinsStart[egBin+1], LWa[2], LWb[2]),
         2)), lty = 3, lwd = 1, col = col2[1])
lines(c(rep(midOfBin, 2), xMasses[2]) + offset,
       c(yLengths[2], rep(massMidOfBin[2,], 2)),
       lty = 1, lwd = 1, col = col2[1])
}
