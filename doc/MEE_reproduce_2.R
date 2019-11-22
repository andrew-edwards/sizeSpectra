## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----setup---------------------------------------------------------------
library(sizeSpectra)

## ------------------------------------------------------------------------
library(sizeSpectra)

## ---- fig.height=7-------------------------------------------------------
list2env(eight.results.default, envir=.GlobalEnv)    # Extract previously saved results


# Plotting
# brange = range( c( Llin.rep, LT.rep, LTplus1.rep, LBmiz.rep,
#    LBbiom.rep, LBNbiom.rep, MLE.rep))
# Comes out as -3.38, 0.233 for seed=42, num.reps = 10000
xrange = c(-3.5, 0.5)         # range of x-axis for histograms - actually to define bins
xbigticks = seq(-3.5, 0.5, 0.5) # NOW THE DEFAULT

xsmallticks = seq(xrange[1], xrange[2], by=0.1)  # ALSO now default
breakshist = seq(xrange[1], xrange[2], by=3/61)
         # bin ends at 0, 2 is centered. width of 4/81 from solving
         #  the mean of a bin equals 2, ((N+1)w + Nw)/2 = 2, and setting number
         #  of the bin that includes 2, N, =40 (which is how many you get with
         #  0.05, which looked good).

# figheight = 7 # 5.6
# figwidth = 5.7    # 5.7 inches for JAE
#postscript("fitting3rep.eps", height = figheight, width = figwidth,
#           horizontal=FALSE,  paper="special")
par(omi = c(0.12, 0.05, 0.12, 0.0))      # outer margins in inches
par(mfrow=c(4,2)) #7,1))

oldmai = par("mai")
par(mai=c(0.5, 0.5, 0.1, 0.3)) # Affects all figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small

vertCol = "red"            # Colour for vertical lines
vertThick = 1              # Thickness for vertical lines

xLim = c(-3.2, -1.5)
xLimLlin = xLim + 2
ylimA = c(0, 5500)
ylimLlin = c(0, 10000)
cexAxis = 0.9      # font size for axes labels to make the y ones fit okay

# Llin.rep has different breakhist so that 0 is a breakpoint. Same width as
#  others, just shifted.
hist(Llin.rep.df$slope, xlim=xLimLlin, breaks=breakshist - breakshist[31],
  xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimLlin)  #  ylim=c(0,1100),
axis(1, at=xbigticks, labels = xbigticks, mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks
axis(1, at=xsmallticks, labels=rep("",length(xsmallticks)), tcl=-0.2)
axis(2, at=c(0, 5000, 10000),
         labels = c(0, 5000, 10000),
         mgp=c(1.7,0.7,0), cex.axis=cexAxis)  # big ticks labelled
axis(2, at=seq(0, 10000, 1000),
         labels = rep("", 11),
         mgp=c(1.7,0.7,0))  # big ticks unlabelled
abline(v=b.known, col=vertCol, lwd=vertThick)
inset = c(-0.08, -0.08)     # inset distance of legend
legend("topleft", "(a) Llin", bty="n", inset=inset)

# TODO think these three lines can be deleted
figlabpos = 0.93    # proportion of x and y axis lengths to put (a) in etc.
xlabpos = 0.75       # just play with a number, as now using pos=4 in text
text( xlabpos, figlabpos * 1100, "(a) Llin", pos=4)

hist(LT.rep.df$slope, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(b) LT", bty="n", inset=inset)

hist(LTplus1.rep.df$slope, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(c) LTplus1", bty="n", inset=inset)

hist(LBmiz.rep.df$slope - 1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(d) LBmiz", bty="n", inset=inset)

hist(LBbiom.rep.df$slope - 2, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(e) LBbiom", bty="n", inset=inset)

hist(LBNbiom.rep.df$slope - 1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(f) LBNbiom", bty="n", inset=inset)


hist(LCD.rep.df$slope - 1, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(g) LCD", bty="n", inset=inset)

hist(MLE.rep.df$b, xlim=xLim, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = ylimA)  #  ylim=c(0,1100),
histAxes()
legend("topleft", "(h) MLE", bty="n", inset=inset)

# text(-2.5, 0, expression(paste("Estimate of ", italic(b))))
mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=-1)
# mtext("hello", side=1, outer=TRUE, line=-1)

## ------------------------------------------------------------------------

par(omi = c(0.12, 0.05, 0.12, 0.0))      # outer margins in inches
par(mfrow=c(2,1))

par(mai=c(0.5, 0.5, 0.1, 0.3))
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small

xLimFix = c(-2.5, -1.5)
yLimFix = c(0, 6000)
hist(MLE.rep.df$b, xlim=xLimFix, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = yLimFix)  #  ylim=c(0,1100),
histAxes()
axis(2, at=6000, labels = 6000, mgp=c(1.7,0.7,0), cex.axis=cexAxis)
                                        # big tick label
legend("topleft", "(a) MLE", bty="n", inset=0.5*inset)

hist(MLEfix.rep.df$b, xlim=xLimFix, breaks=breakshist, xlab="", ylab="Frequency",
  main="", axes=FALSE, ylim = yLimFix)  #  ylim=c(0,1100),
histAxes()
axis(2, at=6000, labels = 6000, mgp=c(1.7,0.7,0), cex.axis=cexAxis)
                                        # big tick label
legend("topleft", "(b) MLEfix", bty="n", inset=0.5*inset)

mtext(expression(paste("Estimate of ", italic(b))), side=1, outer=TRUE, line=-1)

## ---- fig.height=7-------------------------------------------------------
par(omi = c(0.14, 0, 0.1, 0.15))      # outer margins in inches
par(mfrow=c(4,2)) #7,1))

oldmai = par("mai")    #   0.6732 0.5412 0.5412 0.2772  inches I think,
                       #    may be indpt of fig size
par(mai=c(0.3, 0.5, 0.08, 0))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small
vertThick = 1              # Thickness for vertical lines

# Each of these plots a panel for one method. Define xLim if the default
#  (integer-based calculation) is not suitable
Llin.rep.conf.sort = confPlot(Llin.rep.df[ ,-1], legName="(a) Llin",
    xLim = c(-0.25, 0.25)) #, colourCode=FALSE)

LT.rep.conf.sort = confPlot(LT.rep.df[ ,-1], legName="(b) LT", yLab="", yLabels=FALSE)

LTplus1.rep.conf.sort = confPlot(LTplus1.rep.df[ ,-1], legName="(c) LTplus1")

xLimCom = c(-2.65, -1.5)     # common width of axes for the remainder

# range(LBmiz.rep.conf) =  -1.5894281 -0.5268428. Was -1.6, -0.4
LBmiz.rep.conf.sort = confPlot(LBmiz.rep.df[ ,-1] - 1, legName="(d) LBmiz",
    xLim = xLimCom, yLab="", yLabels=FALSE)

# range(LBbiom.rep.conf) = -0.6106510  0.4370617
LBbiom.rep.conf.sort = confPlot(LBbiom.rep.df[ ,-1] - 2, legName="(e) LBbiom",
      xLim = xLimCom)

# range(LBNbiom.rep.conf) = -1.6106510 -0.5629383
LBNbiom.rep.conf.sort = confPlot(LBNbiom.rep.df[ ,-1] - 1, legName="(f) LBNbiom",
     yLab="", yLabels=FALSE, xLim = xLimCom)

# range(LCD.rep.conf) = -1.1744680 -0.8713301   was c(-1.2, -0.8)
LCD.rep.conf.sort = confPlot(LCD.rep.df[ ,-1] - 1, legName="(g) LCD",
    xLim = xLimCom, vertFirst = TRUE)

# c(-2.2, -1.8) - looks good, but better to have consistent
MLE.rep.conf.sort = confPlot(MLE.rep.df[ ,-1], legName="(h) MLE",
    xLim = xLimCom, yLab="", yLabels=FALSE)
mtext(expression(paste("Estimate of ", italic(b)), sep=""),
      side=1, outer=TRUE, line=-0.2, cex=0.8)

## ------------------------------------------------------------------------
# Confidence intervals for MLE and MLEfix methods.
#postscript("fitting3confMLEfix.eps", height = 0.8*figheight,
#           width = 0.8*figwidth,
#           horizontal=FALSE,  paper="special")
par(omi = c(0.12, 0.05, 0.12, 0.0))      # outer margins in inches
par(mfrow=c(2,1)) #7,1))

oldmai = par("mai")    #   0.6732 0.5412 0.5412 0.2772  inches I think,
                       #    may be indpt of fig size
par(mai=c(0.5, 0.5, 0.1, 0.3))  # Affects all four figures if don't change agaiin
par(xaxs="i", yaxs="i")    # Have to define here for hist
par(mgp=c(2.0, 0.5, 0))    # puts axes labels closer I think
par(cex = 0.8)             # With no option all text comes out a bit small

# xrangeMLEfix = c(-2.225, -1.775)  # Range for these two figures
# breakshistMLEfix = seq(xrangeMLEfix[1], xrangeMLEfix[2], by=3/61)

MLE.rep.conf.sort = confPlot(MLE.rep.df[ ,-1], legName="(a) MLE",
    xLim = xLimCom, insetVal = c(-0.04, -0.03), insetVal2 = c(-0.04, 0.05))


MLEfix.rep.conf.sort = confPlot(MLEfix.rep.df[ ,-1], legName="(b) MLEfix",
    xLim = xLimCom, insetVal = c(-0.04, -0.03), insetVal2 = c(-0.04, 0.05))

mtext(expression(paste("Estimate of ", italic(b)), sep=""),
      side=1, outer=TRUE, line=-0.2, cex=0.8)


