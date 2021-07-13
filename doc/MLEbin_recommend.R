## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5.7,
  fig.height = 7
)

## ----setup--------------------------------------------------------------------
library(sizeSpectra)
# library(tibble)  # Else prints all of a tibble

## ----generate-----------------------------------------------------------------
set.seed(42)
x <- rPLB(n = 1000,
          b = -2,
          xmin = 1,
          xmax = 1000)

x.binned <- binData(x,
                    binWidth = "2k")  # 1, 2, 5 or "2k"

head(x.binned$indiv)     # Individual x values and which bin they are assigned to
x.binned$binVals         # Bins and the counts in each bin (with extra columns
                         #  that we won't need here).

## ----likelihood---------------------------------------------------------------
num.bins <- nrow(x.binned$binVals)

# bin breaks are the minima plus the max of the final bin:
binBreaks <- c(dplyr::pull(x.binned$binVals, binMin),
               dplyr::pull(x.binned$binVals, binMax)[num.bins])

binCounts <- dplyr::pull(x.binned$binVals, binCount)

MLEbin.res <-  calcLike(negLL.PLB.binned,
                        p = -1.5,
                        w = binBreaks,
                        d = binCounts,
                        J = length(binCounts),   # = num.bins
                        vecDiff = 1)             # increase this if hit a bound

## ---- plotvec, fig.width = 5.36, fig.height = 8-------------------------------
# fig.width is 0.67 * fig.height (which is 8)
ISD_bin_plot_nonoverlapping(binBreaks = binBreaks,
                            binCounts = binCounts,
                            b.MLE = MLEbin.res$MLE,
                            b.confMin = MLEbin.res$conf[1],
                            b.confMax = MLEbin.res$conf[2])

## ---- plottib, fig.width = 5.36, fig.height = 8-------------------------------
ISD_bin_plot_nonoverlapping(binValsTibble = x.binned$binVals,
                            b.MLE = MLEbin.res$MLE,
                            b.confMin = MLEbin.res$conf[1],
                            b.confMax = MLEbin.res$conf[2],
                            yBig.inc = 500)

## ---- plotLBN, fig.height = 5.7-----------------------------------------------
LBN_bin_plot(binValsTibble = x.binned$binVals,
             b.MLE = MLEbin.res$MLE,
             b.confMin = MLEbin.res$conf[1],
             b.confMax = MLEbin.res$conf[2],
             log.xy = "",
             plot.binned.fitted = FALSE)

## ---- plotLBN2, fig.height = 5.7----------------------------------------------
LBN_bin_plot(binValsTibble = x.binned$binVals,
             b.MLE = MLEbin.res$MLE,
             b.confMin = MLEbin.res$conf[1],
             b.confMax = MLEbin.res$conf[2],
             leg.text = "(b)",
             log.xy = "x",
             plot.binned.fitted = FALSE)

## ---- plotLBN3, fig.height = 5.7----------------------------------------------
LBN_bin_plot(binValsTibble = x.binned$binVals,
             b.MLE = MLEbin.res$MLE,
             b.confMin = MLEbin.res$conf[1],
             b.confMax = MLEbin.res$conf[2],
             leg.text = "(c)",
             log.xy = "xy",
             plot.binned.fitted = FALSE)

## ---- plotLBN4, fig.height = 5.7----------------------------------------------
LBN_bin_plot(binValsTibble = x.binned$binVals,
             b.MLE = MLEbin.res$MLE,
             b.confMin = MLEbin.res$conf[1],
             b.confMax = MLEbin.res$conf[2],
             leg.text = "(c)",
             log.xy = "xy",
             plot.binned.fitted = TRUE)

## ---- plotLBN22, fig.height = 5.7---------------------------------------------
confWidth <- diff(MLEbin.res$conf)
confWidth
LBN_bin_plot(binValsTibble = x.binned$binVals,
             b.MLE = MLEbin.res$MLE,
             b.confMin = MLEbin.res$conf[1] - confWidth/2,
             b.confMax = MLEbin.res$conf[2] + confWidth/2)

