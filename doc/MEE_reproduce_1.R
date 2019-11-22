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
n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
xmax.known = 1000         # known fixed value of xmax
set.seed(42)              # To get the same observations for each run of code.

x = rPLB(n,
         b = b.known,
         xmin = xmin.known,
         xmax = xmax.known)

head(x)

## ------------------------------------------------------------------------
num.bins = 8    # Suggested number of bins for standard histogram and Llin
                #  method; Daan et al. used 8 bins.

# hAAA is the h(istrogram) results for method AAA.
# Llin method - plotting binned data on log-linear axes then fitting regression.
hLlin = Llin.method(x,
                    num.bins = num.bins)

gap.barplot.cust(hLlin$counts,
                 midpoints = hLlin$mids,
                 breakpoints = hLlin$breaks,
                 xlim = c(-10,max(hLlin$breaks)+10),
                 col = rep("grey", length(hLlin$counts)),
                 xaxs = "i",
                 yaxs = "i"
                 )

## ------------------------------------------------------------------------
eight.results = eightMethodsMEE(x, num.bins = num.bins)

## ---- fig.height = 8-----------------------------------------------------
eight.methods.plot(eight.results)

## ---- fig.height=9.4-----------------------------------------------------
MLE.plots.recommend(x = x,
                    b.MLE = eight.results$hMLE.list$b,
                    confVals.MLE = eight.results$hMLE.list$confVals,
                    hLBNbiom.list = eight.results$hLBNbiom.list)

## ---- fig.height = 5.4---------------------------------------------------
par(mai=c(0.8, 0.8, 0.2, 0.3))
MLE.plot(x,
         b = eight.results$hMLE.list$b,
         log="x")

