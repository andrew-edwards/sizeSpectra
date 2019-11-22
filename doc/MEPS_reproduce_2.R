## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5.7,
  fig.height = 7
  )

## ----setup---------------------------------------------------------------
library(sizeSpectra)

## ---- eval=FALSE---------------------------------------------------------
#  MLEbin.MEPS.default <- MLEbin.simulate()

## ------------------------------------------------------------------------
summary(MLEbin.MEPS.default)
dim(MLEbin.MEPS.default$MLE.array)

## ------------------------------------------------------------------------
MLEbin.MEPS.default$MLE.array[ 1:5, , ]

## ------------------------------------------------------------------------
MLEbin.MEPS.default$MLEconf.array[1:5, "Linear 1", "MLEmid", ]

## ------------------------------------------------------------------------
MLEbin.MEPS.default$MLE.array.parameters

## ------------------------------------------------------------------------
MLEmid.MLEbin.hist(MLEbin.MEPS.default)

## ------------------------------------------------------------------------
# These two lines just give a different inset for panel (e):
insetMat = matrix(rep(c(-0.01, -0.04), 4),
                  ncol=2,
                  byrow=TRUE)
insetMat[3, 1] = 0.3

MLEmid.MLEbin.conf(MLEbin.MEPS.default,
                   insetMat = insetMat)

## ---- echo=TRUE, results='asis'------------------------------------------
knitr::kable(MLEmid.MLEbin.table(MLEbin.MEPS.default))

## ---- eval=FALSE---------------------------------------------------------
#  MLEbin.MEPS.xmin16 <- MLEbin.simulate(xmin.known = 16)

## ------------------------------------------------------------------------
MLEmid.MLEbin.hist(MLEbin.MEPS.xmin16)

## ------------------------------------------------------------------------
MLEmid.MLEbin.conf(MLEbin.MEPS.xmin16)

## ---- echo=TRUE, results='asis'------------------------------------------
knitr::kable(MLEmid.MLEbin.table(MLEbin.MEPS.xmin16))

## ---- eval=FALSE---------------------------------------------------------
#  MLEbin.MEPS.cutoff16 <- MLEbin.simulate(cut.off = 16)

## ------------------------------------------------------------------------
MLEmid.MLEbin.hist(MLEbin.MEPS.cutoff16)

## ------------------------------------------------------------------------
MLEmid.MLEbin.conf(MLEbin.MEPS.cutoff16)

## ---- echo=TRUE, results='asis'------------------------------------------
knitr::kable(MLEmid.MLEbin.table(MLEbin.MEPS.cutoff16))

