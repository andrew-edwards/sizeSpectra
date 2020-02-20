## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## -----------------------------------------------------------------------------
library(sizeSpectra)
data = IBTS_data
data

## -----------------------------------------------------------------------------
data
unique = dim(dplyr::summarise(dplyr::group_by(data,
                                              Year,
                                              SpecCode,
                                              LngtClass),
                              count=dplyr::n()))[1]
unique

if( unique != dim(data)[1]) stop("something wrong with 'data'")

## ---- dataSumm----------------------------------------------------------------
dataSumm = dplyr::summarise(dplyr::group_by(data, Year),
                            uniqLngtClass = length(unique(LngtClass)),
                            uniqSpec = length(unique(SpecCode)))
par(mfrow=c(2,1)) #7,1))

plot(dataSumm$Year, dataSumm$uniqLngtClass, xlab="Year",
     ylab="No. unique length classes", type="o",
     ylim=c(0, max(dataSumm$uniqLngtClass)))
plot(dataSumm$Year, dataSumm$uniqSpec, xlab="Year",
     ylab="No. unique species", type="o", ylim=c(0, max(dataSumm$uniqSpec)))

## ----doEachYear, echo=TRUE, eval=FALSE----------------------------------------
#  fullYears = unique(data$Year)
#  fullResults = data.frame()
#  for(ii in fullYears){
#    eightMethodsRes = eightMethods.count(data = data,
#                                         oneYear = ii,
#                                         figName = "nSeaFung")
#    fullResults = rbind(fullResults, eightMethodsRes)
#  }
#  # usethis::use_data(fullResults, overwrite = TRUE)    # Run this to save the data.

## ---- fig.height = 7.5--------------------------------------------------------
trendResults = timeSerPlot.eight()

## ----results='asis'-----------------------------------------------------------
knitr::kable(dplyr::select(trendResults, -adjRsquared),
             digits=c(0, 4, 4, 4, 2, 2))

