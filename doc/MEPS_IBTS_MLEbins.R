## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5.7,
  fig.height = 7
)

## -----------------------------------------------------------------------------
library(sizeSpectra)
library(tibble)  # Else prints all of a tibble
data = IBTS_data
data

## ----speciesNames-------------------------------------------------------------
herringCode = dplyr::filter(specCodeNames, species == "Clupea harengus")$speccode
herringCode
spratCode = dplyr::filter(specCodeNames, species == "Sprattus sprattus")$speccode
spratCode
specCode05 = c(herringCode, spratCode)      # species codes with 0.5cm length bins

## ----dataBin------------------------------------------------------------------
dataBin = dplyr::mutate(data,
                        LngtMax = LngtClass + 1)
aa = which(dataBin$SpecCode %in% specCode05)           # row numbers for herring, sprat
dataBin[aa, "LngtMax"] = dataBin[aa, "LngtMax"] - 0.5  # subtract 0.5 cm to
                                                       # give 0.5-cm wide bins
unique(dataBin$LngtMax - dataBin$LngtClass)            # correctly just has 0.5 and 1
unique( dplyr::filter(dataBin, LngtMax - LngtClass == 0.5)$SpecCode)  # just herring,sprat

dataBin = dplyr::mutate(dataBin, wmax = LWa * LngtMax^LWb)  # calculate max body mass
                                                            # for each bin (min
                                                            # is currently bodyMass)
dataBin = dplyr::rename(dataBin, LngtMin = LngtClass)       # For consistency
dataBin = dplyr::rename(dataBin, wmin = bodyMass)

dataBin = dataBin[ , c("Year", "SpecCode", "LngtMin", "LngtMax",
                       "LWa", "LWb", "wmin", "wmax", "Number")]     # Reorder columns

range(dplyr::mutate(dataBin,
                    wminCheck = LWa * LngtMin^LWb)$wminCheck - dataBin$wmin)
                                              # Verifying that wmin is correct
                                              # (was calculated independently)
length(unique(dataBin$SpecCode))

## ----dataBinsave, eval=FALSE--------------------------------------------------
#  usethis::use_data(dataBin, overwrite = TRUE)

## ----fig.width=7.5, fig.height=6----------------------------------------------
res <- species_bins_plots()

## ----highlight----------------------------------------------------------------
dataHighlight = dplyr::filter(data,
                              SpecCode %in% c(127205, 154675))
dataHighlightSumm = dplyr::summarise(dplyr::group_by(dataHighlight,
                                                     SpecCode),
                                     minLngt = min(LngtClass),
                                     maxLngt = max(LngtClass),
                                     LWa = unique(LWa),
                                     LWb = unique(LWb))
dataHighlightSumm

## ----echo=FALSE, eval=FALSE---------------------------------------------------
#  # Plot of wmax for each species. I think it's a metric used somewhere.
#  plot(1:(dim(dataBinSpecWmax)[1]),
#       dataBinSpecWmax$maxWmax,
#       log="y",
#       xlab="Species index",
#       ylab="Wmax for each species")
#  
#  # Plot of maximum overall body size (max of max(wmax)) each year
#  maxWmaxByYear = dplyr::summarise(dplyr::group_by(dataBin, Year),
#                                   maxWmax = max(wmax))
#  smallTck=0.01
#  yLimWmaxByYear = range(pretty(c(0, max(maxWmaxByYear$maxWmax))))
#  plot(maxWmaxByYear$Year,
#       maxWmaxByYear$maxWmax,
#       xlab="Year",
#       ylab="Wmax for each year",
#       ylim=yLimWmaxByYear)
#  xTicksSmall = maxWmaxByYear$Year
#  axis(1, at = xTicksSmall, labels = rep("", length(xTicksSmall)), tck=-smallTck)
#  yTicksSmall = seq(yLimWmaxByYear[1], yLimWmaxByYear[2], by=2000)
#  axis(2, at = yTicksSmall, labels = rep("", length(yTicksSmall)), tck=-smallTck)
#  # points(2011, dplyr::filter(maxWmaxByYear, Year == 2011)$maxWmax, col="red", pch=19)

## ----MLEbins------------------------------------------------------------------
fullYears = sort(unique(dataBin$Year))
# Do a loop for each year, saving all the results in MLEbins.nSeaFung.new
for(iii in 1:length(fullYears))
  {
    dataBinForLike = dplyr::filter(dataBin,
                                   Year == fullYears[iii])
    dataBinForLike = dplyr::select(dataBinForLike,
                                   SpecCode,
                                   wmin,
                                   wmax,
                                   Number)
    n = sum(dataBinForLike$Number)
    xmin = min(dataBinForLike$wmin)
    xmax = max(dataBinForLike$wmax)

    MLEbins.nSeaFung.oneyear.new  = calcLike(negLL.fn = negLL.PLB.bins.species,
                                             p = -1.9,
                                             suppress.warnings = TRUE,
                                             dataBinForLike = dataBinForLike,
                                             n = n,
                                             xmin = xmin,
                                             xmax = xmax)

    if(iii == 1)
    {
      MLEbins.nSeaFung.new = data.frame(Year = fullYears[iii],
                                        xmin=xmin,
                                        xmax=xmax,
                                        n=n,
                                        b=MLEbins.nSeaFung.oneyear.new$MLE,
                                        confMin=MLEbins.nSeaFung.oneyear.new$conf[1],
                                        confMax=MLEbins.nSeaFung.oneyear.new$conf[2])
    } else {
      MLEbins.nSeaFung.new = rbind(MLEbins.nSeaFung.new,
                                   c(fullYears[iii],
                                     xmin,
                                     xmax,
                                     n,
                                     MLEbins.nSeaFung.oneyear.new$MLE,
                                     MLEbins.nSeaFung.oneyear.new$conf[1],
                                     MLEbins.nSeaFung.oneyear.new$conf[2]))
   }
}

# Need the standard error for weighted linear regression,
#  see eightMethods.count() for details:
MLEbins.nSeaFung.new = dplyr::tbl_df(MLEbins.nSeaFung.new)
MLEbins.nSeaFung.new = dplyr::mutate(MLEbins.nSeaFung.new,
                                     stdErr = (abs(confMin-b) +
                                               abs(confMax-b))/(2*1.96) )
MLEbins.nSeaFung.new

## ----timeseries, fig.width=7.5, fig.height=6----------------------------------
res = timeSerPlot(MLEbins.nSeaFung.new,
                  legName = "(a) MLEbins",
                  yLim = c(-2.2, -0.9),
                  xLab = "Year",
                  method = "MLEbins",
                  legPos = "bottomleft",
                  weightReg = TRUE,
                  xTicksSmallInc = 1,
                  yTicksSmallInc = 0.05)

## -----------------------------------------------------------------------------
trendResultsMLEbinsNew = dplyr::tbl_df(res)
knitr::kable(dplyr::select(trendResultsMLEbinsNew, Method, Low, Trend, High, p, Rsquared),
             digits=c(NA, 4, 4, 4, 2, 2))

## ----fig.width=7.5, fig.height=6.3--------------------------------------------
fullResults.MLEbins = MLEbins.nSeaFung.new  # Should really have just used
                                        # MLEbins..; happened to include nSeaFung early on
trend.MLEbins.new = dplyr::filter(trendResultsMLEbinsNew,
                                  Method == "MLEbins")
fullResults.MLE = dplyr::filter(fullResults, Method == "MLE")

bYears = fullResults.MLE$Year
MLE.col = "blue"
MLEbins.col = "red"
# postscript("nSeaFungCompareTrendsCol.eps", height = 6.3,
#            width = 7.5,
#            horizontal=FALSE,  paper="special")
res.MLE = timeSerPlot(fullResults.MLE,
                      legName = "",
                      xLim = range(bYears),
                      yLim = c(-1.82, -1.35),
                      xLab = "Year",
                      method = "MLE",
                      legPos = "bottomleft",
                      weightReg = TRUE,
                      bCol = MLE.col,
                      confCol = MLE.col,
                      pchVal = 19,
                      regPlot = FALSE,
                      regColNotSig = "lightblue",
                      regColSig = "darkblue",
                      xTicksSmallInc = 1,
                      yTicksSmallInc = 0.02,
                      legExtra = c("MLEbins", "MLE"),
                      legExtraCol = c(MLEbins.col, MLE.col),
                      legExtraPos = "topleft",
                      xJitter = -0.03)       # MLEbins on top as values are higher in figure

res.MLEbins.new = timeSerPlot(fullResults.MLEbins,
                              legName = "",
                              method = "MLEbins",
                              weightReg = TRUE,
                              newPlot = FALSE,
                              bCol = MLEbins.col,
                              confCol = MLEbins.col,
                              pchVal = 19,
                              regPlot = FALSE,
                              regColNotSig = "pink",
                              regColSig = "darkred",
                              xJitter = 0.03)
# dev.off()

## -----------------------------------------------------------------------------
MLEbins.res = MLEbins.nSeaFung.new
MLEbins.res = dplyr::mutate(MLEbins.res,
                            C = (b != -1 ) * (b+1) / ( xmax^(b+1) - xmin^(b+1) ) +
                                (b == -1) * 1 / ( log(xmax) - log(xmin) )
                           )
MLEbins.res = dplyr::select(MLEbins.res, -stdErr)
knitr::kable(dplyr::select(MLEbins.res, Year, xmin, xmax, n, confMin, b,
                           confMax, C),
             digits=c(0, rep(2, 7)))

## ----save, eval=FALSE---------------------------------------------------------
#  usethis::use_data(MLEbins.res, overwrite = TRUE)

