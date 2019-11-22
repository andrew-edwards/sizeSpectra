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
dim(dataOrig)
names(dataOrig)

## ------------------------------------------------------------------------
numAreas = length(unique(dataOrig$Area))
colsKeep = c("Year",
             "AphiaID",
             "LngtClas",
             "CPUE_number_per_hour",
             "a",
             "b",
             "weight_g",
             "CPUE_bio_per_hour")
colsDiscard = setdiff(names(dataOrig), colsKeep)

data = sizeSpectra::s_select(dataOrig, colsKeep)   # uses Sebastian Kranz's s_dplyr_funcs.r

## ----rename--------------------------------------------------------------
if(sum( colsKeep != c("Year", "AphiaID", "LngtClas", "CPUE_number_per_hour",
    "a", "b", "weight_g", "CPUE_bio_per_hour")) > 0)
       { stop("Need to adjust renaming") }
names(data) = c("Year", "SpecCode", "LngtClass", "Number", "LWa", "LWb",
         "bodyMass", "CPUE_bio_per_hour")
# CPUE_bio_per_hour is Number * bodyMass

## ----cm------------------------------------------------------------------
data$LngtClass = data$LngtClass/10

## ----arrange-------------------------------------------------------------
data = dplyr::arrange(data, Year, SpecCode, LngtClass)

## ---- aggregating--------------------------------------------------------
data = dplyr::summarise(dplyr::group_by(data,
                                        Year,
                                        SpecCode,
                                        LngtClass),
                        "Number" = sum(Number)/numAreas,
                        "LWa" = unique(LWa),
                        "LWb" = unique(LWb),
                        "bodyMass" = unique(bodyMass))

## ------------------------------------------------------------------------
range(data$LngtClass)
range(data$bodyMass)
sum(data$bodyMass == 0)
sum(data$bodyMass < 100 )
data = dplyr::filter(data, bodyMass >= 100)
range(data$bodyMass)
data
summary(data)

## ------------------------------------------------------------------------
sum(data$Number)

## ----lengths-------------------------------------------------------------
sort(unique(data$LngtClass))

## ------------------------------------------------------------------------
data = dplyr::ungroup(data)

## ----results='asis'------------------------------------------------------
data_biomass <- dplyr::mutate(data,
                              Biomass = Number * bodyMass)
knitr::kable(rbind(data_biomass[1:6,],
                   data_biomass[(nrow(data_biomass)-5):nrow(data_biomass),
                                     ]),
             digits=c(0, 0, 0, 3, 4, 4, 2, 2))

## ---- dataSumm-----------------------------------------------------------
dataSumm = dplyr::summarise(dplyr::group_by(data, Year),
                            uniqLngtClass = length(unique(LngtClass)),
                            uniqSpec = length(unique(SpecCode)))
par(mfrow=c(2,1)) #7,1))

plot(dataSumm$Year, dataSumm$uniqLngtClass, xlab="Year",
     ylab="No. unique length classes", type="o",
     ylim=c(0, max(dataSumm$uniqLngtClass)))
plot(dataSumm$Year, dataSumm$uniqSpec, xlab="Year",
     ylab="No. unique species", type="o", ylim=c(0, max(dataSumm$uniqSpec)))

## ----speciesNames--------------------------------------------------------
herringCode = dplyr::filter(specCodeNames, species == "Clupea harengus")$speccode
herringCode
spratCode = dplyr::filter(specCodeNames, species == "Sprattus sprattus")$speccode
spratCode
specCode05 = c(herringCode, spratCode)      # species codes with 0.5cm length bins

## ----dataBin-------------------------------------------------------------
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

## ----MLEbins-------------------------------------------------------------
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

## ----timeseries, fig.width=7.5, fig.height=6-----------------------------
# postscript("../IBTS-min-100/trends100.eps",
#           height = 6, width = 7.5, horizontal = FALSE, paper = "special")
res = timeSerPlot(MLEbins.nSeaFung.new,
                  legName = "(a) MLEbins",
                  yLim = c(-3.4, -2.1),
                  xLab = "Year",
                  method = "",
                  legPos = "bottomleft",
                  weightReg = TRUE,
                  xTicksSmallInc = 1,
                  yTicksSmallInc = 0.05)
# dev.off()

## ------------------------------------------------------------------------
trendResultsMLEbinsNew = dplyr::tbl_df(res)
knitr::kable(dplyr::select(trendResultsMLEbinsNew, Method, Low, Trend, High, p, Rsquared),
             digits=c(NA, 4, 4, 4, 2, 2))

## ------------------------------------------------------------------------
MLEbins.res = MLEbins.nSeaFung.new
MLEbins.res = dplyr::mutate(MLEbins.res,
                            C = (b != -1 ) * (b+1) / ( xmax^(b+1) - xmin^(b+1) ) +
                                (b == -1) * 1 / ( log(xmax) - log(xmin) )
                           )
MLEbins.res = dplyr::select(MLEbins.res, -stdErr)
knitr::kable(dplyr::select(MLEbins.res, Year, xmin, xmax, n, confMin, b,
                           confMax, C),
             digits=c(0, rep(2, 7)))

## ----recommended---------------------------------------------------------
dataRecommend.isd = dplyr::select(dataBin,
                                  Year,
                                  wmin,
                                  wmax,
                                  Number)

data.year.list = list()                # to save results for each year
diff.ivec = vector()                   # to save i that have any cumSum !=
                                       # verify TODO
fullYears = sort(unique(dataBin$Year))
for(i in 1:length(fullYears))
  {
    data.year = dplyr::filter(dataRecommend.isd,
                              Year == fullYears[i])

    data.year = dplyr::arrange(data.year,
                               desc(wmin))
    sumNumber = sum(data.year$Number)

    # Have to do not with dplyr:
    wmin.vec = data.year$wmin
    wmax.vec = data.year$wmax
    num.vec  = data.year$Number

    countGTEwmin = rep(NA, length(num.vec)) # to do a manual count
    lowCount = countGTEwmin
    highCount = countGTEwmin

    for(iii in 1:length(countGTEwmin))
    {
      countGTEwmin[iii]    = sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
      lowCount[iii]  = sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
      highCount[iii] = sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
    }
    data.year = cbind(data.year,
                      "countGTEwmin" = countGTEwmin,
                      "lowCount" = lowCount,
                      "highCount" = highCount)
    data.year = dplyr::tbl_df(data.year)

    data.year.list[[i]] = data.year
}


## ----plotparams----------------------------------------------------------
xlim.global = c(min(dataRecommend.isd$wmin),
                max(dataRecommend.isd$wmax))

## ---- eval=FALSE---------------------------------------------------------
#  # ```
#  # {r, animation.hook = 'gifski', interval = 1.5, fig.width = 5.36, fig.height = 8}
#  ## fig.width is 0.67 * fig.height (which is 8)
#  #
#  #for(i in 1:length(fullYears))
#  #  {
#  #  ISD_bin_plot(data.year = data.year.list[[i]],
#  #              b.MLE = dplyr::filter(MLEbins.res, Year == fullYears[i])$b,
#  #              b.confMin = dplyr::filter(MLEbins.res, Year ==
#  #                                                   fullYears[i])$confMin,
#  #              b.confMax = dplyr::filter(MLEbins.res, Year ==
#  #                                                   fullYears[i])$confMax,
#  #              year = fullYears[i],
#  #              xlim = xlim.global,
#  #              xmin = dplyr::filter(MLEbins.res, Year ==
#  #                                                   fullYears[i])$xmin,
#  #              xmax = dplyr::filter(MLEbins.res, Year ==
#  #                                                fullYears[i])$xmax,
#  #             yBig.inc = 500,
#  #              ySmall.inc = 100)
#  # }
#  # ```

## ---- eval=FALSE---------------------------------------------------------
#  for(i in 1:length(fullYears))
#    {
#    postscript(paste0("../IBTS-min-100/IBTS-ISD", fullYears[i], ".eps"),
#             height = 8, width = 5.36,
#             horizontal=FALSE, paper="special")
#  
#    ISD_bin_plot(data.year = data.year.list[[i]],
#                 b.MLE = dplyr::filter(MLEbins.res, Year == fullYears[i])$b,
#                 b.confMin = dplyr::filter(MLEbins.res, Year ==
#                                                       fullYears[i])$confMin,
#                 b.confMax = dplyr::filter(MLEbins.res, Year ==
#                                                       fullYears[i])$confMax,
#                 year = fullYears[i],
#                 xlim = xlim.global,
#                 xmin = dplyr::filter(MLEbins.res, Year ==
#                                                      fullYears[i])$xmin,
#                 xmax = dplyr::filter(MLEbins.res, Year ==
#                                                   fullYears[i])$xmax,
#                 yBig.inc = 500,
#                 ySmall.inc = 100)
#    dev.off()
#  }

