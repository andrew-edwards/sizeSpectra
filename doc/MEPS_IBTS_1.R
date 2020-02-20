## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----setup--------------------------------------------------------------------
library(sizeSpectra)

## -----------------------------------------------------------------------------
dim(dataOrig)
names(dataOrig)
dataOrig[1:5,1:7]
dataOrig[1:5,8:13]

summary(dataOrig)

## -----------------------------------------------------------------------------
numAreas = length(unique(dataOrig$Area))
numAreas
colsKeep = c("Year",
             "AphiaID",
             "LngtClas",
             "CPUE_number_per_hour",
             "a",
             "b",
             "weight_g",
             "CPUE_bio_per_hour")
colsDiscard = setdiff(names(dataOrig), colsKeep)
colsDiscard

## -----------------------------------------------------------------------------
data = sizeSpectra::s_select(dataOrig, colsKeep)   # uses Sebastian Kranz's s_dplyr_funcs.r
data
# str(data)
summary(data)
min(data$CPUE_number_per_hour)

## -----------------------------------------------------------------------------
sum(data$CPUE_number_per_hour == 0)

## ----rename-------------------------------------------------------------------
if(sum( colsKeep != c("Year", "AphiaID", "LngtClas", "CPUE_number_per_hour",
    "a", "b", "weight_g", "CPUE_bio_per_hour")) > 0)
       { stop("Need to adjust renaming") }
names(data) = c("Year", "SpecCode", "LngtClass", "Number", "LWa", "LWb",
         "bodyMass", "CPUE_bio_per_hour")
# CPUE_bio_per_hour is Number * bodyMass

## ----cm-----------------------------------------------------------------------
data$LngtClass = data$LngtClass/10

## ----arrange------------------------------------------------------------------
data = dplyr::arrange(data, Year, SpecCode, LngtClass)
data

## -----------------------------------------------------------------------------
exampleSp = dplyr::filter(data, Year == 1986, SpecCode == 105814, LngtClass == 60)
exampleSp

## ---- aggregating-------------------------------------------------------------
data = dplyr::summarise(dplyr::group_by(data,
                                        Year,
                                        SpecCode,
                                        LngtClass),
                        "Number" = sum(Number)/numAreas,
                        "LWa" = unique(LWa),
                        "LWb" = unique(LWb),
                        "bodyMass" = unique(bodyMass))
data
summary(data)

dplyr::filter(data, SpecCode == 105814, Year == 1986, LngtClass == 60)

## ----bodyMass-----------------------------------------------------------------
data = dplyr::mutate(data,
                     bodyMass2 = LWa * LngtClass^LWb)
if(max(abs(data$bodyMass2 - data$bodyMass)) > 0.0001) stop("Check conversions")
data = dplyr::select(data, -bodyMass2)              # don't keep the confirming column

## -----------------------------------------------------------------------------
range(data$LngtClass)
range(data$bodyMass)
sum(data$bodyMass == 0)   # 2549
sum(data$bodyMass < 4 )   # 6893
data = dplyr::filter(data, bodyMass >= 4)
range(data$bodyMass)
data
summary(data)

## -----------------------------------------------------------------------------
sum(data$Number)

## ----lengths------------------------------------------------------------------
sort(unique(data$LngtClass))
diff(sort(unique(data$LngtClass)))

## ----halfcm-------------------------------------------------------------------
temp = dplyr::filter(data, !(SpecCode %in% c(126417, 126425)))
unique(diff(sort(unique(temp$LngtClass))))

## -----------------------------------------------------------------------------
specCodeNames
length(unique(specCodeNames$speccode))   # checking speccode are unique

## -----------------------------------------------------------------------------
data = dplyr::ungroup(data)

## ----eval=FALSE---------------------------------------------------------------
#  IBTS_data = data
#  usethis::use_data(IBTS_data, overwrite = TRUE)

## ----results='asis'-----------------------------------------------------------
data_biomass <- dplyr::mutate(data,
                              Biomass = Number * bodyMass)
knitr::kable(rbind(data_biomass[1:6,],
                   data_biomass[(nrow(data_biomass)-5):nrow(data_biomass),
                                     ]),
             digits=c(0, 0, 0, 3, 4, 4, 2, 2))

