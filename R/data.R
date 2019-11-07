#
#  Undocumented data sets:
#    'IBTS_data'  'MLEbins.res' 'dataBin' 'dataOrig'
#    'fullResults.MLE' 'fullResults'
#    'specCodeNames' 'trendResults'


#' Results of fitting 10,000 data sets using eight methods (MEE paper)
#'
#' Default simulation results from generating 10,000 data sets from a bounded power law,
#' with $x_min=1, $x_max=1000$ and $b=-2$, and fitting each set using the eight methods, as described in MEE
#' paper. Results are used to produce histograms and confidence intervals
#' of estimated exponent $b$ in MEE Figures 3 and 4, and statistics in Table 2
#' (except for MLEbin method).
#' @format A list with elements:
#' \describe{
#'   \item{Llin.rep.df}{Dataframe with columns `slope`, `confMin` and `confMax`,
#'   giving the estimated slope, min and max of 95% confidence intervals from
#' the Llin method, with each row representing the values for one of the 10,000 simulated data sets}
#'   \item{LT.rep.df, LTplus1.rep.df, LBmiz.rep.df, LBbiom.rep.df, LBNbiom.rep.df, LCD.rep.df}{
#'   Corresponding respective data frames for the next six methods}
#'   \item{MLE.rep.df}{Corresponding dataframe for the MLE method but with the
#' first column estimating $b$ not a slope}
#'   \item{MLEfix.rep.df}{Same as `MLE.rep.df` but for the MLEfix method, fixing $x_{max}$ to the known
#'   value of 1000 for each simulated data set (for Figures A.3, A.4 and A.5 of
#'   MEE).}
#'   \item{MLE.rep.xmax}{The corresponding estimate of $x_max$ for the MLE
#' method (namely `max(x)`) for each simulated data set.}
#'   \item{b}{Default value of `b` (namely -2) for simulations.}
#'   \item{xmin}{Default value of `xmin` (namely 1) for simulations.}
#'   \item{xmax}{Default value of `xmax` (namely 1000) for simulations.}
#' }
#' @source Generated from running `data-raw/simulate-data.R`.
"eight.results.default"

#' Default simulation results of multiple data sets binned using four binning
#' types and fit using MLEmid and MLEbin methods.
#'
#' Simulation results from generating 10,000 data sets from the PLB distribution
#' with default parameter values,
#' binning each data set using linear bins of width 1, 5 and 10, and also bins that
#' progressively double in width, and then fitting each data set using the
#' MLEmid and MLEbin likelihood methods.
#' As in Figures 4 and 5 and Table S.3 of MEPS paper. See 'MEPS_reproduce_2.Rmd'
#' vignette for code for those figures and tables, and  All simulated data sets have
#' the same parameters for PLB and the same sample size 'n'.
#' Individual data sets are not saved as they quickly take up a lot
#' of memory (would be 'num.reps' $\times$ 'n' random numbers, which
#' for the default values is 10^7).

#' @format List containing two arrays and a list of parameters. See
#'   `?sizeSpectra::MLEbin.simulate` for details.
#' @source From running `MLEbin.simulate()` in `data-raw/simulate-data2.R`.
"MLEbin.MEPS.default"

#' Simulation results of multiple data sets (with $x_{min}=16$) binned using four binning
#' types and fit using MLEmid and MLEbin methods.
#'
#' Similar to `MLEbin.MEPS.default` but with $x_{min}=16$. Produces results for
#' Figures S.35 and S.36 and Table S.4 of MEPS paper, as shown in
#' 'MEPS_reproduce_2.Rmd' vignette.
#' @format List containing two arrays and a list of parameters. See
#'   `?sizeSpectra::MLEbin.simulate` for details.
#' @source From running `MLEbin.simulate()` in `data-raw/simulate-data2.R`.
"MLEbin.MEPS.xmin16"

#' Simulation results of multiple data sets (with $x_{min}=1$ but data sampled
#' above 16) binned using four binning
#' types and fit using MLEmid and MLEbin methods.
#'
#' Similar to `MLEbin.MEPS.default` but with only sampling data above the cut
#' off value of 16. Produces results for
#' Figures S.37 and S.38 and Table S.5 of MEPS paper, as shown in
#' 'MEPS_reproduce_2.Rmd' vignette.
#' @format List containing two arrays and a list of parameters. See
#'   `?sizeSpectra::MLEbin.simulate` for details.
#' @source From running `MLEbin.simulate()` in `data-raw/simulate-data2.R`.
"MLEbin.MEPS.cutoff16"

#' Original IBTS dataset downloaded from DATRAS website
#'
#' Julia Blanchard downloaded this dataset as `ibtsQ1cpuelength.RData` that
#' contained a dataframe `q1`, which has been saved here (in
#' `data-raw/IBTS-data.R`) as `dataOrig = tibble::as_tibble(q1)`.
#' See the MEPS paper for details of the queries she used (I have some original
#' code).
#' There is now an R package for extracting database at
#' https://github.com/ices-tools-prod/icesDatras that should probably be used in
#' future for reproducibility.
#' The vignette `MEPS_IBTS_1` shows how `dataOrig` gets simplified (un-needed columns and rows are removed),
#' column names are renamed, and units are standardised, to give the more useful `IBTS_data`.
#' TODO see emails and DATRAS website for permission statement to include
#' TODO change name of dataOrig
#' @format Dataframe with 178,435 rows and 13 columns, which are (I think
#' automatically named from the DATRAS database, with `a` and `b` added from Fung et
#' al. by Julia, though see vignette for units):
##' \describe{
##'   \item{AphiaID}{}
##'   \item{Survey}{}
##'   \item{Year}{}
##'   \item{Quarter}{}
##'   \item{Area}{}
##'   \item{Species}{}
##'   \item{LngtClas}{}
##'   \item{CPUE_number_per_hour}{}
##'   \item{Taxonomic.group}{}
##'   \item{a}{}
##'   \item{b}
##'   \item{weight_g}
##'   \item{CPUE_bio_per_hour}
##' }
##' @source From running `data-raw/IBTS-data.R` on the original data file (that
##' is too big to save in this package).
"dataOrig"

#' The IBTS data after some processing as in vignette `MEPS_IBTS_1`
#'
#' @format Dataframe where each row is a unique combination of `Year`,
#' `SpecCode` and `LngtClass`. Columns are:
##' \describe{
##'   \item{Year}{Year of survey}
##'   \item{SpecCode}{Species code}
##'   \item{LngtClass}{The minimum value (cm) of the 1-cm length bin, or
##'   the 0.5-cm bin for Atlantic Herring and European Sprat}
##'   \item{Number}{Number of individuals per hour of trawling observed for that
##'   combination of `Year`, `SpecCode` and `LngtClass`.
##'   \item{LWa}{Length-weight coefficient $\alpha$ from Fung et al. (2012) for that
##'   species, as per our MEPS equation (1).}
##'   \item{LWb}{Length-weight coefficient $\beta$ from Fung et al. (2012) for that
##'   species, as per our MEPS equation (1).}
##'   \item{bodyMass}{Estimated body mass (g) for an individual of that species,
##'   assuming `LngtClass` to be the length.}
##' }
"IBTS_data"
