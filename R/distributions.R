# distributions.R - functions related to random number generation for bounded
# and unbounded power-law distributions. Also the biomass distribution function
# (MEE equations A.4 and A.8).

##' Bounded and unbounded power-law distributions
##'
##' For the unbounded and bounded power-law distributions (PL and PLB
##' respectively), probability density function (`dPL` and `dPLB`),
##' cumulative distribution function P(X <= x) (`pPL` and `pPLB`),
##' and random generation of values (`rPL` and `rPLB`), with exponent `b`, minimum `xmin` and maximum
##' (for bounded distribution) `xmax` as described in Edwards et al. (2017, Methods in Ecology and
##' Evolution, 8:57-67). Random generation uses the inverse method (e.g. p1215 of Edwards
##' 2008, Journal of Animal Ecology, 77:1212-1222). Unbounded distribution
##' included for completeness but is not used in remaining code.
##' @param x vector of values to compute the density and distribution functions.
##' @param n number of random numbers to be generated (if `length(n) > 1` then
##' generate `length(n)` values)
##' @param b exponent of the distribution (must be <-1 for unbounded)
##' @param xmin minimum bound of the distribution, `xmin > 0`
##' @param xmax maximum bound for bounded distribution, `xmax > xmin`
##' @return `dPL` and `dPLB` return vector of probability density values
##' corresponding to `x`. `pPL` and `pPLB` return vector of cumulative
##' distribution values P(X <= x) corresponding to `x`. `rPL` and `rPLB` return
##' a vector (of length `n`) of independent random draws from the distribution.
##' @name Distributions
NULL

##' @rdname Distributions
##' @export
dPL <- function(x = 1, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in dPL")
    C <- - (b+1) / xmin^(b+1)
    y <- 0 * x     # so have zeros where x < xmin
    y[x >= xmin] <- C * x[x >= xmin]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPL <- function(x = 10, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in qPL")
    y <- 0 * x     # so have zeros where x < xmin
    y[x >= xmin] <- 1 - (x[x >= xmin]/xmin)^(b+1)
    return(y)
  }

##' @rdname Distributions
##' @export
rPL <- function(n = 1, b = -2, xmin = 1)
  {
    if(b >= -1 | xmin <= 0) stop("Parameters out of bounds in rPL")
    u <- runif(n)
    y <- xmin * ( 1 - u ) ^ (1/(b+1))
    return(y)
  }

##' @rdname Distributions
##' @export
dPLB <- function(x = 1, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in dPLB")
    if(b != -1)
        { C <- (b+1) / ( xmax^(b+1) - xmin^(b+1) )
        } else
        { C <- 1/ ( log(xmax) - log(xmin) )
        }
    y <- 0 * x     # so have zeros where x < xmin or x > xmax
    y[x >= xmin & x <= xmax] <- C * x[x >= xmin & x <= xmax]^b
    return(y)
  }

##' @rdname Distributions
##' @export
pPLB <- function(x = 10, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y <- 0 * x        # so have zeros where x < xmin
    y[x > xmax] <- 1
    if(b != -1)
        {  xmintobplus1 <- xmin^(b+1)
           denom <- xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] <-
                                        ( x[x >= xmin & x <= xmax]^(b + 1) -
                                          xmintobplus1 ) / denom
        } else
        {  logxmin <- log(xmin)
           denom <- log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)
  }

##' @rdname Distributions
##' @export
rPLB <- function(n = 1, b = -2, xmin = 1, xmax = 100)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in rPLB")
    u <- runif(n)
    if(b != -1)
        { y <- ( u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
        } else
        { y <- xmax^u * xmin^(1-u)
        }
    return(y)
  }

##' Biomass distribution function from MEE equations A.4 and A.8
##'
##' Total biomass between `xmin` and `x`, assuming a bounded power-law
##' distribution of body masses between `xmin` and `xmax` and a given value of
##' exponent `b`, and a total of `n` individuals.
##' Given by MEE equations A.4 and A.8. Can then be called by `pBiomassBins` to
##' give total biomass (and normalised biomass) in each bin.
##'
##' @param x vector of values for which to calculate the total biomass between
##' `xmin` and the value
##' @param b exponent of the PLB distribution
##' @param xmin minimum bound of the distribution, `xmin > 0`
##' @param xmax maximum bound for bounded distribution, `xmax > xmin`
##' @param n
##' @return returns of vector total biomass between `xmin` and each value of `x`
##' @export
##' @author Andrew Edwards
##' @examples
##' @donttest{
##' pBiomass(x = c(1, 5, 10, 20, 50, 100),
##'          b = -2,
##'          xmin = 1,
##'          xmax = 100,
##'          n = 1000)
##' }
pBiomass <- function(x = NULL,
                     b = NULL,
                     xmin = NULL,
                     xmax = NULL,
                     n = NULL){

  if(xmin <= 0 | xmin >= xmax | min(x) < xmin | max(x) > xmax | n <= 0){
    stop("Parameters out of bounds in pBiomass")
  }

  if(b != -1){
    C <- (b+1) / ( xmax^(b+1) - xmin^(b+1) )
  } else {
    C <- 1/ ( log(xmax) - log(xmin) )
  }

  if(b != -2){
    biomass <- n * C * (x^(b+2) - xmin^(b+2)) / (b + 2)
  } else {
    biomass <- n * C * (log(x) - log(xmin))
  }

  return(biomass)
}


##' Total and normalised biomass in each bin for a fitted distribution and given
##' bin breaks
##'
##' Bin breaks are input as EITHER a single tibble `binValsTibble`
##'  with each row representing a bin, OR as a vector `binBreaks` of breaks.
##'
##' @param ... extra arguments passed to `bBiomass`: `b`, `xmin`, `xmax, `n`.
##'   Needs to go first to avoid
##'   partial matching of `b` (to be part of `...`) with whichever of
##'   `binValsTibble` or `binBreaks` is not supplied in the call to
##'   `pBiomassBins()`. Will likely want to specify these in every call (do not
##'   rely on the default).
##' @param binValsTibble tibble of binned data with each row representing a bin
##'   and with columns `wmin` and `wmax` (min and max break of each bin), or
##'   columns named `binMin` and `binMax`.
##'   Extra columns are ignored. Similar to `LBN_bin_plot()`.
##' @param binBreaks vector of bin breaks
##' @return tibble with each row corresponding to a bin, and columns `wmin`,
##'   `wmax`,  `binWidth`, `estBiomass`, and `estBiomassNorm`. If the input is a tibble with
##'   columns `wmin` and `wmax` (or `binMin` and `binMax`), then
##'   columns `estBiomass` and `estBiomassNorm` are appended to the input tibble.
##' @export
##' @author Andrew Edwards
##' @examples
##' @donttest{
##' binBreaks = c(1, 10, 20, 50, 100)
##' pBiomassBins(binBreaks = binBreaks) # uses default pBiomass() values
##'
##' # Same data as a tibble:
##' testTibble <- dplyr::tibble(binMin = binBreaks[-length(binBreaks)],
##'                             binMax = binBreaks[-1])
##' pBiomassBins(binValsTibble = testTibble)
##' @}
pBiomassBins <- function(...,
                         binValsTibble = NULL,
                         binBreaks = NULL,
                         xmin = NULL,
                         xmax = NULL){

  stopifnot(
    "Need binValsTibble OR binBreaksto be NULL" =
    (!is.null(binValsTibble) & is.null(binBreaks) ) |
    (is.null(binValsTibble) & !is.null(binBreaks) ) )

  ifelse(!is.null(binValsTibble),
         binTibble <- binValsTibble,
         # Create tibble from the vector binBreaks:
         binTibble <- dplyr::tibble(wmin = binBreaks[-length(binBreaks)],
                                    wmax = binBreaks[-1])
         )

  # Should really insist on one format of input tibble; trying to be flexible
  # and back compatible
  if("wmin" %in% names(binTibble)){
    if(is.null(xmin)){
      xmin <- min(binTibble$wmin)
    }
    if(is.null(xmax)){
      xmax <- max(binTibble$wmax)
    }
  } else {
    if(is.null(xmin)){
      xmin <- min(binTibble$binMin)
    }
    if(is.null(xmax)){
      xmax <- max(binTibble$binMax)
    }
  }


  if("wmin" %in% names(binTibble)){
    binTibble <- dplyr::mutate(binTibble,
                               binWidth = wmax - wmin,
                               estBiomass = pBiomass(x = binTibble$wmax,
                                                     xmin = xmin,
                                                     xmax = xmax,
                                                     ...) -
                                 pBiomass(x = binTibble$wmin,
                                          xmin = xmin,
                                          xmax = xmax,
                                          ...),
                               estBiomassNorm = estBiomass / binWidth)
  } else {
    binTibble <- dplyr::mutate(binTibble,
                               binWidth = binMax - binMin,
                               estBiomass = pBiomass(x = binTibble$binMax,
                                                     xmin = xmin,
                                                     xmax = xmax,
                                                     ...) -
                                 pBiomass(x = binTibble$binMin,
                                          xmin = xmin,
                                          xmax = xmax,
                                          ...),
                               estBiomassNorm = estBiomass / binWidth)
  }
  return(binTibble)
}

##' Wrapper to call `pBiomassBins()` for three values of b (MLE and conf limits)
##'
##' Only takes a tibble as the input data (unlike `pBiomassBins()` and
##' `pBiomass()`. Currently needs `wmin`, `wmax` and `Number` as columns.
##'
##' @param ... extra arguments passed to `bBiomassBins()`: `b`, `xmin`, `xmax, `n`.
##' @param binValsTibble tibble of binned data in the form required for `pBiomassBins()`
##' @param b.MLE maximum likelihood estimate of *b* (ideally from the MLEbin method)
##' @param b.confMin lower 95\% confidence limits of *b*
##' @param b.confMax upper 95\% confidence limits of *b*
##' @return `binValsTibble` with extra columns `estBiomassMLE` and
##'   `estBiomassNormMLE` for the estimated biomass and normalised biomass for
##'   `b.MLE`, extra columns `estBiomassConfMin` and `estBiomassNormConfMin` for
##'   the same but using `b.confMin`, `estBiomassConfMax` and
##'   `estBiomassNormConfMax` for `b.confMax`.
##' @export
##' @author Andrew Edwards
##' @examples
##' @donttest{
##' # see `MLEbin_recommend` vignette
##' @}
pBiomassBinsConfs <- function(...,
                              binValsTibble,
                              b.MLE,
                              b.confMin,
                              b.confMax){

  # MLE value
  binTibbleConfs <- pBiomassBins(binValsTibble = binValsTibble,
                                 b = b.MLE,
                                 n = sum(binValsTibble$Number),
                                 ...) %>%
    dplyr::rename(estBiomassMLE = estBiomass,
                  estBiomassNormMLE = estBiomassNorm)

  # Minimum of confidence interval for b
  binTibbleConfs <- pBiomassBins(binValsTibble = binTibbleConfs,
#                                 xmin = min(binTibbleConfs$wmin),
#                                 xmax = max(binTibbleConfs$wmax),
                                 b = b.confMin,
                                 n = sum(binTibbleConfs$Number),
                                 ...) %>%
    dplyr::rename(estBiomassConfMin = estBiomass,
                  estBiomassNormConfMin = estBiomassNorm)

  # Maximum of confidence interval for b
  binTibbleConfs <- pBiomassBins(binValsTibble = binTibbleConfs,
#                                 xmin = min(binTibbleConfs$wmin),
#                                 xmax = max(binTibbleConfs$wmax),
                                 b = b.confMax,
                                 n = sum(binTibbleConfs$Number),
                                 ...) %>%
    dplyr::rename(estBiomassConfMax = estBiomass,
                  estBiomassNormConfMax = estBiomassNorm)

  invisible(binTibbleConfs)
}
