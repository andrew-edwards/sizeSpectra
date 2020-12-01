# likelihood.R - likelihood-related functions for statistical fitting.

##' Calculate maximum likelihood estimate and 95\% confidence
##'      interval
##'
##' Calculate maximum likelihood estimate of a parameter and its 95\% confidence
##'  interval using the profile log-likelihood method, for a given
##'  negative log-likelihood function and its arguments (other parameters and data).
##' Will likely give warnings that can safely be ignored (see
##'  `suppress.warnings` description below).
##'
##' @param negLL.fn negative log-likelihood function that take arguments
##'  (parameters and data) in ... and returns a negative
##'  log-likelihood value.
##' @param p starting point to calculate the maximum likelihood estimate
##' @param vecDiff the range over which to test the negative log-likelihood
##'  to construct the confidence interval. Default is 0.5 and a symmetric
##'  range is tested for fitting size spectra, since for movement data
##'  sets in Table 2 of Edwards (2011; 92(6):1247-1257) the intervals were
##'  symmetric, so symmetric seems a good start.
##' @param vecInc increments to try, the accuracy of the resulting bounds
##'  will depend on this. Note that a resulting interval of, say,
##'  (-2.123, -1.987) means that that interval is contained within the
##'  true 95\% interval, which is itself contained within (-2.124, -1.986).
##'  The true bounds lie between the stated lower bounds and between
##'  the stated upper bounds. So reduce `vecInc` if further accuracy is needed.
##' @param suppress.warnings If TRUE then suppress warnings from the `nlm()`
##'   calculations; for the `MEPS_IBTS_MLEbins` vignette these occur a lot, and
##'   are always:
##'   `Warning in nlm(f = negLL.fn, p = p, ...) :
##'    NA/Inf replaced by maximum positive value`. The same warning often happens in other
##'   situations also. It is due to the likelihood function blowing up,
##'   presumably when searching some very very very unlikely region of paramater space.
##'
##' @param ... further arguments (including parameters and data) to `negLL.fn()`
##' @return       list containing:
##'   * MLE: the maximum likelihood estimate
##'   * conf: the 95\% confidence interval of the MLE
##' @export
##' @author Andrew Edwards
calcLike = function(negLL.fn,
                    p,
                    vecDiff=0.5,
                    vecInc=0.001,
                    suppress.warnings = FALSE,
                    ...)
  {
    if(suppress.warnings){
      minLL = suppressWarnings(nlm(f=negLL.fn,
                                   p=p, ...)
                               )
      } else {
    minLL = nlm(f = negLL.fn,
                p = p,
                ...)
  }

    MLE = minLL$estimate
    conf = profLike(negLL.fn=negLL.fn,
                    MLE=MLE,
                    minNegLL=minLL$minimum,
                    vecDiff = vecDiff,
                    ...)
    res = list(MLE = MLE, conf = conf)
    return(res)
}

##' Profile log-likelihood method to calculate 95\% confidence interval
##'
##' Profile log-likelihood method to calculate 95\% confidence interval
##'  for a given negative log-likelihood function, maximum likelihood estimate
##'  (MLE) and minimum of of the negative log-likelihood function. Based on
##'  Hilborn and Mangel (1997), The Ecological Detective, p162.
##'
##' @param negLL.fn negative log-likelihood function that take arguments
##'   (parameters and data) in ... and returns a negative log-likelihood value
##' @param MLE maximum likelihood estimate (already calculated)
##' @param minNegLL the minimum of the negative log-likelihood function, at the
##'   MLE (by definition)
##' @param vecDiff value defining  range over which to test the negative log-likelihood
##'   to construct the confidence interval; range is `MLE` \eqn{\pm} `vecDiff`. Default is 0.5 and a symmetric
##'   range is tested for fitting size spectra, since for movement data
##'   sets in Table 2 of Edwards (2011; 92(6):1247-1257) the intervals were
##'   symmetric, so symmetric seems a good start.
##' @param vecInc increments to try, the accuracy of the resulting bounds
##'      will depend on this. Note that a resulting interval of, say,
##'     (-2.123, -1.987) means that that interval is contained within the
##'     true 95\% interval, which is itself contained within (-2.124, -1.986).
##'     The true bounds lie between the stated lower bounds and between
##'     the stated upper bounds. So reduce vecInc if further accuracy is needed.
##' @param ... further arguments to `negLL.fn()`
##' @return two-component vector of the 95\% confidence interval
##' @export
##' @author Andrew Edwards
profLike = function(negLL.fn, MLE, minNegLL, vecDiff=0.5, vecInc=0.001, ...)
    {
    vec = seq(MLE - vecDiff, MLE + vecDiff, vecInc)
                 # Values of parameter to test to obtain confidence interval

    # LLvals = vector(length=length(bvec))
    LLvals = sapply(X=vec, FUN=negLL.fn, ...)
    critVal = minNegLL  + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
    vecIn95 = vec[ LLvals < critVal ]
                      # values in 95% confidence interval
    conf = c(min(vecIn95), max(vecIn95))
    if(conf[1] == min(vec) | conf[2] == max(vec))
      { dev.new()
        plot(vec, LLvals)
        abline(h = critVal, col="red")
        stop("Need to make vecDiff larger - see R window")   # Could automate
      }
    return(conf)
}

##' Calculate MLEs and 95\% confidence intervals of `b` for separate body-mass segments
##'  of a data set
##'
##' Rather than trying to estimate a single MLE of `b` across an entire data
##' set, it may make more sense fitting the PLB distribution independently
##' across segments (distinct ranges) of the data. For example, if the full
##' range of body sizes are not collected using the same sampling protocols,
##' such as when combining phytoplankton and zooplankton data.
##'
##' @param p initial value of `b` for the numerical optimisation
##' @param w vector of length `J+1` giving the bin breaks `w_1, w_2, ..., w_{J+1}`
##' @param d vector of length `J` giving the count in each bin; must have `d_1,
##'   d_J > 0`
##' @param negLL.fn negative log-likelihood function to use - currently only
##'   works for negLL.PLB.binned
##' @param segmentBreaks the breakpoints to use to separate the range of body
##'   masses into distinct segments to be fit separately
##' @param ... further inputs to negLL.PLB.binned ***
##'
##' @return tibble with the original data with columns
##' @export
##' @author Andrew Edwards
##' @examples
##' @donttest{
##' @
##' @}
calcLikeSegments <- function(p = -1.5,
                             w,
                             d,
                             negLL.fn = negLL.PLB.binned,
                             segmentBreaks,
                             ...){
  # make data into tibble, with a row for each bin, like other examples; or
  # allow a tibble to be input
  # assign each bin a segment, according to segmentBreaks
  # fit using this, implies may need to input a vector of vecDiff's and vecInc's
  # and maybe more
#  MLEbin.res <-  calcLike(negLL.PLB.binned,
#                        p = -2.1, # PL.bMLE.binned,
#                        w = binbreaks,
#                        d = bincounts,
#                        J = length(bincounts),   # = num.bins
#                        vecDiff = 1,             # increase this if hit a bound
#                        vecInc = 1e-10)

  # save a tibble of results with columns segment, b.MLE, MLElow, MLEhigh or
  # whatever I've called it elsewhere.
  # Test on data.
  # Test on simulated data.

}




##' Calculate negative log-likelihood for the bounded power-law
##'   distribution
##'
##' Calculate the negative log-likelihood of the parameters `b`, `xmin` and
##'   `xmax` given data `x` for the PLB model. Returns the negative
##'   log-likelihood. Will be called by `nlm()` or similar. `xmin` and `xmax`
##'   are just estimated as the min and max of the data, not numerically using likelihood.
##' @param b value of `b` for which to calculate the negative log-likelihood
##' @param x vector of values of data (e.g. masses of individual fish)
##' @param n `length(x)`, have as an input to avoid repeatedly calculating it when
##'   function is called multiple times in an optimization routine
##' @param xmin minimum value of `x` to avoid repeatedly calculating
##' @param xmax maximum value of `x` to avoid repeatedly calculating
##' @param sumlogx `sum(log(x))` to avoid repeatedly calculating
##' @return negative log-likelihood of the parameters given the data
##' @export
##' @author Andrew Edwards
negLL.PLB = function(b, x, n, xmin, xmax, sumlogx)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in negLL.PLB")
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumlogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumlogx
      }
    return(neglogLL)
}

##' Calculate negative log-likelihood for the bounded power-law
##'   distribution given count data
##'
##' Calculate the negative log-likelihood of the parameters `b`, `xmin` and `xmax`
##' given count data for the PLB model. Returns the negative log-likelihood.
##' Will be called by `nlm()` or similar, but `xmin` and `xmax` will then be estimated
##' as the min of lowest bin and max of the largest, not numerically using
##' likelihood.
##'
##' For testing the MLEmid methods (using midpoints of bins), then
##' give `xmin` and `xmax` explicitly as the lowest and highest bin breaks
##' because the `x` values correspond to bins. But if `x` just represents
##' counts of discrete values then no need to specify `xmin` and `xmax`, they
##' will be automatically determined as `min(x)` and `max(x)`, respectively,
##' although it can be good to specify them to avoid repeated calculation.
##'
##' @param b value of `b` for which to calculate the negative log-likelihood.
##' @param x vector of length K corresponding to data values `x_k`, with a
##'   corresponding count `c` being the number of times that `x_k` is repeated
##' @param c vector of length `K` giving the counts `c_k` for each `k=1, 2, 3, ..., K`.
##'    Must have `c[1]>0` and `c[K]>0`, i.e. non-zero counts for first and last
##'    `x_k`. Note that the `c_k` do not have to be integer-valued.
##' @param K number of `c_k` values (length of `c`).
##' @param xmin minimum value of `x_k`, as an input to avoid repeatedly calculating.
##' @param xmax maximum value of `x_k`, as an input to avoid repeatedly calculating.
##' @param sumclogx `sum( c * log(x) )`, to avoid repeatedly calculating.
##' @return negative log-likelihood of the parameters given the data
##' @export
##' @author Andrew Edwards
negLL.PLB.counts = function(b, x, c, K=length(c), xmin=min(x), xmax=max(x),
    sumclogx = sum(c * log(x)))
  {
    if(xmin <= 0 | xmin >= xmax | length(x) != K | length(c) != K |
         c[1] == 0 | c[K] == 0 | min(c) < 0)
         stop("Parameters out of bounds in negLL.PLB.counts")
    n = sum(c)
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumclogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumclogx
      }
    return(neglogLL)
  }


##' Calculate negative log-likelihood for the bounded power-law
##'   distribution given binned data
##'
##' Calculates the negative log-likelihood of the parameters `b`, `xmin` and
##' `xmax` given binned data for the PLB model. Returns the negative
##' log-likelihood. Will be called by `nlm()` or similar, but `xmin` and `xmax`
##' will just be estimated as the minimum of lowest bin and maximum of the
##' largest bin, respectively, (since they are the MLEs), no need to do
##' numerically. Specifically this is the negative of the log-likelihood
##' function given in (A.70) of the MEE paper and (S.27) of the MEPS paper,
##' where the latter fixed a minor error in (A.75) of the MEE paper.
##'
##' @param b value of b for which to calculate the negative log-likelihood
##' @param w vector of length `J+1` giving the bin breaks `w_1, w_2, ..., w_{J+1}`
##' @param d vector of length `J` giving the count in each bin; must have `d_1,
##'   d_J > 0`
##' @param J number of bins (length of `d`)
##' @param xmin minimum value of bins, as an input to avoid repeatedly calculating
##' @param xmax maximum value of bins, as an input to avoid repeatedly calculating
##' @return negative log-likelihood of the parameters given the data
##' @export
##' @author Andrew Edwards
negLL.PLB.binned = function(b,
                            w,
                            d,
                            J = length(d),
                            xmin = min(w),
                            xmax = max(w))

  {
   # Ideally should fix those J=length(d) type things in args, that messed me up in
   # another project.
    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
         d[1] == 0 | d[J] == 0 | min(d) < 0)
         stop("Parameters out of bounds in negLL.PLB")
    n = sum(d)
    if(b != -1)
      { neglogLL = n * log( abs( w[J+1]^(b+1) - w[1]^(b+1) ) ) -
            sum( d * log( abs( w[-1]^(b+1) - w[-(J+1)]^(b+1) ) ) )
      } else
      { neglogLL = n * log( log(w[J+1]) - log(w[1]) ) -
            sum( d * log( ( log(w[-1]) - log(w[-(J+1)]) ) ) )
      }
    return(neglogLL)
  }

##' Calculate the negative log-likelihood of `b` for the PLB model, given
##'  species-specific binned data (MLEbins method). DEPRECATED -- `use negLL.PLB.bins.species()`
##'
##' Calculate the negative log-likelihood of *b* for the PLB model,
##'  given binned data where the bins can be different for each species.
##'  Returns the negative log-likelihood. This was Andy's original likelihood
##'  function but we replaced it with Mike's simpler one, so use
##'  `negLL.PLB.bins.species()`. Documented this one first, so keeping.
##'  Will be called by `nlm()` or similar, but `xmin` and `xmax` will just be estimated
##'  as the min of lowest bin and max of the largest bin (i.e. their MLEs),
##'  no need to do numerically.
##'
##' @param b value of `b` for which to calculate the negative log-likelihood
##' @param dataBinForLike table data frame (tbl_df) where each row is the count in a bin
##'  of a species, and columns (and corresponding mathematical notation in MEPS Appendix) are:
##'  * `SpecCode`: code for each species, `s`
##'  * `wmin`: lower bound of the bin, `w_\{sj\}` where `j` is the bin number
##'  * `wmax`: upper bound of the bin, `w_\{s,j+1\}`
##'  * `Number`: count in that bin for that species, `d_\{sj\}`
##'  For each species the first and last bins must be non-empty, i.e.
##'   `w_\{s1\}, w_\{s,J_s +1\} > 0`. Should ideally write code to check that before
##'    calling this function (since this gets repeatedly called).
##' @param dataBinForLikeSummary tbl_df with one row for each species, giving
##'  the minimum lower bound `w_\{s1\}` and maximum upper bound `w_\{s,J_s +1\}`
##'  and the number of counts for that species, where `J_s` is the number of
##'  bins for species `s` (won't need to explicilty specify). Columns are:
##'  * `SpecCode`: code for each species `s`
##'  * `wminSpecies`: minimum lower bound `w_\{s1\}`
##'  * `wmaxSpecies: maximum upper bound `w_\{s,J_s +1\}`
##'  * `n_s`: total number of counts for species `s` `n_s`
##' @return  negative log-likelihood of the parameters given the data
##' @author Andrew Edwards
negLL.PLB.binned.species = function(b, dataBinForLike, dataBinForLikeSummary)
  {
  # Should put these in a pre-process function:
  #    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
  #       d[1] == 0 | d[J] == 0 | min(d) < 0)
  #       stop("Parameters out of bounds in negLL.PLB")
  if(b != -1)
      {  # First component of (A.55) and then sum:
         temp1 = dplyr::mutate(dataBinForLikeSummary, comp1 = - n_s *
             log( abs( wmaxSpecies^(b+1) - wminSpecies^(b+1) ) ) )
         comp1Sum = sum(temp1$comp1)
         #
         temp2 = dplyr::mutate(dataBinForLike, comp2 = Number *
             log( abs( wmax^(b+1) - wmin^(b+1) ) ) )
         comp2Sum = sum(temp2$comp2)
         #
         neglogLL = - comp1Sum - comp2Sum
      } else
      {
        stop("NOT DONE b=-1 yet; adapt, but first correct, negLL.PLB.binned;
              though negLL.PLB.binned.species is deprecated, use
              negLL.PLB.bins.species")
      }
    return(neglogLL)
  }

##' Calculate the negative log-likelihood of `b` for the PLB model, given
##'  species-specific binned data (MLEbins method)
##'
##' Calculate the negative log-likelihood of *b* for the PLB model,
##'  given binned data where the bins can be different for each species, namely
##'  the MLEbins method derived as equations (S.18) and (S.26) in MEPS paper.
##'  Returns the negative log-likelihood.
##'  Will be called by `nlm()` or similar, but `xmin` and `xmax` will just be estimated
##'  as the min of lowest bin and max of the largest bin (i.e. their MLEs),
##'  no need to do numerically. See Supplementary Material of MEPS paper for derivation, and
##'  the vignettes for example use.
##'
##' @param b value of `b` for which to calculate the negative log-likelihood
##' @param dataBinForLike table data frame (tbl_df) where each row is the count in a bin
##'  of a species, and columns (and corresponding mathematical notation in MEPS
##'   Supplementary Material) are:
##'  * `SpecCode`: code for each species, `s`
##'  * `wmin`: lower bound of the bin, `w_\{sj\}` where `j` is the bin number
##'  * `wmax`: upper bound of the bin, `w_\{s, j+1\}`
##'  * `Number`: count in that bin for that species, `d_\{sj\}`
##'  For each species the first and last bins must be non-empty, i.e.
##'   `w_\{s1\}, w_\{s,J_s +1\} > 0`.
##' @param n total number of counts `n = \sum_\{sj\} d_\{sj\}` over all `s` and `j`
##' @param xmin maximum likelihood estimate for `xmin`, `xmin = min_\{sj\}
##'   w_\{s, 1\}`
##' @param xmax maximum likelihood estimate for `xmax`, `xmax = max_\{sj\}
##'   w_\{s, J_s+1\}`
##' @return  negative log-likelihood of the parameters given the data
##' @author Andrew Edwards
##' @export
negLL.PLB.bins.species = function(b, dataBinForLike, n, xmin, xmax)
  {
  # Would be useful to put into a pre-processing function:
  #    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
  #       d[1] == 0 | d[J] == 0 | min(d) < 0)
  #       stop("Parameters out of bounds in negLL.PLB.bins.species")
  if(b != -1)
      {  # From MEPS equation (S.18), first calculate each component in the
         # summations and then sum them:
         temp2 = dplyr::mutate(dataBinForLike,
                        comp2 = Number * log( abs( wmax^(b+1) - wmin^(b+1) ) ) )
         comp2Sum = sum(temp2$comp2)

         logLL = - n * log( abs( xmax^(b+1) - xmin^(b+1) ) ) + comp2Sum
         neglogLL = - logLL      # Negative log-likelihood
      } else
      {
        # Not fully tested, but should work:
         temp2 = dplyr::mutate(dataBinForLike,
                        comp2 = Number * log( log(wmax) - log(wmin) ) )
         comp2Sum = sum(temp2$comp2)

         logLL = - n * log( log(xmax) - log(xmin) ) + comp2Sum
         neglogLL = - logLL      # Negative log-likelihood
      }
    return(neglogLL)
  }
