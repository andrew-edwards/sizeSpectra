# likelihood.R - likelihood-related functions for statistical fitting.

##' Calculate maximum likelihood estimate and 95\% confidence
##'      interval
##'
##' Calculate maximum likelihood estimate of a parameter and its 95\% confidence
##'  interval using the profile log-likelihood method, for a given
##'  negative log-likelihood function and its arguments (other parameters and data).
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
##'   situations also.
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
##' function given in (A.70) of Edwards et al. (2017) and (S.27) of Edwards et
##' al. (TODO), where the latter fixed a minor error in (A.75) of Edwards et
##' al. (2017).
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
negLL.PLB.binned = function(b, w, d, J=length(d), xmin=min(w), xmax=max(w))
                                        # sumlogx)
  {
   # TODO fix those J=length(d) type things in args, that messed me up in
   # another project
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
##'  species-specific binned data (MLEbins method)
##'
##' Calculate the negative log-likelihood of *b* for the PLB model,
##'  given binned data where the bins can be different for each species.
##'  Returns the negative log-likelihood. TODO understand this: USES ANDY'S ORIGINAL LIKELIHOOD
##'  FUNCTION -- SEE `negLL.PLB.bins.species()` to use Mike's simpler one.
##'  Will be called by `nlm()` or similar, but `xmin` and `xmax` will just be estimated
##'  as the min of lowest bin and max of the largest bin (i.e. their MLEs),
##'  no need to do numerically. See Appendix of MEPS paper for derivation.
##'
##' @param b value of `b` for which to calculate the negative log-likelihood
##' @param dataBinForLike table data frame (tbl_df) where each row is the count in a bin
##'  of a species, and columns (and corresponding mathematical notation in MEPS Appendix) are:
##'  * `SpecCode`: code for each species, `s`
##'  * `wmin`: lower bound of the bin, `w_\{sj\}` where `j` is the bin number
##'  * `wmax`: upper bound of the bin, `w_\{s,j+1\}`
##'  * `Number`: count in that bin for that species, `d_\{sj\}`
##'  For each species the first and last bins must be non-empty, i.e.
##'   `w_\{s1\}, w_\{s,J_s +1\} > 0`. TODO-CHECK: Write code to check that before
##'    calling this function (since this gets repeatedly called). TODO math properly
##' @param dataBinForLikeSummary tbl_df with one row for each species, giving
##'  the minimum lower bound `w_\{s1\}` and maximum upper bound `w_\{s,J_s +1\}`
##'  and the number of counts for that species, where `J_s` is the number of
##'  bins for species `s` (won't need to explicilty specify). Columns are:
##'  * `SpecCode`: code for each species `s`
##'  * `wminSpecies`: minimum lower bound `w_\{s1\}`
##'  * `wmaxSpecies: maximum upper bound `w_\{s,J_s +1\}`
##'  * `n_s`: total number of counts for species `s` `n_s`
##' @return  negative log-likelihood of the parameters given the data
##' @export
##' @author Andrew Edwards
negLL.PLB.binned.species = function(b, dataBinForLike, dataBinForLikeSummary)
  {
  # TODO - think done: COPIED FROM negLL.PLB.binned for now, for which b=-1 needs correcting
  #
  # TODO - check if done MOVE THESE TO PRE-PROCESS function:
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
      { stop("NOT DONE b=-1 yet; adapt, but first correct, negLL.PLB.binned")
      }
    return(neglogLL)
  }
##' TODO DOCUMENT ME IF NECESSARY - YES, is called from MLEbins vignette
##'
##' TODO Need to clarify how this differs from `negLL.PLB.binned.species()`
##'
##' @param b temp
##' @param dataBinForLike temp
##' @param n temp
##' @param xmin temp
##' @param xmax temp
##' @return temp
##' @export
##' @author Andrew Edwards
negLL.PLB.bins.species = function(b, dataBinForLike, n, xmin, xmax)
  {
  # Calculates the negative log-likelihood of b for the PLB model,
  #  given binned data where the bins can be different for each species -- i.e.
  #  the MLEbins method. Returns the negative log-likelihood.
  #  Will be called by nlm or similar, but xmin and xmax will just be estimated
  #  as the min of lowest bin and max of the largest bin (that is their MLEs),
  #  no need to do numerically. See Appendix of second manuscript for derivation.
  #
  # COPIED FROM negLL.PLB.binned.species, for which b=-1 may need correcting
  #
  # Args:
  #   b: value of b for which to calculate the negative log-likelihood
  #   dataBinForLike: table data frame where each row is the count in a bin
  #      of a species, where the columns (and corresponding mathematical
  #      notation in Appendix) are:
  #          SpecCode: code for each species, $s$
  #          wmin: lower bound of the bin, $w_{sj}$ where $j$ is the bin number
  #          wmax: upper bound of the bin, $w_{s,j+1}$
  #          Number: count in that bin for that species, $d_{sj}.
  #   n: total number of counts $n = \Sum_{sj} d_{sj}$ over all $s$ and $j$
  #   xmin: maximum likelihood estimate for xmin, xmin = $min_{sj} w_{s,1}$
  #   xmax: maximum likelihood estimate for xmax, xmax = $max_{sj} w_{s,J_s+1}$
  #
  # Returns:
  #   negative log-likelihood of the parameters given the data.
  #
  # Useful to have in pre-processing function:
  #    if(xmin <= 0 | xmin >= xmax | length(d) != J | length(w) != J+1 |
  #       d[1] == 0 | d[J] == 0 | min(d) < 0)
  #       stop("Parameters out of bounds in negLL.PLB")
  if(b != -1)
      {  # From updated equation (A.63** number will change):

         temp2 = dplyr::mutate(dataBinForLike,
                        comp2 = Number * log( abs( wmax^(b+1) - wmin^(b+1) ) ) )
         comp2Sum = sum(temp2$comp2)

         logLL = - n * log( abs( xmax^(b+1) - xmin^(b+1) ) ) + comp2Sum
         neglogLL = - logLL      # Negative log-likelihood
      } else
      { stop("NOT DONE b=-1 yet; adapt, but first correct, negLL.PLB.binned")
      }
    return(neglogLL)
  }
