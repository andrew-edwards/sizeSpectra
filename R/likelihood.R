# likelihood.R - likelihood related functions for statistical fitting.

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
