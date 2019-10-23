##' Simulate, bin and fit data using four different binning methods and two
##'  likelihood approaches
##'
##' Simulate multiple data sets from a known individual size distribution (the
##' PLB distribution), bin them using linear bins of width 1, 5 and 10, and
##' using bins that progressively double in width, and then fit each data set
##' using the MLEmid and MLEbin likelihood methods. As in Figures 4, 5, and
##' S.35-S.38, and Tables S.3-S.5 of MEPS
##' paper. See `MEPS_reproduce_2.Rmd` vignette for code for those figures and tables.
##' All simulated data sets have the same parameters for PLB and the same sample
##' size `n`. Individual data sets are not saved as they quickly take up a lot
##' of memory (would be `num.reps` $\times$ `n` random numbers, which for the
##' default values is 10^7).
##'
##' @param n sample size of each simulated data set (numeric)
##' @param b.known known fixed value of b for all simulations
##' @param xmin.known known fixed value of xmin (minimum allowable x value); currently needs to
##'   be a power of two (since makes it simpler to define the bin widths that
##'   double in size).
##' @param xmax.known known fixed value of xmax (maximum allowable x value)
##' @param num.reps number of random samples to draw, where each sample is a set
##'   of `n` random numbers (like throwing `n` PLB dice `num.reps` times)
##' @param seed seed for random number generator (default is the same as for MEE paper)
##' @param binType list containing numeric values for linear bin widths and/or
##'   "2k" (the only other option for now) for bins that double in size. Values
##'   other than the defaults have not yet been tested but should work.
##' @param vecDiffVal value to go into `profLike()` to compute confidence intervals.
##' @return list containing:
##'
##'   * MLE.array: three-dimensional array with element `[i, j, k]` representing
##'     the estimate of *b* obtained from random sample `i`, bin type `j`, and MLE
##'     method `k`. Size is `num.reps` $\times$ `length(binType)` $\times$ 2.
##'
##'   * MLEconf.array: four-dimensional array with vector
##'     `MLEconf.array[i, j, k, ]` being the confidence interval
##'     `c(confMin, confMax)` for random sample `i`, bin type `j`, and MLE method
##'     `k`.
##'
##'   * MLE.array.parameters: list containing values of
##'     `n`,
##'     `b.known`,
##'     `xmin.known`,
##'     `xmax.known`,
##'     `num.reps`,
##'     `binType`,
##'     `binTypes`,
##'     `binType.name`
##'
##' @export
##' @author Andrew Edwards
MLEbin.simulate = function(n = 1000,
                           b.known = -2,
                           xmin.known = 1,
                           xmax.known = 1000,
                           num.reps = 10000,
                           seed = 42,
                           binType = list(1, 5, 10, "2k"),
                           vecDiffVal = 0.5
                           )
{
  if(!is.wholenumber(log2(xmin.known)) | !is.wholenumber(xmin.known))
    { stop("If want xmin.known to not be an integer and not be an integer
             power of 2 then need to edit binData; may just need to think
             about whether can just set startInteger = FALSE, but will
             need to edit binData to be able to define the binWidth for
             2k method. Will need to think about this for real data; maybe
             best to just remove the first bin (and a fit a range that is
             encompassed by the binned data).")
}

  set.seed(seed)
  binTypes = length(binType)
  binType.name = binType
  binType.name[!binType %in% "2k"] = paste0("Linear ",
                                            binType[which(!binType %in% "2k")])

  MLEmethod.name = c("MLEmid", "MLEbin")    # If these change then need to change likelihood
  MLEmethods = length(MLEmethod.name)       #  calls below, and table output, so not automatic.

  # Do 3-dimensional array for MLEs and then an array for confMin and confMax
  MLE.array = array(NA,
                   dim=c(num.reps, binTypes, MLEmethods),
                   dimnames=list(1:num.reps,
                                 unlist(binType.name),
                                 MLEmethod.name))
  # No need to name the rows, just index by simulation number
  # MLE.array[i,j,k] is random sample i, bin type j, MLE method k

  # Record the confidence intervals
  MLEconf.array = array(NA,
                        dim=c(num.reps, binTypes, MLEmethods, 2),
                        dimnames=list(1:num.reps,
                                      unlist(binType.name),
                                      MLEmethod.name,
                                      c("confMin", "confMax")))
  # MLEconf.array[i,j,k, ] is confidence interval [c(confMin, confMax)] for
  #   random sample i, bin type j, MLE method k

  # Main loop for doing the fitting num.reps times
  for(i in 1:num.reps)
    {
    if(num.reps > 1000)
      {
      if(i %in% seq(1000, num.reps, 1000)) {print(paste("i = ", i))}
                                # to show progress
      }

    x = rPLB(n,
             b = b.known,
             xmin = xmin.known,
             xmax = xmax.known)

    for(j in 1:binTypes)                    # Loop over binning type
      {
        bins.list = binData(x,
                            binWidth=binType[[j]])  # 1, 2, 5 or "2k"
        num.bins = dim(bins.list$binVals)[1]

        binBreaks = bins.list$binVals[,"binMin"]$binMin   # Loses the column names
        maxOfMaxBin = bins.list$binVals[num.bins, "binMax"]$binMax
        binBreaks = c(binBreaks, maxOfMaxBin)             # Append endpoint of final bin

        binCounts = bins.list$binVals[,"binCount"]$binCount
        binMids = bins.list$binVals[,"binMid"]$binMid     # Midpoints of bins

        if(sum(!is.wholenumber(binCounts)) > 0)
          { stop("Need to adapt code for noninteger counts")
          }
                                # May be needed by someone, but not yet. Though I think
                                #  negLL.PLB.counts can handle it.

        sumCntLogMids = sum(binCounts * log(binMids))

        # MLEmid (maximum likelihood using midpoints) calculations.

        # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
        #  as a starting point for nlm for MLE of b for PLB model.
        PL.bMLE = 1/( log(min(binBreaks)) - sumCntLogMids/sum(binCounts) ) - 1

        # Note that using the min and max of binBreaks for xmin and xmax,
        #  because this method acknowledges that data are binned (but
        #  negLL.PLB.counts can act on just counts of discrete values). The
        #  min and max of binBreaks are the MLE's of xmin and xmax.
        MLEmid.res = calcLike(negLL.PLB.counts,
                              p = PL.bMLE,
                              x = binMids,
                              c = binCounts,
                              K = num.bins,
                              xmin = min(binBreaks),
                              xmax = max(binBreaks),
                              sumclogx = sumCntLogMids,
                              vecDiff = vecDiffVal)

        MLE.array[i, j, "MLEmid"] = MLEmid.res$MLE
        MLEconf.array[i, j, "MLEmid", ] = MLEmid.res$conf

        # MLEbin (maximum likelihood on binned data) calculations.

        MLEbin.res = calcLike(negLL.PLB.binned,
                              p = PL.bMLE,
                              w = binBreaks,
                              d = binCounts,
                              J = length(binCounts),
                              vecDiff = vecDiffVal)

        MLE.array[i, j, "MLEbin"] = MLEbin.res$MLE
        MLEconf.array[i, j, "MLEbin", ] = MLEbin.res$conf

    }   # End of for(j in 1:binTypes) loop

  }  # End of for(i in 1:num.reps) loop

  MLE.array.parameters = list("n" = n,
                              "b.known" = b.known,
                              "xmin.known" = xmin.known,
                              "xmax.known" = xmax.known,
                              "num.reps" = num.reps,
                              "binType" = binType,
                              "binTypes" = binTypes,
                              "binType.name" = binType.name)
  return(list(MLE.array,
              MLEconf.array,
              MLE.array.parameters))
}
