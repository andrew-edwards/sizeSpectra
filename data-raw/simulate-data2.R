# Figures 4 and 5 for MEPS paper TODO... that take too long to produce when building a vignette. Code here ensures reproducibility.
# Change filename of data when decided on a format for multiple data sets.


n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin. Make it a power of 2
if(!is.wholenumber(log2(xmin.known)) | !is.wholenumber(xmin.known))
     { stop("If want xmin.known to not be an integer and not be an integer
             power of 2 then need to edit binData; may just need to think
             about whether can just set startInteger = FALSE, but will
             need to edit binData to be able to define the binWidth for
             2k method. Will need to think about this for real data; maybe
             best to just remove the first bin (and a fit a range that is
             encompassed by the binned data).")
}
xmax.known = 1000         # known fixed value of xmax

num.reps = 10000          # number of times to draw sets of n random numbers.
                          #  (throwing n PLB dice num.reps times)
set.seed(42)              # Same seed as for original simulations in
                          #  first manuscript.

binType = list(1, 5, 10, "2k")   # Must be numeric, for linear bins, or "2k"
                                 #  so define as a list
binTypes = length(binType)
binType.name = binType
binType.name[!binType %in% "2k"] = paste0("Linear ", binType[which(!binType %in% "2k")])

MLEmethod.name = c("MLEmid", "MLEbin")    # If these change then need to change likelihood
MLEmethods = length(MLEmethod.name)       #  calls below, and table output; not automatic.

# Do 3-dimensinoal array for MLEs and then an array for confMin and confMax
MLE.array = array(NA,
                  dim=c(num.reps, binTypes, MLEmethods),
                  dimnames=list(1:num.reps,
                                unlist(binType.name),
                                MLEmethod.name))
  # no need to name the rows, just index by simulation number
  # MLE.array[i,j,k] is random sample i, bin type j, MLE method j

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
                                          # show progress, can be handy
    }

  x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)

  for(j in 1:binTypes)                    # Loop over binning type
    {
      bins.list = binData(x, binWidth=binType[[j]])  # 1, 2, 5 or "2k", etc.
      num.bins = dim(bins.list$binVals)[1]

      binBreaks = bins.list$binVals[,"binMin"]$binMin   # Loses the column names
      maxOfMaxBin = bins.list$binVals[num.bins, "binMax"]$binMax
      binBreaks = c(binBreaks, maxOfMaxBin)             # Append endpoint of final bin

      binCounts = bins.list$binVals[,"binCount"]$binCount
      binMids = bins.list$binVals[,"binMid"]$binMid     # Midpoints of bins

      if(sum(!is.wholenumber(binCounts)) > 0)
        { stop("Need to adapt code for noninteger counts")
        }
                                # May need to do at some point. Though I think
                                #  negLL.PLB.counts can handle it, but don't
                                #  worry about until everything more
                                #  functionalised.

      sumCntLogMids = sum(binCounts * log(binMids))

      # MLEmid (maximum likelihood using midpoints) calculations.

      # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
      #  as a starting point for nlm for MLE of b for PLB model.
      PL.bMLE = 1/( log(min(binBreaks)) - sumCntLogMids/sum(binCounts)) - 1

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
                            sumclogx = sumCntLogMids)

      MLE.array[i,j,"MLEmid"] = MLEmid.res$MLE
      MLEconf.array[i,j,"MLEmid", ] = MLEmid.res$conf

      # MLEbin (maximum likelihood on binned data) calculations.

      MLEbin.res = calcLike(negLL.PLB.binned,
                            p = PL.bMLE,
                            w = binBreaks,
                            d = binCounts,
                            J = length(binCounts))

      MLE.array[i,j,"MLEbin"] = MLEbin.res$MLE
      MLEconf.array[i,j,"MLEbin", ] = MLEbin.res$conf

  }   # End of for(j in 1:binTypes) loop

}  # End of for(i in 1:num.reps) loop

MLE.array.known = list("n" = n,
                       "b.known" = b.known,
                       "xmin.known" = xmin.known,
                       "xmax.known" = xmax.known,
                       "num.reps" = num.reps,
                       "binType" = binType,
                       "binTypes" = binTypes,
                       "binType.name" = binType.name)

usethis::use_data(MLE.array, overwrite = TRUE)
usethis::use_data(MLEconf.array, overwrite = TRUE)
usethis::use_data(MLE.array.known, overwrite = TRUE)
