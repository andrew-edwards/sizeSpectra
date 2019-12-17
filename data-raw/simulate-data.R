# Creates default simulation and fitting results (the list data object
#  eight.results.default) for 10,000 data sets for the MEE paper data. Running
#  as part of the vignette takes too long when building the vignette.
#  See ?sizeSpectra::eight.results.default for details.

n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
xmax.known = 1000         # known fixed value of xmax

num.reps = 10000          # number of times to draw sets of n random numbers.
                          #  (throwing n PLB dice num.reps times)
set.seed(42)

# Changing Llin.rep etc to be data frame with three columns, less to save. Call
# ...df for now then remove that.
NA.vec = numeric(num.reps)*NA
Llin.rep.df = data.frame(slope = NA.vec, confMin = NA.vec, confMax = NA.vec)
LT.rep.df = Llin.rep.df
LTplus1.rep.df = Llin.rep.df
LBmiz.rep.df = Llin.rep.df
LBbiom.rep.df = Llin.rep.df
LBNbiom.rep.df = Llin.rep.df
LCD.rep.df = Llin.rep.df
MLE.rep.df = data.frame(b = NA.vec, confMin = NA.vec, confMax = NA.vec)
MLEfix.rep.df = MLE.rep.df  # Adding in MLE calculations when we fix xmax=xmax.known

MLE.rep.xmax = NA.vec      # Also save the xmax =max(x) for each run, to see how
                           #  correlates with estimate of b.

num.bins = 8    # number of bins for standard histogram and Llin method, though
                #  this is only a suggestion (and can get overridden). Daan used
                #  8 bins.
hLBmiz.num.bins = num.bins    # for mizer method
# Main loop for doing the fitting num.reps times
for(iii in 1:num.reps)
{
  # if(iii %in% seq(1000, num.reps, 1000)) print(paste("iii = ", iii))
                                   #  to show progress, but commented out here

  x = rPLB(n, b = b.known, xmin = xmin.known, xmax = xmax.known)

  log.x = log(x)                      # to avoid keep calculating
  sum.log.x = sum( log.x )
  xmin = min(x)
  xmax = max(x)

  eight.results = sizeSpectra::eightMethodsMEE(x,
                                               num.bins = num.bins,
                                               b.only = TRUE)

  # Could do more efficiently, but just modify existing code for now

  Llin.rep.df[iii, ]       = eight.results$hLlin[1:3]
#  Llin.rep.conf[iii,] = eight.results$hLlin[2:3]

  LT.rep.df[iii, ]         = eight.results$hLT[1:3]
#  LT.rep.conf[iii,]   = eight.results$hLT[2:3]

  LTplus1.rep.df[iii, ]    = eight.results$hLTplus1[1:3]
#  LTplus1.rep.conf[iii,] = eight.results$hLT[2:3]

  LBmiz.rep.df[iii, ]      = eight.results$hLBmiz[1:3]
#  LBmiz.rep.conf[iii,]= eight.results$hLBmiz[2:3]

  LBbiom.rep.df[iii, ]     = eight.results$hLBbiom[1:3]
#   LBbiom.rep.conf[iii, ]= eight.results$hLBbiom[2:3]

  LBNbiom.rep.df[iii, ]    = eight.results$hLBNbiom[1:3]
#  LBNbiom.rep.conf[iii, ] = eight.results$hLBNbiom[2:3]

  LCD.rep.df[iii, ]        = eight.results$hLCD[1:3]
#  LCD.rep.conf[iii,]  = eight.results$hLCD[2:3]

  MLE.rep.df[iii, ]        = eight.results$hMLE[1:3]
#  MLE.rep.conf[iii,]  = eight.results$hMLE[2:3]

  MLE.rep.xmax[iii] = xmax

  # MLE (maximum likelihood method) calculations, but fix xmax=xmax.known
  # Think there is now a function for MLE, so could take it out of eightMethodsMEE function

  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1

  PLBfix.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
      xmin=xmin, xmax=xmax.known, sumlogx=sum.log.x) #, print.level=2 )

  PLBfix.bMLE = PLBfix.minLL$estimate


  # 95% confidence intervals for MLE method.
  PLBfix.minNegLL = PLBfix.minLL$minimum

  # Values of b to test to obtain confidence interval. For the movement data
  #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  #  symmetric interval here.
  bvec = seq(PLBfix.bMLE - 0.5, PLBfix.bMLE + 0.5, 0.001)

  PLBfix.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
    {
        PLBfix.LLvals[i] = negLL.PLB(bvec[i], x=x, n=length(x), xmin=xmin,
            xmax=xmax.known, sumlogx=sum.log.x)
    }
  critVal = PLBfix.minNegLL  + qchisq(0.95,1)/2
                    # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLBfix.LLvals < critVal ]
                    # b values in 95% confidence interval
  PLBfix.MLE.bConf = c(min(bIn95), max(bIn95))
  if(PLBfix.MLE.bConf[1] == min(bvec) | PLBfix.MLE.bConf[2] == max(bvec))
    { windows()
      plot(bvec, PLBfix.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger for PLBfix - see R window")   # Could automate
    }

  MLEfix.rep.df[iii, ] = c(PLBfix.bMLE, PLBfix.MLE.bConf)
#  MLEfix.rep.conf[iii,] = c(PLBfix.MLE.bConf[1], PLBfix.MLE.bConf[2])

}  # End for for(iii in 1:num.reps) loop

eight.results.default <- list(Llin.rep.df    = Llin.rep.df,
                              LT.rep.df      = LT.rep.df,
                              LTplus1.rep.df = LTplus1.rep.df,
                              LBmiz.rep.df   = LBmiz.rep.df,
                              LBbiom.rep.df  = LBbiom.rep.df,
                              LBNbiom.rep.df = LBNbiom.rep.df,
                              LCD.rep.df     = LCD.rep.df,
                              MLE.rep.df     = MLE.rep.df,
                              MLEfix.rep.df  = MLEfix.rep.df,
                              MLE.rep.xmax   = MLE.rep.xmax,
                              b.known        = b.known,
                              xmin           = xmin,
                              xmax           = xmax)

usethis::use_data(eight.results.default, overwrite = TRUE)
