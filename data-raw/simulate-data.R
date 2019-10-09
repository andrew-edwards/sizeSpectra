## Run each line once to create results to bundle with package
##  that take too long to produce when building a vignette. Code here ensures reproducibility.

n = 1000                  # sample size
b.known = -2              # known fixed value of b
xmin.known = 1            # known fixed value of xmin
xmax.known = 1000         # known fixed value of xmax

num.reps = 10 #000          # number of times to draw sets of n random numbers.
                          #  (throwing n PLB dice num.reps times)
set.seed(42)

# Record the slope or b for each method
Llin.rep = numeric(num.reps)*NA
LT.rep = Llin.rep
LTplus1.rep = Llin.rep
LBmiz.rep = Llin.rep
LBbiom.rep = Llin.rep
LBNbiom.rep = Llin.rep
LCD.rep = Llin.rep
MLE.rep = Llin.rep
MLEfix.rep = Llin.rep  # Adding in MLE calculations when we fix xmax=xmax.known

# Record the confidence intervals (in hindsight could have maybe done lists)
Llin.rep.conf = data.frame(confMin = Llin.rep, confMax = Llin.rep)
LT.rep.conf = Llin.rep.conf
LTplus1.rep.conf = Llin.rep.conf
LBmiz.rep.conf = Llin.rep.conf
LBbiom.rep.conf = Llin.rep.conf
LBNbiom.rep.conf = Llin.rep.conf
LCD.rep.conf = Llin.rep.conf
MLE.rep.conf = Llin.rep.conf
MLEfix.rep.conf = Llin.rep.conf

MLE.rep.xmax = Llin.rep   # Also save the xmax =max(x) for each run, to see how
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
  xmax = max(x)   #TODO check, should go in function presumably

  eight.results = sizeSpectra::eightMethodsMEE(x,
                                               num.bins = num.bins,
                                               b.only = TRUE)

  # Could do more efficiently, but just modify existing code for now

  Llin.rep[iii]       = eight.results$hLlin[1]
  Llin.rep.conf[iii,] = eight.results$hLlin[2:3]

  LT.rep[iii]         = eight.results$hLT[1]
  LT.rep.conf[iii,]   = eight.results$hLT[2:3]

  LTplus1.rep[iii]    = eight.results$hLTplus1[1]
  LTplus1.rep.conf[iii,] = eight.results$hLT[2:3]

  LBmiz.rep[iii]      = eight.results$hLBmiz[1]
  LBmiz.rep.conf[iii,]= eight.results$hLBmiz[2:3]

  LBbiom.rep[iii]     = eight.results$hLBbiom[1]
  LBbiom.rep.conf[iii, ]= eight.results$hLBbiom[2:3]

  LBNbiom.rep[iii]    = eight.results$hLBNbiom[1]
  LBNbiom.rep.conf[iii, ] = eight.results$hLBNbiom[2:3]

  LCD.rep[iii]        = eight.results$hLCD[1]
  LCD.rep.conf[iii,]  = eight.results$hLCD[2:3]

  MLE.rep[iii]        = eight.results$hMLE[1]
  MLE.rep.conf[iii,]  = eight.results$hMLE[2:3]

  MLE.rep.xmax[iii] = xmax

  # MLE (maximum likelihood method) calculations, but fix xmax=xmax.known
  # TODO make a function for MLE, then take it out of eightMethodsMEE function

  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1

  PLBfix.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
      xmin=xmin, xmax=xmax.known, sumlogx=sum.log.x) #, print.level=2 )

  PLBfix.bMLE = PLBfix.minLL$estimate
  MLEfix.rep[iii] = PLBfix.bMLE

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
  MLEfix.rep.conf[iii,] = c(PLBfix.MLE.bConf[1], PLBfix.MLE.bConf[2])

}  # End for for(iii in 1:num.reps) loop

usethis::use_data("DATASET")
