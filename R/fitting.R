# fitting.R - non-likelihood related functions for statistical fitting.

##' Tests if each value in a vector is a whole number
##'
##' Test if each value in a vector is a whole number. Taken from ?is.integer and
##' documentation is added here.
##'
##' @param x vector
##' @param tol tolerance
##' @return vector of TRUE/FALSE corresponding to each element of `x`
##' @export
##' @author Andrew Edwards
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
    {
    abs(x - round(x)) < tol
    }

##' Fit size spectrum using Llin method
##'
##' Fit size spectrum using Llin method, which is plotting binned counts on
##' log-linear axes and then fitting a linear regression.
##' @param bodyMass vector of individual body masses
##' @param num.bins  number of bins to be used, though this is only a suggestion
##'   since can get over-ridden by `hist()`
##' @param binBreaks breaks for the bins to be used to bin the data and
##'   then fit the regression
##' @return list containing:
##'
##'   * mids: midpoint of bins
##'
##'   * log.counts: log(counts) in each bin
##'
##'   * counts: counts in each bin
##'
##'   * lm: results of the linear regression of `log.counts ~ mids`
##'
##'   * slope: slope of the linear regression fit
##'
##'   * breaks: bin breaks
##'
##'   * confVals: 95\% confidence interval of the fitted slope
##' @export
##' @author Andrew Edwards
Llin.method = function(bodyMass,
                       num.bins = NULL,
                       binBreaks = NULL)
    {
        if(!is.vector(bodyMass)) stop("bodyMass not a vector in Llin.method")
        if(anyNA(bodyMass)) stop("bodyMass contains NA's in Llin.method")
        if(min(bodyMass) <= 0) stop("bodyMass needs to be >0 in Llin.method")
        x = bodyMass
        #
        if(!is.null(binBreaks))
           {
            if(min(diff(binBreaks)) < 0) stop("binBreaks need to be increasing")
            breaks = binBreaks
           }  else { breaks = num.bins }
        hLlin = hist(x, breaks=breaks, plot=FALSE) # was breaks=num.bins, 4/11/15
        hLlin.mids = hLlin$mids

        hLlin.log.counts = log(hLlin$counts)
        hLlin.log.counts[ is.infinite(hLlin.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

        hLlin.lm = lm( hLlin.log.counts ~ hLlin.mids, na.action=na.omit)

        y = list(mids = hLlin.mids,
                 log.counts = hLlin.log.counts,
                 counts = hLlin$counts,
                 lm = hLlin.lm,
                 slope = hLlin.lm$coeff[2],
                 breaks = hLlin$breaks,
                 confVals = confint(hLlin.lm, "hLlin.mids",0.95))
        return(y)
    }

##' Fit biomass size spectrum with the LBNbiom and LBbiom methods
##'
##' Use the log-binning with normalisation technique (LBNbiom method) to
##' calculate the slope of the biomass size spectra. Slope is from fitting
##' a linear regression of log10(normalised biomass in bin) against
##' log10(midpoint of bin). Bins can be defined by user, else are created to
##' double in size. Also calculates slope for the LBbiom method, for which
##'   biomasses are not normalised.
##'
##' @param bodyMass vector of individual body masses
##' @param counts dataframe (or array) with first column being a body mass
##'   value, and second column being the counts of the number of individuals for
##'   that body mass. Only bodyMass or counts can be specified.
##' @param binBreaks breaks for the bins to be used to bin the data and then fit
##'   the regression. If not provided then it calculates them as bin widths that
##'   double in size (with breaks that are powers of 2) that encompass the data,
##'   resulting in binBreaks   ..., 0.25, 0.5, 1, 2, 4, 8, 16,.... as
##'   necessary.
##'
##' @param lowerCutOff body mass representing the lower cut off for the range being
##'   fit
##' @return list containing:
##'   * `indiv`: dataframe with a row for each `bodyMass` value and with columns:
##'     + 'x': original `bodyMass` values.
##'     +  `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum, maximum,
##'        and width, respectively, of the bin within which the `x` value falls.
##'        If `indiv` has >=10^6 rows then it is not saved.
##'        If `counts` was specified then, for now, an equivalent `bodyMass`
##'        vector was created and is column `x` (i.e. body masses are repeated).
##'   * binVals: dataframe with a row for each bin and columns:
##'     + `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum, maximum, and
##'        width, respectively, of the bin.
##'     + `totalBiom`: total biomass in that bin
##'     + `totalBiomNorm`: normalised total biomass in that bin, defined as `totalBiom / binWidth`
##'     +  `log10....`: `log10` of some of the above quantities
##'   * norm.lm: `lm()` result of the linear regression fit using normalised
##'       biomass in each bin
##'   * norm.slope: slope of the linear regression fit using normalised
##'       biomass in each bin
##'   * unNorm.lm: `lm()` result of the linear regression fit when not
##'       normalising the biomass in each bin
##'   * unNorm.slope: slope of the linear regression fit when not
##'       normalising the biomass in each bin
##' @export
##' @author Andrew Edwards
LBNbiom.method = function(bodyMass = NULL,
                          counts = NULL,
                          binBreaks = NULL,
                          lowerCutOff = 0)
    {
        if(!is.null(bodyMass) & !is.null(counts)) {
            stop("need only one of bodyMass or counts in LBNbiom.method") }
        if(is.null(bodyMass) & is.null(counts)) {
            stop("need bodyMass or counts in LBNbiom.method") }
        if(!is.null(bodyMass)) {
          if(!is.vector(bodyMass))stop("bodyMass not a vector in LBNbiom.method")
          if(anyNA(bodyMass)) stop("bodyMass contains NA's in LBNbiom.method")
          if(min(bodyMass) <= 0)stop("bodyMass needs to be >0 in LBNbiom.method")
          }
        if(!is.null(counts))  {
          if(dim(counts)[2] != 2)stop("counts needs two cols in LBNbiom.method")
          if(min(counts[,1]) < 0) {
              stop("body masses in counts need to be >= 0 in LBNbiom.method") }
          if(min(counts[,2]) < 0) {
              stop("numbers in counts need to be >= 0 in LBNbiom.method") }
          }
        # First wrote code for x = bodyMass, then wanted to add in the option
        #  to have counts as an input. Could write code that would make
        #  use of the counts dataframe explicitly, but actually quite easy
        #  to just create the longer vector x (though may be slightly slower
        #  computationally), to save writing extensive new code.
        if(!is.null(bodyMass)) {x = bodyMass} else
           {x = rep(counts[,1], counts[,2]) }
        #
        if(is.null(binBreaks))
           {
            binBreaks = 2^( floor(log2(min(x))) : ceiling(log2(max(x))) )
           } else
           {
            if(min(diff(binBreaks)) < 0) { stop("binBreaks need to be increasing
                                              in LBNbiom.method")}
           }
        if(!(lowerCutOff %in% c(0, binBreaks))) { stop("need lowerCutOff
                           to be 0 or one of the binBreaks in LBNbiom.method") }
           # We could allow a user to specify a lowerCutoff that was not
           #  a binBreak, but it just makes defining the binBreaks a bit more
           #  fiddly -- code could be modified if a user wished. Although in
           #  practice people would plot the binned data and then choose which
           #  points (binned counts) to ignore when fitting the regression.
        indiv = data.frame(x)       # dataframe with one row for each individual
        indiv$binMid =cut(x,
                          breaks=binBreaks,
                          right=FALSE,
                          include.lowest=TRUE,
                          labels = binBreaks[-length(binBreaks)] + 0.5*diff(binBreaks))
        indiv$binMin =cut(x,
                          breaks=binBreaks,
                          right=FALSE,
                          include.lowest=TRUE,
                          labels = binBreaks[-length(binBreaks)])
        indiv$binMax =cut(x,
                          breaks=binBreaks,
                          right=FALSE,
                          include.lowest=TRUE,
                          labels = binBreaks[-1])
        # indiv$binWidth =cut(x, breaks=binBreaks, right=FALSE,
        #    include.lowest=TRUE, labels = diff(binBreaks))
        # indiv = dplyr::mutate(indiv, binWidth = binMax - binMin)
           # Above commands avoid any problems with bins with 0 counts.
           # Don't really need all of them, but include for completeness.
        indiv$binMid = as.numeric(as.character(indiv$binMid))
        indiv$binMin = as.numeric(as.character(indiv$binMin))
        indiv$binMax = as.numeric(as.character(indiv$binMax))
           # Now calculate biomass in each bin class:
        binVals = dplyr::summarise(dplyr::group_by(indiv, binMid),
                                   binMin = unique(binMin),
                                   binMax = unique(binMax),
                                   binWidth = binMax - binMin,
                                   totalBiom = sum(x),
                                   totalBiomNorm = totalBiom / binWidth )
           # binWidth uses new columns binMax and binMin
        binVals = binVals[order(binVals$binMid),]   # order by binMid
        #
        binVals = dplyr::mutate(binVals,
                                log10binMid = log10(binMid),
                                log10totalBiom = log10(totalBiom),
                                log10totalBiomNorm = log10(totalBiomNorm))
        binVals[ is.infinite(binVals$log10totalBiom), "log10totalBiom"] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
        binVals[ is.infinite(binVals$log10totalBiomNorm), "log10totalBiomNorm"] = NA
        binVals = dplyr::mutate(binVals, aboveCutOff = (binMid > lowerCutOff))
                  # aboveCutOff is TRUE/FALSE for the regression
        unNorm.lm = lm(log10totalBiom ~ log10binMid,
                       data = dplyr::filter(binVals, aboveCutOff),
                       na.action = na.omit)
        unNorm.slope = unNorm.lm$coeff[2]
        unNorm.conf = confint(unNorm.lm,
                              "log10binMid",
                              0.95)

        norm.lm = lm( log10totalBiomNorm ~ log10binMid,
                     data = dplyr::filter(binVals, aboveCutOff),
                     na.action=na.omit)
        norm.slope = norm.lm$coeff[2]
        norm.conf = confint(norm.lm,
                            "log10binMid",
                            0.95)
        if(dim(indiv)[1] < 10^6) {       # only save indiv if not too big
          y = list(indiv = indiv,
                   binVals = binVals,
                   unNorm.lm = unNorm.lm,
                   unNorm.slope = unNorm.slope,
                   unNorm.conf = unNorm.conf,
                   norm.lm = norm.lm,
                   norm.slope = norm.slope,
                   norm.conf = norm.conf,
                   lowerCutOff = lowerCutOff)
          } else {
            y = list(binVals = binVals,
                     unNorm.lm = unNorm.lm,
                     unNorm.slope = unNorm.slope,
                     unNorm.conf = unNorm.conf,
                     norm.lm = norm.lm,
                     norm.slope = norm.slope,
                     norm.conf = norm.conf,
                     lowerCutOff = lowerCutOff)
          }
        return(y)
    }

##' Calculate the bin breaks for the LBmiz method
##'
##' Calculate the bin breaks for the LBmiz (log binning as done by the package
##'   mizer), given `xmin` and `xmax` (min and max of data) and the number of
##'   bins, `k`. To be minimised by `nlm()` to calculate `beta`.
##'
##' @param beta to be calculated, `log10(beta)` is the constant binwidth on
##'   log10 scale. `beta` is solution to
##' ```
##' 0 = beta^(k-2) * (2 * beta-1) - xmax/xmin
##'
##'   = 2 * beta^(k-1) - beta^(k-2) - xmax/xmin
##' ```
##' @param xmin minimum of data (lower bound of lowest bin)
##' @param xmax maximum of data (upper bound of highest bin)
##' @param k number of bins
##' @return value to be minimised by `nlm()`
##' @export
##' @author Andrew Edwards
LBmizbinsFun = function(beta, xmin, xmax, k)
    {
    # if( beta < 1) stop("beta needs to be >1; try reducing the number of
    #                       requested bins (k)")
    if( xmin <= 0 | xmin >= xmax | xmin >= xmax | k < 2 )
        { stop("Parameters out of bounds in LBmizbinsFun") }
    #fun = abs(beta^(k-2) * (2 * beta - 1) - xmax/xmin)   # to minimise
    # Use log scale for numerical stability
    fun = abs( (k-2)*log(beta) + log(2 * beta - 1) - log(xmax) + log(xmin))
                                         # to minimise
    return(fun)
    }

##' Construct bins that double in size and encompass the data,
##'
##' Construct bins that double in size and encompass the data, resulting in
##' `binBreaks` ..., 0.25, 0.5, 1, 2, 4, 8, 16,.... as necessary, where the
##' breaks are powers of 2. Adapting from `LBNbiom.method()`. See `binData()`
##' for a generalised version.
##'
##' @param x vector of individual values (e.g. body masses)
##' @param counts dataframe (or array) with first column being an `x` value
##'  (e.g. body mass), and second column being the `counts` of the number of
##'   individuals for that value. Only `x` or `counts` can be specified.
##' @return       list containing:
##' * indiv: dataframe with a row for each x value, with columns:
##'   + `x`: original `x` values
##'   + `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum, maximum,
##'   and width, respectively, of the bin within which the `x` value falls.   If
##'   `indiv` has >=10^6 rows then it isn't saved.  If `counts` was specified
##'   then an equivalent `x` vector is created and is column `x` (i.e. `x`
##'   values are repeated). May not be the most efficient way, but it easiest to
##'   program. May not need `indiv` returned, but it needs to be calculated
##'   anyway.
##' * `binVals`: dataframe with a row for each new bin, with columns:
##'   + `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum, maximum, and width,
##'     respectively, of the bin
##'   + `binCount`: total number of individuals in that bin
##'   + `binCountNorm`: `binCount / binWidth`
##'   +  `binSum`: sum of individual values in that bin (appropriate if `x`
##'      represents biomass, but not length)
##'   + `binSumNorm`: `binSum / binWidth`
##'   + `log10...`: `log10()` of some of the above quantities.
##' @export
##' @author Andrew Edwards
log2bins = function(x = NULL, counts = NULL)
    {
        if(!is.null(x) & !is.null(counts)) {
            stop("need only one of x or counts in log2bins") }
        if(is.null(x) & is.null(counts)) {
            stop("need x or counts in log2bins") }
        if(!is.null(x)) {
          if(!is.vector(x))stop("x not a vector in log2bins")
          if(anyNA(x)) stop("x contains NA's in log2bins")
          if(min(x) <= 0)stop("x needs to be >0 in log2bins")
          }
        if(!is.null(counts))  {
          if(dim(counts)[2] != 2)stop("counts needs two cols in log2bins")
          if(min(counts[,1]) < 0) {
              stop("x values in counts need to be >= 0 in log2bins") }
          if(min(counts[,2]) < 0) {
              stop("numbers in counts need to be >= 0 in log2bins") }
          }
        # As for LBNbiom.method(), could write code that would make
        #  use of the counts dataframe explicitly, but actually quite easy
        #  to just create the longer vector x (though may be slightly slower
        #  computationally), to save writing extensive new code.
        if(is.null(x))
           {x = rep(counts[,1], counts[,2]) }
        #
        binBreaks = 2^( floor(log2(min(x))) : ceiling(log2(max(x))) )

        indiv = data.frame(x)       # dataframe with one row for each individual
        indiv$binMid =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)] + 0.5*diff(binBreaks))
        indiv$binMin =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)])
        indiv$binMax =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-1])
        # indiv$binWidth =cut(x, breaks=binBreaks, right=FALSE,
        #    include.lowest=TRUE, labels = diff(binBreaks))
        # indiv = mutate(indiv, binWidth = binMax - binMin)
           # Above commands avoid any problems with bins with 0 counts.
           # Don't really need all of them, but include for completeness.
        indiv$binMid = as.numeric(as.character(indiv$binMid))
        indiv$binMin = as.numeric(as.character(indiv$binMin))
        indiv$binMax = as.numeric(as.character(indiv$binMax))
           # Now calculate biomass in each bin class:
        binVals = dplyr::summarise(dplyr::group_by(indiv, binMid),
            binMin = unique(binMin),
            binMax = unique(binMax),
            binWidth = binMax - binMin,
            binCount = length(x),
            binCountNorm = binCount / binWidth,
            binSum = sum(x),
            binSumNorm = binSum / binWidth )
           # binWidth uses new columns binMax and binMin
        binVals = binVals[order(binVals$binMid),]   # order by binMid
        #
        binVals = dplyr::mutate(binVals, log10binMid = log10(binMid),
            log10binCount = log10(binCount),
            log10binSum = log10(binSum),
            log10binSumNorm = log10(binSumNorm))
        binVals[ is.infinite(binVals$log10binCount),
                  "log10binCount"] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
        binVals[ is.infinite(binVals$log10binCountNorm),
                  "log10binCountNorm"] = NA
        binVals[ is.infinite(binVals$log10binSum),
                  "log10binSum"] = NA
        binVals[ is.infinite(binVals$log10binSumNorm),
                  "log10binSumNorm"] = NA
        if(dim(indiv)[1] < 10^6) {       # only save indiv if not too big
          y = list(indiv = indiv, binVals = binVals)
          } else
          {
          y = list(binVals = binVals)
          }
        return(y)
       }
##' Compute size-spectra exponents for a dataset using all eight methods
##'
##' TODO: Developed and is called in `nSea15analysis.Snw`, which was modified
##' from what was in fitting2.r. May be specific to the IBTS data set.
##'
##' @param oneYear the year of data to use, from that in the multiple years contained
##'   in `data`. TODO: that should presumably be `dataForYear`
##' @param dataForYear  TODO: I had this for `data`: local data frame that has a
##'   unique row for every combination of `Year`, `SpecCode` and `LngtClass`,
##'   with columns:
##' * `Number`: number of observed individuals of that species in that
##'   length class in that year
##' * `bodyMass`: body mass representative of such an individual, as calculated
##'   previously by `LWa * LngtClass ^ LWb`.
##' @param figName figure name, will get appended by `-Year` for each year, to
##'   create a `.png` file for each year.
##' @return data frame with one row for each method, with columns:
##' * `Year`
##' * `Method`
##' * `b`: estimate of b from that method
##' * `confMin`: lower end of 95\% confidence interval of `b` for that method
##' * `confMax`: upper end of 95\% confidence interval of `b` for that method.
##'
##' Also saves a `.png` figure called `paste(figName, "-", oneYear, ".png")`, of the
##'   fits for each of the eight methods.
##' @export
##' @author Andrew Edwards
eightMethods = function(oneYear = 1980,
  dataForYear = dplyr::filter(data, Year == oneYear), figName = "eightMethods" )
  {
  # TODO: move this to help if necessary
  # Need x, a vector of individual fish sizes (lengths or weights), which is how
  #  the original methods functions are written.
  # Now adding explicit methods in eightMethods.counts, such as
  #  Llin.method.counts, to deal explicitly with counts, and that
  #  should also work for non-integer counts. For integer
  #  counts, the expansion to give x should give the same results.
  x = rep(dataForYear$bodyMass, dataForYear$Number)
              # 3.3 million for 1980 for old nSea15, but it was quick
  log.x = log(x)
  sum.log.x = sum( log.x )
  xmin = min(x)
  xmax = max(x)

  figheight = 7     # inches, 5.6 For 4x2 figure
  figwidth = 5.7    #

  num.bins = 8   # number of bins for standard histogram and Llin method, though
                 #  this is only a suggestion (and can get overridden). Daan used
                 #  8 bins.

  # postscript("nSea1980-fitting2a.eps", height = figheight,
  #           width = figwidth, horizontal=FALSE, paper="special")
  # Postscript plots all 3million points, ends up >150Mb file. So try .png:
    png(paste(figName,
              "-",
              oneYear,
              ".png",
              sep=""),
        height = figheight,
        width = figwidth,
        res=300,,
        units="in")

  par(mfrow=c(4,2))
  oldmai = par("mai")    #  0.95625 0.76875 0.76875 0.39375  inches I think,
                         #   think may be indpt of fig size
  par(mai=c(0.4, 0.5, 0.16, 0.3))  # Affects all figures if don't change again
                                  #  Changed 3rd from 0.05 to 0.13
  mgpVals = c(1.6,0.5,0)            # mgp values   2.0, 0.5, 0

  # Notation:
  # hAAA - h(istrogram) for method AAA.

  # Llin method - plotting binned data on log-linear axes then fitting regression
  #  as done by Daan et al. 2005.
  hLlin.list = Llin.method(x, num.bins = num.bins)

  #eightMethodsRes = data.frame("Year"=oneYear, "Method"="Llin",
  #    "b" = hLlin.list$slope,
  #    "confMin"= hLlin.list$confVals[1], "confMax"= hLlin.list$confVals[2])
  # That gives hLlin.mids as the name of the row, so does this though
  hLlin.b = hLlin.list$slope
  hLlin.confMin = hLlin.list$confVals[1]
  hLlin.confMax = hLlin.list$confVals[2]
  eightMethodsRes = data.frame("Year"=oneYear,
                               "Method"="Llin",
                               "b" = hLlin.b,
                               "confMin" = hLlin.confMin,
                               "confMax"= hLlin.confMax,
                               row.names=NULL)

  plot( hLlin.list$mids, hLlin.list$log.counts,
     xlab=expression(paste("Bin midpoints for data ", italic(x))),
     ylab = "Log (count)", mgp=mgpVals)   # xlim=c(0, 400),

  lm.line(hLlin.list$mids, hLlin.list$lm)
  inset = c(0, -0.04)     # inset distance of legend
  legend("topright", paste("(a) Llin slope=", signif(hLlin.list$slope, 3)),
       bty="n", inset=inset)

  mtext(
   paste("                                                           ",
         oneYear) )
  # LT method - plotting binned data on log-log axes then fitting regression,
  #  as done by Boldt et al. 2005, natural log of counts plotted against natural
  #  log of size-class midpoints.

  # Use Llin method's binning.
  hLT.log.mids = log(hLlin.list$mids)
  hLT.log.counts = log(hLlin.list$counts)
  hLT.log.counts[ is.infinite(hLT.log.counts) ] = NA
                    # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLT.lm = lm( hLT.log.counts ~ hLT.log.mids, na.action=na.omit)
  hLT.slope = hLT.lm$coeff[2]
  hLT.conf = confint(hLT.lm, "hLT.log.mids", 0.95)

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LT"),
      "b" = hLT.slope, "confMin"= hLT.conf[1], "confMax"= hLT.conf[2],
      row.names=NULL))

  plot( hLT.log.mids, hLT.log.counts,
       xlab=expression(paste("Log (bin midpoints for data ", italic(x), ")")),
       ylab = "Log (count)", mgp=mgpVals)

  lm.line(hLT.log.mids, hLT.lm)
  legend("topright", paste("(b) LT b=", signif(hLT.slope, 3)), bty="n",
         inset=inset)

  # LTplus1 method - plotting linearly binned data on log-log axes then fitting
  #  regression of log10(counts+1) vs log10(midpoint of bins), as done by
  #  Dulvy et al. (2004).

  # Use Llin method's binning.
  hLTplus1.log10.mids = log10(hLlin.list$mids)
  hLTplus1.log10.counts = log10(hLlin.list$counts + 1)
  hLTplus1.log10.counts[ is.infinite(hLTplus1.log10.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
                  #  but the + 1 avoids this issue here
  hLTplus1.lm = lm( hLTplus1.log10.counts ~ hLTplus1.log10.mids,
      na.action=na.omit)
  hLTplus1.slope = hLTplus1.lm$coeff[2]

  hLTplus1.conf = confint(hLTplus1.lm, "hLTplus1.log10.mids", 0.95)
  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LTplus1"),
      "b" = hLTplus1.slope, "confMin"= hLTplus1.conf[1],
      "confMax"= hLTplus1.conf[2], row.names=NULL))

  plot( hLTplus1.log10.mids, hLTplus1.log10.counts,
       xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
       ylab = "Log10 (count+1)", mgp=mgpVals)

  lm.line(hLTplus1.log10.mids, hLTplus1.lm)
  legend("topright", paste("(c) LTplus1 b=", signif(hLTplus1.slope, 3)),
       bty="n", inset=inset)

  # LBmiz method - binning data using log10 bins, plotting results on natural
  #  log axes (as in mizer). Mizer does abundance size spectrum or biomass
  #  size spectrum - the latter multiplies abundance by the min of each bin
  #  (see below).

  # Construction of bins is as follows, from Finlay
  #  Scott:
  # The bins dimensions can be specified by the user by passing in min_w, max_w
  #  [min values for the lowest and highest bins] and no_w arguments [number of
  #  bins]. These are then used:
  #    w <- 10^(seq(from=log10(min_w), to=log10(max_w), length.out=no_w))
  #    dw <- diff(w)
  #    dw[no_w] <- dw[no_w-1] # Set final dw as same as penultimate bin
  #
  #The w values are the break points of the bins (the start of the bin).
  # Regarding the regression, x and w will have the same length since x
  # is just the abundance (numbers or biomass) at size w.

  hLBmiz.num.bins = num.bins

  beta = nlm(LBmizbinsFun, 2, xmin=xmin, xmax=xmax, k=hLBmiz.num.bins)$est

  # hLBmiz.bins = c(beta^(0:(k-1)) * xmin, xmax)
  hLBmiz.bins = c(beta^(0:(hLBmiz.num.bins-1)) * min(x), max(x))
   # Mizer bin specification, with final bin being same width as penultimate bin
  hLBmiz = hist(x, breaks=hLBmiz.bins, plot=FALSE)     # linear scale. Only for counts.

  # From mizer's getCommunitySlopeCode.r:
  #  "Calculates the slope of the community abundance through time by
  #  performing a linear regression on the logged total numerical abundance
  #  at weight and logged weights (natural logs, not log to base 10,
  #  are used)."  So regress log(total counts) against log(weights)
  #  (not log10 and not normalised). And it's actually
  #   on the minima of the bins (their w).

  hLBmiz.log.min.of.bins = log(hLBmiz.bins[-length(hLBmiz.bins)])# min of bins
  hLBmiz.log.counts = log(hLBmiz$counts)
  hLBmiz.log.counts[ is.infinite(hLBmiz.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLBmiz.lm = lm( hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, na.action=na.omit)

  # hLBmiz.slope = hLBmiz.lm$coeff[2]
  # Need to subtract 1, since want to work in terms of b not slopes now
  hLBmiz.b = hLBmiz.lm$coeff[2] - 1

  hLBmiz.conf = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95) - 1

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LBmiz"),
      "b" = hLBmiz.b, "confMin"= hLBmiz.conf[1],
      "confMax"= hLBmiz.conf[2], row.names=NULL))

  plot( hLBmiz.log.min.of.bins, hLBmiz.log.counts,
     xlab=expression(paste("Log (minima of bins for data ", italic(x), ")")),
     ylab = "Log (count)", mgp=mgpVals)
     # axes=FALSE, xaxs="i", yaxs="i", , xlim=c(log10(1), log10(650)),
     #  ylim=c(log10(0.7), log10(1100)))  # So axes are logged

  lm.line(hLBmiz.log.min.of.bins, hLBmiz.lm)
  legend("bottomleft", paste("(d) LBmiz b=", signif(hLBmiz.b, 3)),
         bty="n", inset=inset)

  # mizer biomass size spectra - see option (not using) in mizerBiom.eps below.


  # LBbiom method - binning data using log2 bins, calculating biomass not counts
  #  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
  #  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

  hLBNbiom.list = LBNbiom.method(x)    # Does this method and the next.

  hLBbiom.b = hLBNbiom.list[["unNorm.slope"]] - 2
  hLBbiom.conf = hLBNbiom.list[["unNorm.conf"]] - 2

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LBbiom"),
      "b" = hLBbiom.b, "confMin"= hLBbiom.conf[1],
      "confMax"= hLBbiom.conf[2], row.names=NULL))

  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiom,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["unNorm.lm"]])
  legend("bottomleft", paste("(e) LBbiom b=",
     signif(hLBbiom.b, 3)), bty="n", inset=c(-0.08, -0.04))

  # LBNbiom method - on biomass, not counts, as per Julia's 2005 paper.
  #  log2 bins of bodymass, sum the total biomass in each bin, normalise
  #  biomasses by binwidths, fit regression to log10(normalised biomass) v
  #  log10(midpoint of bin).

  # hLBNbiom.list = LBNbiom.method(x) - already done above

  hLBNbiom.b = hLBNbiom.list[["norm.slope"]] - 1
  hLBNbiom.conf = hLBNbiom.list[["norm.conf"]] - 1

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LBNbiom"),
      "b" = hLBNbiom.b, "confMin"= hLBNbiom.conf[1],
      "confMax"= hLBNbiom.conf[2], row.names=NULL))

  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiomNorm,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (normalised biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["norm.lm"]])
  legend("bottomleft", paste("(f) LBNbiom b=",
     signif(hLBNbiom.b, 3)), bty="n", inset=inset)

  # Cumulative Distribution, LCD method
  x.sorted = sort(x, decreasing=TRUE)
  logSorted = log(x.sorted)
  logProp = log((1:length(x))/length(x))

  hLCD.lm = lm(logProp ~ logSorted)   # plot(fitsortedlog10) shows
                                                   #  residuals not good
  hLCD.slope = hLCD.lm$coeff[2]

  hLCD.b = hLCD.lm$coeff[2] - 1
  hLCD.conf = confint(hLCD.lm, "logSorted", 0.95) - 1

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LCD"),
      "b" = hLCD.b, "confMin"= hLCD.conf[1],
      "confMax"= hLCD.conf[2], row.names=NULL))

  plot(logSorted, logProp, main="",
     xlab=expression(paste("Log ", italic(x))),
     ylab=expression( paste("Log (prop. of ", values > italic(x), ")")),
     mgp=mgpVals) # , axes=FALSE)
     #xlim=c(0.8, 1000), xaxs="i", ylim=c(0.0008, 1), yaxs="i",

  lm.line(logSorted, hLCD.lm, col="red")
  # murankfreq = 1 - fitsortedlog10$coeff[2]       # mu = 1 - slope
  legend("bottomleft", paste("(g) LCD b=", signif(hLCD.b, 3)), bty="n",
       inset=inset)

  # MLE (maximum likelihood method) calculations.

  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1

  PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
      xmin=xmin, xmax=xmax, sumlogx=sum.log.x) #, print.level=2 )

  PLB.bMLE = PLB.minLL$estimate

  # 95% confidence intervals for MLE method.

  PLB.minNegLL = PLB.minLL$minimum

  # Values of b to test to obtain confidence interval. For the real movement data
  #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  #  symmetric interval here.

  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.001)  # If make 0.0001 then do
                              # get an interval for raw 1980 data

  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
      {
          PLB.LLvals[i] = negLL.PLB(bvec[i], x=x, n=length(x), xmin=xmin,
              xmax=xmax, sumlogx=sum.log.x)
      }
  critVal = PLB.minNegLL  + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
                      # b values in 95% confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if(PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == max(bvec))
    { windows()
      plot(bvec, PLB.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger - see R window")   # Could automate
    }

  # MLE.rep.xmax[iii] = xmax   - not storing xmax for now
  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("MLE"),
      "b" = PLB.bMLE, "confMin"= PLB.MLE.bConf[1],
      "confMax"= PLB.MLE.bConf[2], row.names=NULL))

  # To plot rank/frequency style plot:
  plot(sort(x, decreasing=TRUE), 1:length(x), log="xy",
       xlab=expression(paste("Values, ", italic(x))),
       ylab=expression( paste("Number of ", values >= x)), mgp=mgpVals)
    # , xlim=c(2, 500), ylim=c(0.8, 40000), axes=FALSE, xaxs="i", yaxs="i",
    #  mgp=mgpVals)

  x.PLB = seq(min(x), max(x), length=1000)     # x values to plot PLB
  y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                    xmax = max(x.PLB))) * length(x)
  lines(x.PLB, y.PLB, col="red") #, lty=5)

  legend("bottomleft", paste("(h) MLE b=", signif(PLB.bMLE, 3)), bty="n",
         inset=inset)

  # To add the curves at the limits of the 95% confidence interval:
  #for(i in c(1, length(bIn95)))   # for(i in 1:length(bIn95))  to see all vals
  #    {
  #      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
  #                  xmax = max(x.PLB))) * length(x), col="red", lty=3)
  #    }

  dev.off()
  return(eightMethodsRes)
}

# Copying from .Rmd the calculations for all eight methods for MEE
# approach. Need to make a function to then use in other simulations and look at
# sensitivities. And separate from the plotting
##' Use all eight methods to fit a simple vector of body masses
##'
##' Use all eight of the methods described in MEE paper to fit size spectra to a
##' vector of body masses. Data could be lengths but then the LBmiz, LBbiom and
##' LBNbiom methods are meaningless (TODO think).
##'
##' @param x vector of individual body masses
##' @param num.bins suggested number of bins for Llin, LT and LTplus1 methods
##' @param b.only TRUE only returns slope or b plus confidence intervals, else
##'   return full details
##' @return TODO Bunch of lists
##' @export
##' @author Andrew Edwards
eightMethodsMEE <- function(x,
                            num.bins = 8,
                            b.only = FALSE){
  # Notation:
  # hAAA - h(istrogram) for method AAA.

  # Llin method
  hLlin.list = Llin.method(x,
                           num.bins = num.bins)    #**RETURN ME**TODO

  # LT method - plotting binned data on log-log axes then fitting regression,
  #  as done by Boldt et al. 2005, natural log of counts plotted against natural
  #  log of size-class midpoints.

  # Use Llin method's binning.
  hLT.log.mids = log(hLlin.list$mids)
  hLT.log.counts = log(hLlin.list$counts)
  hLT.log.counts[ is.infinite(hLT.log.counts) ] = NA
                    # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLT.lm = lm(hLT.log.counts ~ hLT.log.mids, na.action=na.omit)

  hLT.list = list(log.mids = hLT.log.mids,
                  log.counts = hLT.log.counts,
                  lm = hLT.lm,
                  slope = hLT.lm$coeff[2],
                  confVals = confint(hLT.lm, "hLT.log.mids", 0.95))
                  # breaks = hLlin$breaks,
  # LTplus1 method - plotting linearly binned data on log-log axes then fitting
  #  regression of log10(counts+1) vs log10(midpoint of bins), as done by
  #  Dulvy et al. (2004).

  # Use Llin method's binning.
  hLTplus1.log10.mids = log10(hLlin.list$mids)
  hLTplus1.log10.counts = log10(hLlin.list$counts + 1)
  hLTplus1.log10.counts[ is.infinite(hLTplus1.log10.counts) ] = NA
                    # lm can't cope with -Inf, which appear if 0 counts in a bin
                    #  though the + 1 avoids this issue here
  hLTplus1.lm = lm( hLTplus1.log10.counts ~ hLTplus1.log10.mids, na.action=na.omit)

  hLTplus1.list = list(log10.mids = hLTplus1.log10.mids,
                       log10.counts = hLTplus1.log10.counts,
                       lm = hLTplus1.lm,
                       slope = hLTplus1.lm$coeff[2],
                       confVals = confint(hLTplus1.lm, "hLTplus1.log10.mids", 0.95))
                       # breaks = hLlin$breaks,

  # LBmiz method - binning data using log10 bins, plotting results on natural
  #  log axes (as in mizer). Mizer does abundance size spectrum or biomass
  #  size spectrum - the latter multiplies abundance by the min of each bin
  #  (see below).

  # Construction of bins is as follows, from Finlay Scott:
  # The bins dimensions can be specified by the user by passing in min_w, max_w
  #  (min values for the lowest and highest bins) and no_w (number of bins)
  #  arguments. These are then used:
  #    w <- 10^(seq(from=log10(min_w), to=log10(max_w), length.out=no_w))
  #    dw <- diff(w)
  #    dw[no_w] <- dw[no_w-1] # Set final dw as same as penultimate bin
  #
  # The w values are the break points of the bins (the start of the bin).

  hLBmiz.num.bins = num.bins

  beta = nlm(LBmizbinsFun, 2, xmin=xmin, xmax=xmax, k=hLBmiz.num.bins)$est

  hLBmiz.bins = c(beta^(0:(hLBmiz.num.bins-1)) * min(x), max(x))
     # Mizer bin specification, with final bin being same width as penultimate bin

  hLBmiz = hist(x, breaks=hLBmiz.bins, plot=FALSE)     # linear scale

  # From mizer's getCommunitySlopeCode.r:
  #  "Calculates the slope of the community abundance through time by
  #  performing a linear regression on the logged total numerical abundance
  #  at weight and logged weights (natural logs, not log to base 10, are used)."
  #  So regress log(total counts) against log(weights) (not log10 and not
  #  normalised). And it's actually on the minima of the bins (their w).

  hLBmiz.log.min.of.bins = log(hLBmiz.bins[-length(hLBmiz.bins)]) # min of each bin
  hLBmiz.log.counts = log(hLBmiz$counts)
  hLBmiz.log.counts[ is.infinite(hLBmiz.log.counts) ] = NA
                    # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLBmiz.lm = lm( hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, na.action=na.omit)

  hLBmiz.list = list(log.min.of.bins = hLBmiz.log.min.of.bins,
                     log.counts = hLBmiz.log.counts,
                     lm = hLBmiz.lm,
                     slope = hLBmiz.lm$coeff[2],
                     confVals = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95))
                    # breaks = hLlin$breaks,

  # LBbiom method - binning data using log2 bins, calculating biomass (not counts)
  #  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
  #  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

  hLBNbiom.list = LBNbiom.method(x)    # Does this LBbiom and LBNbiom methods.

  # LBNbiom method - on biomass, not counts, as per Julia Blanchard's 2005 paper.
  #  log2 bins of bodymass, sum the total biomass in each bin, normalise
  #  biomasses by binwidths, fit regression to log10(normalised biomass) v
  #  log10(midpoint of bin).

  # Gets done in LBNbiom.method(x) call above


  # Cumulative Distribution, LCD method
  x.sorted = sort(x, decreasing=TRUE)
  logSorted = log(x.sorted)
  logProp = log((1:length(x))/length(x))

  hLCD.lm = lm(logProp ~ logSorted)   # plot(fitsortedlog10) shows
                                      #  residuals not good

  hLCD.list = list(logSorted = logSorted,
                   logProp = logProp,
                   lm = hLCD.lm,
                   slope = hLCD.lm$coeff[2],
                   confVals = confint(hLCD.lm, "logSorted", 0.95))
                    # breaks = hLlin$breaks,

  # MLE (maximum likelihood method) calculations.

  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  PL.bMLE = 1/( log(min(x)) - sum.log.x/length(x)) - 1

  PLB.minLL =  nlm(negLL.PLB,
                   p=PL.bMLE,
                   x=x,
                   n=length(x),
                   xmin=xmin,
                   xmax=xmax,
                   sumlogx=sum.log.x) #, print.level=2 )

  PLB.bMLE = PLB.minLL$estimate

  # 95% confidence intervals for MLE method.

  PLB.minNegLL = PLB.minLL$minimum

  # Values of b to test to obtain confidence interval. For the real movement data
  #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  #  symmetric interval here.

  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.00001)

  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
      {
        PLB.LLvals[i] = negLL.PLB(bvec[i],
                                  x=x,
                                  n=length(x),
                                  xmin=xmin,
                                  xmax=xmax,
                                  sumlogx=sum.log.x)
      }
  critVal = PLB.minNegLL  + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
                      # b values in 95% confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))

  if(PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == max(bvec))
    { windows()
      plot(bvec, PLB.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger - see R window")   # Could automate
    }

  hMLE.list = list(b = PLB.bMLE,
                   confVals = PLB.MLE.bConf)

  if(b.only){
    return(                # slope (or b), conf interval lower and upper
                           # bounds. Note that columns are unnamed (though first
                           # value seems to be named sometimes).
           list(hLlin    = c(hLlin.list$slope, hLlin.list$confVals),
                hLT      = c(hLT.list$slope, hLT.list$confVals),
                hLTplus1 = c(hLTplus1.list$slope, hLTplus1.list$confVals),
                hLBmiz   = c(hLBmiz.list$slope, hLBmiz.list$confVals),
                hLBbiom  = c(hLBNbiom.list[["unNorm.slope"]], hLBNbiom.list[["unNorm.conf"]]),
                hLBNbiom = c(hLBNbiom.list[["norm.slope"]], hLBNbiom.list[["norm.conf"]]),
                hLCD     = c(hLCD.list$slope, hLCD.list$confVals),
                hMLE     = c(hMLE.list$b, hMLE.list$confVals)))
  } else {
    return(list(hLlin.list = hLlin.list,
                hLT.list = hLT.list,
                hLTplus1.list = hLTplus1.list,
                hLBmiz.list = hLBmiz.list,
                hLBNbiom.list = hLBNbiom.list,
                hLCD.list = hLCD.list,
                hMLE.list = hMLE.list))
  }
}

##' Use the LBN and LB methods to calculate the slope of the biomass size spectra for count data
##'
##' Use the log-binning with normalisation technique to calculate the slope of
##' the biomass size spectra, for count data. Slope is from fitting
##' a linear regression of `log10(normalised biomass in bin)`
##' against `log10(midpoint of bin)`. Bins can be defined by user,
##' else are created to double in size. Also calculates slope
##' for biomasses not being normalised (LB method).
##'
##' @param valCounts valCounts: data.frame (or tbl_df) with columns `bodyMass`
##'   and `Number`, which is the count for each body mass. `bodyMass` can
##'   represent midpoints, say, of existing bins, or be the actual
##'   species-specific converted-to-bodyMass values. Number can be non-integer,
##'   which can arise from standardising, say, trawl data scaled to be per hour.
##' @param binBreaks breaks for the bins to be used to bin the data and
##'   then fit the regression. If not provided then it calculates
##'   them as bin widths that double in size that encompass the data,
##'   resulting in `binBreaks` that are powers of 2, giving  ..., 0.25, 0.5, 1,
##'   2, 4, 8, 16,...  as necessary.
##' @param lowerCutOff body mass value representing the lower cut off for
##'   the range being fit.
##' @return list containing:
##'
##'   * `valCounts2`: dataframe `valCounts` with extra columns `binMin`, the
##'       minimum of the bin into which that `bodyMass` falls, and `biomass`,
##'       the biomass corresponding to `bodyMass * Number`.
##'
##'   * `binVals`: dataframe with a row for each bin and columns:
##'     + `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum,
##'        maximum, and width, respectively, of the bin
##'     +  `totalBiom`: total biomass in that bin
##'     +  `totalBiomNorm`: normalised total biomass in that bin,
##'             defined as `totalBiom / binWidth`
##'     + `log10....`: `log10` of some of the above quantities
##'   * `norm.lm`: `lm()` result of the linear regression fit using normalised
##'       biomass in each bin
##'   * `norm.slope`: slope of the linear regression fit using normalised
##'       biomass in each bin
##'   * `unNorm.lm`: `lm()` result of the linear regression fit when not
##'       normalising the biomass in each bin
##'   *  `unNorm.slope:` slope of the linear regression fit when not
##'       normalising the biomass in each bin
##' @export
##' @author Andrew Edwards
LBNbiom.method.counts = function(valCounts, binBreaks = NULL, lowerCutOff = 0)
    {
        if(!is.data.frame(valCounts))
            { stop("valCounts not a data.frame in LBNbiom.method.counts")}
        if(anyNA(valCounts))
            { stop("valCounts contains NA's in LBNbiom.method.counts") }
        if(min(valCounts$bodyMass) <= 0)
            { stop("valCountsbodyMass needs to be >0 in LBNbiom.method.counts") }
        #   x = bodyMass
        xmin = min(valCounts$bodyMass)
        xmax = max(valCounts$bodyMass)
        #
        if(is.null(binBreaks))
           {
            binBreaks = 2^( floor(log2(xmin)) : ceiling(log2(xmax)) )
           } else
           {
            if(min(diff(binBreaks)) < 0) { stop("binBreaks need to be increasing
                                              in LBNbiom.method.counts")}
           }
        if(!(lowerCutOff %in% c(0, binBreaks))) { stop("need lowerCutOff
                           to be 0 or one of the binBreaks in LBNbiom.method.counts") }
           # We could allow a user to specify a lowerCutoff that was not
           #  a binBreak, but it just makes defining the binBreaks a bit more
           #  fiddly -- code could be modified if a user wished. Although in
           #  practice people would plot the binned data and then choose which
           #  points (binned counts) to ignore when fitting the regression.
        valCounts2 = dplyr::mutate(valCounts,
                                   binMin = binBreaks[findInterval(bodyMass,
                                                                   binBreaks,
                                                                   rightmost.close=TRUE)],
                                   biomass = bodyMass * Number)
                                        # need total biomass for each row
        if(max(valCounts2$binMin)  == binBreaks[length(binBreaks)])
           { stop("check binBreaks in LBNbiom.method.counts") }   # Shouldn't occur

        binVals = dplyr::summarise(dplyr::group_by(valCounts2, binMin),
                                   binCount = sum(Number),
                                   totalBiom = sum(biomass))
        # No guarantee that every binMin value hLlBNbiom.mids shows up here
        missing1 = setdiff(binBreaks[-length(binBreaks)], binVals$binMin)
        missing = cbind(binMin = missing1,
                        binCount = rep(0, length(missing1)),
                        totalBiom = rep(0, length(missing1)))
        binVals = rbind(binVals, missing)
        binVals = tbl_df(binVals)
        binVals = dplyr::arrange(binVals, binMin)
        if( max( abs( binBreaks[-length(binBreaks)] - binVals$binMin) ) > 0)
            { stop("check binVals in LBNbiom.method.counts") }
        # So all the breaks except final show up as binMin, and have been ordered.
        #  Therefore to add binMax just do:
        binVals = dplyr::mutate(binVals,
                                binMax = binBreaks[-1],
                                binWidth = binMax - binMin,
                                binMid = binMin + binWidth/2,
                                totalBiomNorm = totalBiom / binWidth )
        binVals = dplyr::mutate(binVals,
                                log10binMid = log10(binMid),
                                log10totalBiom = log10(totalBiom),
                                log10totalBiomNorm = log10(totalBiomNorm))
        binVals[ is.infinite(binVals$log10totalBiom),
                "log10totalBiom"] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
        binVals[ is.infinite(binVals$log10totalBiomNorm),
                "log10totalBiomNorm"] = NA
        binVals = dplyr::mutate(binVals, aboveCutOff = (binMid > lowerCutOff))
                  # aboveCutOff is TRUE/FALSE for the regression
        unNorm.lm = lm(log10totalBiom ~ log10binMid,
                       data = dplyr::filter(binVals, aboveCutOff),
                       na.action=na.omit)
        unNorm.slope = unNorm.lm$coeff[2]
        unNorm.conf = confint(unNorm.lm, "log10binMid", 0.95)

        norm.lm = lm(log10totalBiomNorm ~ log10binMid,
                     data = dplyr::filter(binVals, aboveCutOff),
                     na.action=na.omit)
        norm.slope = norm.lm$coeff[2]
        norm.conf = confint(norm.lm, "log10binMid", 0.95)
        y = list(valCounts2 = valCounts2,
                 binVals = binVals,
                 unNorm.lm = unNorm.lm,
                 unNorm.slope = unNorm.slope,
                 unNorm.conf = unNorm.conf,
                 norm.lm = norm.lm,
                 norm.slope = norm.slope,
                 norm.conf = norm.conf,
                 lowerCutOff = lowerCutOff)
        return(y)
}

##' The Llin method for count data

##' The Llin method, which is plotting binned counts on log-linear
##'  axes and then fitting a regression, for count data.
##'
##' @param valCounts   valCounts: data.frame (can be tbl_df) with columns `bodyMass`
##'    and `Number` (which is the count for each body mass). `bodyMass` can
##'    represent midpoints, say, of existing bins, or be the actual
##'    species-specific converted-to-bodyMass values. `Number` can be
##'    non-integer, which can arise from standardising, say, trawl data to be per
##'    hour.
##'

##' @param num.bins number of bins to be used, though this is only a suggestion
##'   since can get over-ridden by `hist()`. Need to specify `num.bins` OR `binBreaks`.
##' @param binBreaks breaks for the bins to be used to bin the data and
##'   then fit the regression.
##' @return list containing:
##'   * `mids`: midpoint of bins
##'   * `log.counts`: `log(counts)` in each bin
##'   * `counts`: counts in each bin
##'   * `lm`: results of the linear regression of `log.counts ~ mids`
##'   * `slope`: slope of the linear regression fit
##'   * `breaks`: bin breaks
##'   * `confVals`: 95\% confidence interval of the fitted slope
##'   * `stdError`: standard error of the fitted slope
##'
##' @export
##' @author Andrew Edwards
Llin.method.counts = function(valCounts, num.bins = NULL, binBreaks = NULL)
    {
        if(!is.data.frame(valCounts))
            { stop("valCounts not a data.frame in Llin.method.counts")}
        if(anyNA(valCounts))
            { stop("valCounts contains NA's in Llin.method.counts") }
        if(min(valCounts$bodyMass) <= 0)
            { stop("valCountsbodyMass needs to be >0 in Llin.method.counts") }
        if(is.null(binBreaks) & is.null(num.bins))
            { stop("need binBreaks or num.bins in Llin.method.counts") }
        if(!is.null(binBreaks) & !is.null(num.bins))
            { stop("need one of binBreaks OR num.bins in Llin.method.counts") }
        if(!is.null(binBreaks))      # use to get the bin breaks
           {
            if(min(diff(binBreaks)) < 0) stop("binBreaks need to be increasing")
            hLlin.temp = hist(valCounts$bodyMass, breaks = binBreaks, plot=FALSE)
                      # will give error if binBreaks don't span bodyMass values
           }  else    # { breaks = num.bins }
           {
            hLlin.temp = hist(valCounts$bodyMass, breaks = num.bins, plot=FALSE)
           }          # Still lets hist select the breaks; num.bins is a guide
        breaks = hLlin.temp$breaks
        hLlin.mids = hLlin.temp$mids

        valCounts2 = dplyr::mutate(valCounts,
            binMid = hLlin.mids[findInterval(bodyMass, breaks,
                rightmost.close=TRUE)])
        binVals = dplyr::summarise(dplyr::group_by(valCounts2, binMid), binCount = sum(Number))
        # No guarantee that every binMid hLlin.mids shows up here (unlikely,
        #  given linear binning and long-tailed data). Even though need to
        #  exclude 0 counts (since get logged) from lm, still best to return
        #  values for all bins, especially since the bins can be specified.
        missing1 = setdiff(hLlin.mids, binVals$binMid)
        missing = cbind(binMid = missing1, binCount = rep(0, length(missing1)))
        binVals = rbind(binVals, missing)
        binVals = dplyr::tbl_df(binVals)
        binVals = dplyr::arrange(binVals, binMid)
        if( max( abs( hLlin.mids - binVals$binMid) ) > 0)
            { stop("check binVals in Llin.method.counts") }
        hLlin.log.counts = log(binVals$binCount)
                  # binVals$binCount was hLlin$counts
        hLlin.log.counts[ is.infinite(hLlin.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

        hLlin.lm = lm( hLlin.log.counts ~ hLlin.mids, na.action=na.omit)

        y = list(mids = hLlin.mids, log.counts = hLlin.log.counts,
            counts = binVals$binCount, lm = hLlin.lm, slope = hLlin.lm$coeff[2],
            breaks = breaks,
            confVals = confint(hLlin.lm, "hLlin.mids",0.95),
            stdErr = coef(summary(hLlin.lm))["hLlin.mids", "Std. Error"])
        return(y)
       }


##' Compute exponents for the IBTS dataset using all eight methods, explicitly
##'   using the counts which can be non-integer.
##'
##' Compute exponents for the IBTS dataset using all eight methods, explicitly
##'   using the counts (number of each bodyMass), which can be non-integer.
##'
##' In `eightMethods()`, had to expand data to get a vector, `x`, of individual
##' fish sizes (lengths or weights), which is how the original methods functions
##' are written.
##' Now adding explicit methods here, such as
##'  `Llin.method.counts()`, to deal explicitly with counts, and that
##'  should also work for non-integer counts. For integer
##'  counts, the original expansion to give `x` should give the same results.
##' @param data local data frame that has a unique row for every combination of
##' `Year`, `SpecCode` and `LngtClass`. The `Number` column is  the number of
##' observed individuals of that species in that length class in that
##' year. `bodyMass` is the body mass representative of such an individual, as
##' calculated previously by `LWa * LngtClass ^ LWb`.
##' @param oneYear the year of data to use, from that in the multiple years
##'   contained in data
##' @param figName figure name, will get appended by `-Year` for each year, to create
##'   a `.png` plot for each year.
##' @return data frame with one row for each method, with columns:
##'   * `Year`
##'   * `Method`
##'   * `b` (estimate of *b* from that method)
##'   * `confMin` (lower end of 95\% confidence interval of *b* for that method)
##'   * `confMax` (upper end of 95\% confidence interval of *b* for that method).
##'  a .png figure `paste(figName, "-", oneYear, ".png")` of the fits for each of the eight methods
##'
##' @export
##' @author Andrew Edwards
eightMethods.count = function(data = data, oneYear = 1980,
  figName = "eightMCnt" )
  {
  dataForYear = dplyr::filter(data, Year == oneYear)
  valCounts = dplyr::select(dataForYear, bodyMass, Number)  # Just Number of individuals
                            # with each bodyMass (Number can be non-integer).
                            # Don't worry about species here.
  valCounts = ungroup(valCounts)      # Else it retains group info, affecting
                                      #  results for, e.g., findInterval(),
                                      #  so safer to ungroup.
  xmin = min(valCounts$bodyMass)
  xmax = max(valCounts$bodyMass)

  figheight = 7 # 5.6     For 4x2 figure
  figwidth = 5.7    # 5.7 inches for JAE

  num.bins = 8   # number of bins for standard histogram and Llin method, though
                 #  this is only a suggestion (and can get overridden). Daan used
                 #  8 bins.

  png(paste(figName, "Eight-", oneYear, ".png", sep=""), height = figheight,
             width = figwidth, res=300,, units="in")

  par(mfrow=c(4,2))
  oldmai = par("mai")
  par(mai=c(0.4, 0.5, 0.16, 0.3))
  mgpVals = c(1.6,0.5,0)

  # Notation:
  # hAAA - h(istrogram) for method AAA.

  # Llin method - plotting binned data on log-linear axes then fitting regression
  #  as done by Daan et al. 2005.
  # hLlin.list = Llin.method(x, num.bins = num.bins)
  hLlin.list = Llin.method.counts(valCounts, num.bins = num.bins)

  #eightMethodsRes = data.frame("Year"=oneYear, "Method"="Llin",
  #    "b" = hLlin.list$slope,
  #    "confMin"= hLlin.list$confVals[1], "confMax"= hLlin.list$confVals[2])
  # That gives hLlin.mids as the name of the row, so does this though
  hLlin.b = hLlin.list$slope
  hLlin.confMin = hLlin.list$confVals[1]
  hLlin.confMax = hLlin.list$confVals[2]
  hLlin.stdErr = hLlin.list$stdErr
  eightMethodsRes = data.frame("Year"=oneYear,
                               "Method"="Llin",
                               "b" = hLlin.b,
                               "confMin"= hLlin.confMin,
                               "confMax"= hLlin.confMax,
                               "stdErr" = hLlin.stdErr,
                               row.names=NULL)

  plot( hLlin.list$mids,
       hLlin.list$log.counts,
       xlab=expression(paste("Bin midpoints for data ", italic(x))),
       ylab = "Log (count)",
       mgp=mgpVals)

  lm.line(hLlin.list$mids, hLlin.list$lm)
  inset = c(0, -0.04)     # inset distance of legend
  legend("topright",
         paste("(a) Llin slope=", signif(hLlin.list$slope, 3)),
         bty="n",
         inset=inset)

  mtext(
   paste("                                                           ",
         oneYear) )

  # LT method - plotting binned data on log-log axes then fitting regression,
  #  as done by Boldt et al. 2005, natural log of counts plotted against natural
  #  log of size-class midpoints.

  # Use Llin method's binning.
  hLT.log.mids = log(hLlin.list$mids)
  hLT.log.counts = log(hLlin.list$counts)
  hLT.log.counts[ is.infinite(hLT.log.counts) ] = NA
                    # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLT.lm = lm( hLT.log.counts ~ hLT.log.mids, na.action=na.omit)
  hLT.slope = hLT.lm$coeff[2]
  hLT.conf = confint(hLT.lm, "hLT.log.mids", 0.95)
  hLT.stdErr = coef(summary(hLT.lm))["hLT.log.mids", "Std. Error"]

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LT"),
      "b" = hLT.slope, "confMin"= hLT.conf[1], "confMax"= hLT.conf[2],
      "stdErr" = hLT.stdErr, row.names=NULL))

  plot( hLT.log.mids, hLT.log.counts,
       xlab=expression(paste("Log (bin midpoints for data ", italic(x), ")")),
       ylab = "Log (count)", mgp=mgpVals)

  lm.line(hLT.log.mids, hLT.lm)
  legend("topright", paste("(b) LT b=", signif(hLT.slope, 3)), bty="n",
         inset=inset)

  # LTplus1 method - plotting linearly binned data on log-log axes then fitting
  #  regression of log10(counts+1) vs log10(midpoint of bins), as done by
  #  Dulvy et al. (2004).

  # Use Llin method's binning.
  hLTplus1.log10.mids = log10(hLlin.list$mids)
  hLTplus1.log10.counts = log10(hLlin.list$counts + 1)
  hLTplus1.log10.counts[ is.infinite(hLTplus1.log10.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
                  #  but the + 1 avoids this issue here
  hLTplus1.lm = lm( hLTplus1.log10.counts ~ hLTplus1.log10.mids,
      na.action=na.omit)
  hLTplus1.slope = hLTplus1.lm$coeff[2]
  hLTplus1.stdErr = coef(summary(hLTplus1.lm))[
      "hLTplus1.log10.mids", "Std. Error"]

  hLTplus1.conf = confint(hLTplus1.lm, "hLTplus1.log10.mids", 0.95)
  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LTplus1"),
      "b" = hLTplus1.slope, "confMin"= hLTplus1.conf[1],
      "confMax"= hLTplus1.conf[2], "stdErr" = hLTplus1.stdErr, row.names=NULL))

  plot( hLTplus1.log10.mids, hLTplus1.log10.counts,
       xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
       ylab = "Log10 (count+1)", mgp=mgpVals)

  lm.line(hLTplus1.log10.mids, hLTplus1.lm)
  legend("topright", paste("(c) LTplus1 b=", signif(hLTplus1.slope, 3)),
       bty="n", inset=inset)

  # LBmiz method

  hLBmiz.num.bins = num.bins

  beta = nlm(LBmizbinsFun, 2, xmin=xmin, xmax=xmax, k=hLBmiz.num.bins)$est

  # hLBmiz.bins = c(beta^(0:(k-1)) * xmin, xmax)
  hLBmiz.bins = c(beta^(0:(hLBmiz.num.bins-1)) * xmin, xmax)
  hLBmiz.mins = hLBmiz.bins[-length(hLBmiz.bins)]    # min of each bin
   # Mizer bin specification, with final bin being same width as penultimate bin

  # Adapting from Llin.method.counts:
  LBmiz.valCounts = dplyr::mutate(valCounts,
                binMin = hLBmiz.mins[findInterval(bodyMass, hLBmiz.bins,
                rightmost.closed=TRUE)])    # Note that this would retain groups
  LBmiz.binVals = dplyr::summarise(dplyr::group_by(LBmiz.valCounts, binMin),
                binCount = sum(Number))
        # No guarantee that every binMid hLlin.mids shows up here (unlikely,
        #  given linear binning and long-tailed data). Even though need to
        #  exclude 0 counts (since get logged) from lm, still best to return
        #  values for all bins, especially since the bins can be specified.
  LBmiz.missing1 = setdiff(hLBmiz.mins, LBmiz.binVals$binMin)
  LBmiz.missing = cbind(binMin = LBmiz.missing1,
      binCount = rep(0, length(LBmiz.missing1)))
  LBmiz.binVals = rbind(LBmiz.binVals, LBmiz.missing)
                    # works even if LBmiz.missing has no rows
  LBmiz.binVals = dplyr::tbl_df(LBmiz.binVals)
  LBmiz.binVals = dplyr::arrange(LBmiz.binVals, binMin)
  if( max( abs( hLBmiz.mins - LBmiz.binVals$binMin) ) > 0)
            { stop("check LBmiz.binVals$binMin in eightMethods.count") }

  hLBmiz.log.min.of.bins = log(LBmiz.binVals$binMin)    # min of bins
  hLBmiz.log.counts = log(LBmiz.binVals$binCount)
  hLBmiz.log.counts[ is.infinite(hLBmiz.log.counts) ] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin

  hLBmiz.lm = lm( hLBmiz.log.counts ~ hLBmiz.log.min.of.bins, na.action=na.omit)

  # hLBmiz.slope = hLBmiz.lm$coeff[2]
  # Need to subtract 1, since want to work in terms of b not slopes now
  hLBmiz.b = hLBmiz.lm$coeff[2] - 1

  hLBmiz.conf = confint(hLBmiz.lm, "hLBmiz.log.min.of.bins", 0.95) - 1
  hLBmiz.stdErr = coef(summary(hLBmiz.lm))[
      "hLBmiz.log.min.of.bins", "Std. Error"]

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LBmiz"),
      "b" = hLBmiz.b, "confMin"= hLBmiz.conf[1],
      "confMax"= hLBmiz.conf[2], "stdErr" = hLBmiz.stdErr, row.names=NULL))

  plot( hLBmiz.log.min.of.bins, hLBmiz.log.counts,
     xlab=expression(paste("Log (minima of bins for data ", italic(x), ")")),
     ylab = "Log (count)", mgp=mgpVals)
     # axes=FALSE, xaxs="i", yaxs="i", , xlim=c(log10(1), log10(650)),
     #  ylim=c(log10(0.7), log10(1100)))  # So axes are logged

  lm.line(hLBmiz.log.min.of.bins, hLBmiz.lm)
  legend("bottomleft", paste("(d) LBmiz b=", signif(hLBmiz.b, 3)),
         bty="n", inset=inset)

  # LBbiom method - binning data using log2 bins, calculating biomass not counts
  #  in each bin, plotting log10(biomass in bin) vs log10(midpoint of bin)
  #  as done by Jennings et al. (2007), who used bins defined by a log2 scale.

  hLBNbiom.list = LBNbiom.method.counts(valCounts)    # Does this method and the next.

  hLBbiom.b = hLBNbiom.list[["unNorm.slope"]] - 2
  hLBbiom.conf = hLBNbiom.list[["unNorm.conf"]] - 2
  hLBbiom.stdErr = coef(summary(hLBNbiom.list$unNorm.lm))[
      "log10binMid", "Std. Error"]

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LBbiom"),
      "b" = hLBbiom.b, "confMin"= hLBbiom.conf[1],
      "confMax"= hLBbiom.conf[2], "stdErr" = hLBbiom.stdErr, row.names=NULL))

  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiom,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["unNorm.lm"]])
  legend("bottomleft", paste("(e) LBbiom b=",
     signif(hLBbiom.b, 3)), bty="n", inset=c(-0.08, -0.04))

  # LBNbiom method - on biomass, not counts, as per Julia Blanchard's 2005 paper.
  #  log2 bins of bodymass, sum the total biomass in each bin, normalise
  #  biomasses by binwidths, fit regression to log10(normalised biomass) v
  #  log10(midpoint of bin).

  hLBNbiom.b = hLBNbiom.list[["norm.slope"]] - 1
  hLBNbiom.conf = hLBNbiom.list[["norm.conf"]] - 1
  hLBNbiom.stdErr = coef(summary(hLBNbiom.list$norm.lm))[
      "log10binMid", "Std. Error"]

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LBNbiom"),
      "b" = hLBNbiom.b, "confMin"= hLBNbiom.conf[1],
      "confMax"= hLBNbiom.conf[2], "stdErr" = hLBNbiom.stdErr, row.names=NULL))

  plot(hLBNbiom.list[["binVals"]]$log10binMid,
     hLBNbiom.list[["binVals"]]$log10totalBiomNorm,
     xlab=expression(paste("Log10 (bin midpoints for data ", italic(x), ")")),
     ylab = "Log10 (normalised biomass)", mgp=mgpVals)

  lm.line(hLBNbiom.list[["binVals"]]$log10binMid, hLBNbiom.list[["norm.lm"]])
  legend("bottomleft", paste("(f) LBNbiom b=",
     signif(hLBNbiom.b, 3)), bty="n", inset=inset)

  # Cumulative Distribution, LCD method
  # logProp = log((1:length(x))/length(x))               # x equivalent:
  # This method should really split up the cumProp into, say, 1000 values
  #  since for the regression each point gets weighted the same but this
  #  isn't quite right. But this is how someone would likely do it.
  #  To plot the results for the MLE method below I'm generating 1000 values.
  LCD.valCounts = dplyr::select(valCounts, bodyMass, Number)
  LCD.valCounts = dplyr::arrange(LCD.valCounts, desc(bodyMass))
                                        # x.sorted = sort(x, decreasing=TRUE)
  sumNumber = sum(LCD.valCounts$Number)
  LCD.valCounts = dplyr::mutate(LCD.valCounts, cumSum = cumsum(Number))
                                        # 1:length(x)
  LCD.valCounts = dplyr::mutate(LCD.valCounts, cumProp = cumSum / sumNumber)
                                        # logProp = log((1:length(x))/length(x))

  LCD.valCounts = dplyr::mutate(LCD.valCounts, logBodyMass = log(bodyMass),
                            logCumProp = log(cumProp))
                                       # logSorted = log(x.sorted)
  hLCD.lm = lm(logCumProp ~ logBodyMass, data = LCD.valCounts)
                                        # hLCD.lm = lm(logProp ~ logSorted)
  hLCD.slope = hLCD.lm$coeff[2]

  hLCD.b = hLCD.lm$coeff[2] - 1
  hLCD.conf = confint(hLCD.lm, "logBodyMass", 0.95) - 1
  hLCD.stdErr = coef(summary(hLCD.lm))[
      "logBodyMass", "Std. Error"]

  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("LCD"),
      "b" = hLCD.b, "confMin"= hLCD.conf[1],
      "confMax"= hLCD.conf[2], "stdErr" = hLCD.stdErr, row.names=NULL))

   plot(LCD.valCounts$logBodyMass, LCD.valCounts$logCumProp, main="",
     xlab=expression(paste("Log ", italic(x))),
     ylab=expression( paste("Log (prop. of ", values >= italic(x), ")")),
     mgp=mgpVals)

  lm.line(LCD.valCounts$logBodyMass, hLCD.lm, col="red")

  legend("bottomleft", paste("(g) LCD b=", signif(hLCD.b, 3)), bty="n",
       inset=inset)


  # MLE (maximum likelihood method) calculations.
  MLE.valCounts = dplyr::select(valCounts, bodyMass, Number)
  # Should be faster to group repeated values, just in case that's not done:
  MLE.valCounts = dplyr::summarise(dplyr::group_by(valCounts, bodyMass), Count = sum(Number))
  MLE.valCounts = dplyr::arrange(MLE.valCounts, desc(bodyMass))

  sumCounts = sum(MLE.valCounts$Count)
  if(abs( sumCounts - sumNumber) > 0.001)
      { stop("Check sumCounts in eightMethods.count()") }
  MLE.K = dim(MLE.valCounts)[1]         # Number of bodyMass values
  if(MLE.valCounts[1, "Count"] == 0 |
     MLE.valCounts[MLE.K, "Count"] == 0)
      { stop("Need first and last counts to be zero in
          eightMethods.count() for MLE method")}

  # Can adapt the code, but better to just to take such zero counts out of data
  # Adapting (for counts) analytical value of MLE b for PL model
  #  (Box 1, Edwards et al. 2007)
  #  as a starting point for nlm for MLE of b for PLB model.
  MLE.xmin = min(MLE.valCounts$bodyMass)
  MLE.xmax = max(MLE.valCounts$bodyMass)
  MLE.sumCntLogMass = sum(MLE.valCounts$Count * log(MLE.valCounts$bodyMass))
  PL.bMLE.counts = 1/( log(MLE.xmin) - MLE.sumCntLogMass/sumCounts) - 1

  PLB.minLL =  nlm(negLL.PLB.counts, p=PL.bMLE.counts,
      x=MLE.valCounts$bodyMass, c=MLE.valCounts$Count, K=MLE.K,
      xmin=MLE.xmin, xmax=MLE.xmax, sumclogx=MLE.sumCntLogMass)
                                        #, print.level=2 )

  PLB.bMLE = PLB.minLL$estimate

  # 95% confidence intervals for MLE method.

  PLB.minNegLL = PLB.minLL$minimum

  # Values of b to test to obtain confidence interval. For the real movement data
  #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  #  symmetric interval here.

  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.00001)

  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
      {
          PLB.LLvals[i] = negLL.PLB.counts(bvec[i], x=MLE.valCounts$bodyMass,
           c=MLE.valCounts$Count, K=MLE.K,
           xmin=MLE.xmin, xmax=MLE.xmax, sumclogx=MLE.sumCntLogMass)
      }
  critVal = PLB.minNegLL  + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
                      # b values in 95% confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if(PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == max(bvec))
    { windows()
      plot(bvec, PLB.LLvals)
      abline(h = critVal, col="red")
      stop("Need to make bvec larger - see R window")   # Could automate
    }

  #  Assume for now that confidence interval is PLB.bMLE +/- 1.96 * stdErr
  #   just to quickly calc stdErr. So stdErr = -(confMin - PLB.bMLE)/1.96.
  #      PLB.bMLE - PLB.MLE.bConf - is symmetric anyway, implying quadratic
  #      likelihood profile, implying normal approximation is okay. So use
  #      it now to go backwards. To properly calculate should do the
  #      Fisher information. stdErr = 1 / (sqrt( d^2 logLik /db^2 at MLE))
  #      though that is fiddly to derive, where d is derivative.
  PLB.MLE.stdErr = mean(abs((PLB.MLE.bConf - PLB.bMLE)/1.96))

  # MLE.rep.xmax[iii] = xmax   - not storing xmax for now
  eightMethodsRes = rbind(eightMethodsRes,
      data.frame("Year"=oneYear, "Method"=as.factor("MLE"),
      "b" = PLB.bMLE, "confMin"= PLB.MLE.bConf[1],
      "confMax"= PLB.MLE.bConf[2], "stdErr" = PLB.MLE.stdErr, row.names=NULL))
  # To plot rank/frequency style plot, because of non-integer counts want to
  #  calculate cumulative proportions (as in LCD), and then generate 1,000
  #  individual body masses values that 'represent' the observed distribution,
  #  and then plot y axis as Proportion >= x, with 1,000 points plotted.
  MLE.valCounts = dplyr::mutate(MLE.valCounts, cumSum = cumsum(Count))
                                        # x equivalent:  1:length(x)
  MLE.valCounts = dplyr::mutate(MLE.valCounts, cumProp = cumSum / sumCounts)
                                        # logProp = log((1:length(x))/length(x))
  # So if you had a sample of 1,000 that includes xmin and xmax, what would
  #  the remaining 998 body masses be?:
  if(MLE.valCounts[dim(MLE.valCounts)[1],"cumProp"] != 1)
      { stop("Check MLE.valcounts in eightMethods.count()")  }
  MLE.sim = dplyr::tbl_df(data.frame(cumPropSim =
      seq(as.numeric(MLE.valCounts[1,"cumProp"]), 1, length=ceiling(sumCounts))))
                                        # simulated cumulative proportions
                                        # or maybe do on sumSum, though props
                                        #  better capture endpoints I think.
                                        # length does ceiling anyway, so make
                                        #  explicit here. Did try just 1000
                                        #  example fish, but that doesn't capture
                                        #  the details in the tail (since
                                        #  actual sample size is 33,592.02), so
                                        #  this does what a sample of 33,593
                                        #  would likely look like.
  # Had to update this due to change in dplyr; this approach wouldn't work, or
  #  using select
  #MLE.sim = dplyr::mutate(MLE.sim,
  #                 bodyMassSim = as.numeric(MLE.valCounts[findInterval(
  #                                  cumPropSim, MLE.valCounts$cumProp),
  #                                 "bodyMass"]))

  MLE.sim = dplyr::mutate(MLE.sim,
                          bodyMassSim = MLE.valCounts[findInterval(cumPropSim,
                                                                   MLE.valCounts$cumProp), ]$bodyMass)

      # findInterval():
      #  vec = MLE.valCounts$cumProp     # breakpoints
      # Find interval containing each of the  elements of
      #  x = MLE.sim$cumPropSim  (x in findInterval() terminology)

      #        LCD.valCounts = dplyr::mutate(LCD.valCounts, logBodyMass = log(bodyMass),
      #                      logCumProp = log(cumProp))
                                       # logSorted = log(x.sorted)
  plot(MLE.sim$bodyMassSim, MLE.sim$cumPropSim, log="xy",
       xlab=expression(paste("Values, ", italic(x))),
       ylab=expression( paste("Proportion of ", values >= x)), mgp=mgpVals)

  x.PLB = seq(xmin, xmax, length=1000)     # x values to plot PLB
  y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
                    xmax = max(x.PLB))) # * sumCounts
  lines(x.PLB, y.PLB, col="red") #, lty=5)

  legend("bottomleft", paste("(h) MLE b=", signif(PLB.bMLE, 3)), bty="n",
         inset=inset)

  # To add the curves at the limits of the 95% confidence interval:
  #for(i in c(1, length(bIn95)))   # for(i in 1:length(bIn95))  to see all vals
  #    {
  #      lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
  #                  xmax = max(x.PLB))) * length(x), col="red", lty=3)
  #    }

  dev.off()
  return(eightMethodsRes)
}

##' Construct bins either double in size or are of equal width, and encompass
##'  the data
##'
##' Construct bins that start from `floor(min(x))` or `min(x)` and either double
##'    in size or are of equal width, and encompass the data. More generalised
##'    version of `log2bins()`.
##'
##' @param x vector of individual values (e.g. body masses). Only `x` or
##'   `counts` can be specified.
##' @param counts dataframe (or array) with first column being an x value
##'  (e.g. body mass), and second column being the counts of the
##'  number of individuals for that value. Only `x` (the vector) or `counts` can
##'   be specified.
##' @param binWidth type of bins to use:
##'   * `"2k"` will result in `binBreaks` that:
##'     + with `startInteger=TRUE` are powers of 2, i.e. ..., 0.25, 0.5, 1, 2, 4, 8, 16,....
##'     + with `startInteger=FALSE` are bins that double in size and  start with
##'       `min(x)`; not yet implemented, since have to think about what the width of
##'       the first bin should be.
##'   * numeric value (call it `a`) will result in binBreaks are separated by `a` and span the
##'       data, that:
##'     + with `startInteger=TRUE` start from `z = floor(min(x))` and are then
##'          `z, z+a, z+2a, z+3a, ....`   (if `z = 0` then power-law cannot be fit
##'        so then need to use `startInteger=FALSE`)
##'     + with `startInteger=FALSE` start from `z = min(x)` and are then
##'           `z, z+a, z+2a, z+3a, ....`
##'   * only `binWidth` or `binBreaks` can be specified.
##' @param binBreaks pre-defined bin breaks as a vector. Only `binWidth`
##'   or `binBreaks` can be specified.
##' @param startInteger TRUE or FALSE, whether to start the bin breaks at an integer
##'   power of 2 (for method `"2k"`) or an integer. See `binWidth` above.
##'   `startInteger` is ignored if `binBreaks` is specified.
##' @return list containing:
##'   * indiv: dataframe with a row for each `x` value, with columns:
##'      + `x`: original `x` values
##'      + `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum,
##'      maximum, and width, respectively, of the bin within
##'      which the `x` value falls.  If indiv has `>=10^6` rows then it isn't saved.
##'      If `counts` was specified then an equivalent `x`
##'      vector is created and is column `x` (i.e. `x` values are repeated). May
##'      not be the most efficient way, but it easiest to program.
##'   * binVals: dataframe with a row for each new bin and columns:
##'      + `binMid`, `binMin`, `binMax`, `binWidth`: midpoint, minimum,
##'         maximum, and width, respectively, of the bin
##'      + `binCount`: total number of individuals in that bin
##'      + `binCountNorm`: normalised bin count, `binCount / binWidth`
##'      + `binSum`: sum of individual values in that bin (appropriate if `x`
##'         represents biomass, but not length)
##'      + `binSumNorm`: `binSum / binWidth`
##'      + `log10....` - `log10()` of some of the above quantities
##' @export
##' @author Andrew Edwards
binData = function(x = NULL, counts = NULL, binWidth = NULL, binBreaks = NULL,
                   startInteger = TRUE)
    {
        if(!is.null(x) & !is.null(counts)) {
            stop("need only one of x or counts in binData") }
        if(is.null(x) & is.null(counts)) {
            stop("need x or counts in binData") }
        if(!is.null(x)) {
          if(!is.vector(x))stop("x not a vector in binData")
          if(anyNA(x)) stop("x contains NA's in binData")
          if(min(x) <= 0)stop("x needs to be >0 in binData")
          }
        if(!is.null(counts))  {
          if(dim(counts)[2] != 2)stop("counts needs two cols in binData")
          if(min(counts[,1]) < 0) {
              stop("x values in counts need to be >= 0 in binData") }
          if(min(counts[,2]) < 0) {
              stop("numbers in counts need to be >= 0 in binData") }
          if(sum(!is.wholenumber(counts[,2])) > 0) {
              stop("numbers in counts need to be integers in binData;
                      for non-integer count see a new function. Currently,
                      such a new function has no name [so says Jaqen H'ghar]. Or it may be easier to
                      adapt binData.") }
          }
        if(is.null(binWidth) & is.null(binBreaks)) {
            stop("need one of binWidth or binBreaks in binData") }
        if(!is.null(binWidth) & !is.null(binBreaks)) {
            stop("need only one of binWidth or binBreaks in binData") }
        if(startInteger != "TRUE" & startInteger != "FALSE"){
            stop("startInteger must be TRUE or FALSE in binData") }
        # As for LBNbiom.method(), could write code that would make
        #  use of the counts dataframe explicitly, but actually quite easy
        #  to just create the longer vector x (though may be slightly slower
        #  computationally), to save writing extensive new code. But do this
        #  at some point for noninteger counts.
        if(is.null(x))
           { x = rep(counts[,1], counts[,2])
             minx = min(counts[,1])
             maxx = max(counts[,1])
           } else
           { minx = min(x)
             maxx = max(x)
           }
        #
        if(!is.null(binBreaks))
           {
           if(minx < min(binBreaks) | maxx > max(binBreaks) )
             { stop("binBreaks do not span data in binData")
             }
           } else           # create binBreaks based on binWidth
           {
           if(binWidth == "2k")
             {
             if(startInteger)
               { binBreaks = 2^( floor(log2(minx)) : ceiling(log2(maxx)) )
               } else
               { stop("startInteger currently needs to be TRUE when
                   binWidth = 2k")
               }
             } else     # If not "2k"
             {
             if(!is.numeric(binWidth))
               { stop("binWidth must be 2k or a number (in quotes is okay
                         in quotes) in binData().")
               }
             # startInteger says whether to start from an integer value or
             #  start from min(x),
             z = floor(minx) * startInteger + minx * !startInteger
             binBreaks = seq( z, by=binWidth,
                        length = ceiling( (maxx - z)/binWidth) + 1)
             }
           }

        indiv = data.frame(x)       # dataframe with one row for each individual
        indiv$binMid =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)] + 0.5*diff(binBreaks))
        indiv$binMin =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-length(binBreaks)])
        indiv$binMax =cut(x, breaks=binBreaks, right=FALSE, include.lowest=TRUE,
            labels = binBreaks[-1])
        #
        indiv$binMid = as.numeric(as.character(indiv$binMid))
        indiv$binMin = as.numeric(as.character(indiv$binMin))
        indiv$binMax = as.numeric(as.character(indiv$binMax))
        # Now calculate biomass in each bin class:
        binVals = dplyr::summarise(dplyr::group_by(indiv, binMid),
            binMin = unique(binMin),
            binMax = unique(binMax),
            binWidth = binMax - binMin,
            binCount = length(x),
            binCountNorm = binCount / binWidth,
            binSum = sum(x),
            binSumNorm = binSum / binWidth )
           # binWidth uses new columns binMax and binMin
        # Indices for minima of bins that have zero counts and so do not
        #  appear in binVals yet:
        emptyBinMinInd = !(binBreaks[-length(binBreaks)] %in% binVals$binMin)
        emptyBinMin = binBreaks[emptyBinMinInd]
        empties = length(emptyBinMin)
        emptyBinMax = binBreaks[-1][emptyBinMinInd]
        emptyBinWidth = emptyBinMax - emptyBinMin
        emptyBinMid = emptyBinMin + emptyBinWidth/2

        emptyVals = as.data.frame(cbind(emptyBinMid,
                                        emptyBinMin,
                                        emptyBinMax,
                                        emptyBinWidth,
                                        matrix(0, nrow=empties, ncol=4)))
        names(emptyVals) = names(binVals)
        binVals = rbind(binVals, emptyVals)         # still a local df

        binVals = binVals[order(binVals$binMid), ]   # order by binMid

        binVals = dplyr::mutate(binVals,
                                log10binMid = log10(binMid),
                                log10binCount = log10(binCount),
                                log10binSum = log10(binSum),
                                log10binSumNorm = log10(binSumNorm))
        binVals[is.infinite(binVals$log10binCount),
                "log10binCount"] = NA
                  # lm can't cope with -Inf, which appear if 0 counts in a bin
        binVals[is.infinite(binVals$log10binCountNorm),
                "log10binCountNorm"] = NA
        binVals[is.infinite(binVals$log10binSum),
                "log10binSum"] = NA
        binVals[is.infinite(binVals$log10binSumNorm),
                "log10binSumNorm"] = NA
        if(dim(indiv)[1] < 10^6) {       # only save indiv if not too big
          y = list(indiv = indiv, binVals = binVals)
          } else
          {
          y = list(binVals = binVals)
          }
        return(y)
       }
