# utility functions

##' Convert a vector of body lengths to a vector of body masses using length-weight
##' relationship
##'
##' Convert a vector of body lengths to a vector of body masses using the formula
##'  `mass = LWa * lengths^LWb`. `LWa` and `LWb` will be species-specific.
##'  User must ensure that units are consistent.
##' @param lengths vector of lengths
##' @param LWa multiplicative coefficient in equation
##' @param LWb exponent in equation
##' @return vector of body masses
##' @export
##' @author Andrew Edwards
lengthToMass = function(lengths, LWa, LWb)
    {
        if(min(c(lengths, LWa, LWb)) < 0)
            { stop("Need positive arguments in lengthToMass")  }
        return(LWa * lengths^LWb)
    }

##' Convert a vector of body masses to a vector of body lengths using length-weight
##' relationship (for simulations)
##'
##' Convert a vector of body masses to a vector of body lengths using the
##' formula
##'  `length = (mass / LWa)^(1/b)`, which is the inverse of the standard
##'  `mass = LWa * lengths^LWb`. `LWa` and `LWb` will be species-specific.
##'  User must ensure that units are consistent.
##' @param masses vector of masses
##' @param LWa multiplicative coefficient in equation
##' @param LWb exponent in equation
##' @return vector of body lengths
##' @export
##' @author Andrew Edwards
massToLength = function(masses, LWa, LWb)
    {
      if(min(c(masses, LWa, LWb)) < 0){
        stop("Need positive arguments in massToLength")}
      return( (masses/LWa)^(1/LWb) )
    }


##' Format x to have the specified number of decimal places
##'
##' Format x to have supplied number of decimal places, have
##'  thousands seperated by commas.
##' Taken from our Pacific Hake assessment:
##' https://github.com/cgrandin/hake-assessment/blob/master/doc/r/r-functions/utilities.r
##'
##' @param x valuesTODO
##' @param dec.points number of decimal places
##' @return x with the specified number of decimal places, and comma for thousands
##' @export
##' @author Chris Grandin
f <- function(x, dec.points = 0){
  return(format(round(x,dec.points), big.mark = ",", nsmall = dec.points))
}
