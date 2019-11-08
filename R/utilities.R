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

##' Format x to have the specified number of decimal places
##'
##' TODO: change points to places presumably
##' TODO: give a more descriptive name
##' Format x to have supplied number of decimal places, have
##'  thousands seperated by commas.
##' Taken from our Pacific Hake assessment:
##' https://github.com/cgrandin/hake-assessment/blob/master/doc/r/r-functions/utilities.r
##'
##' @param x TODO
##' @param dec.points  TODO
##' @return x with the specified number of decimal places, and comma for thousands
##' @export
##' @author Chris Grandin
f <- function(x, dec.points = 0){
  return(format(round(x,dec.points), big.mark = ",", nsmall = dec.points))
}
