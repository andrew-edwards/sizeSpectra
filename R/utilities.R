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
