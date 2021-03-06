% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{pBiomass}
\alias{pBiomass}
\title{Biomass distribution function from MEE equations A.4 and A.8}
\usage{
pBiomass(x, b, xmin, xmax, n)
}
\arguments{
\item{x}{vector of values for which to calculate the total biomass between
\code{xmin} and the value}

\item{b}{exponent of the PLB distribution}

\item{xmin}{minimum bound of the distribution, \code{xmin > 0}}

\item{xmax}{maximum bound for bounded distribution, \code{xmax > xmin}}
}
\value{
returns of vector total biomass between \code{xmin} and each value of \code{x}
}
\description{
Total biomass between \code{xmin} and \code{x}, assuming a bounded power-law
distribution of body masses between \code{xmin} and \code{xmax} and a given value of
exponent \code{b}, and a total of \code{n} individuals.
Given by MEE equations A.4 and A.8. Can then be called by \code{pBiomassBins} to
give total biomass (and normalised biomass) in each bin.
}
\author{
Andrew Edwards
}
