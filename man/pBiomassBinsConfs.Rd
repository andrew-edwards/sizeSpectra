% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{pBiomassBinsConfs}
\alias{pBiomassBinsConfs}
\title{Wrapper to call \code{pBiomassBins()} for three values of b (MLE and conf limits)}
\usage{
pBiomassBinsConfs(..., binValsTibble, b.MLE, b.confMin, b.confMax)
}
\arguments{
\item{...}{extra arguments passed to \code{bBiomassBins()}: \code{b}, \code{xmin}, \verb{xmax, }n`.}

\item{binValsTibble}{tibble of binned data in the form required for \code{pBiomassBins()}}

\item{b.MLE}{maximum likelihood estimate of \emph{b} (ideally from the MLEbin method)}

\item{b.confMin}{lower 95\\% confidence limits of \emph{b}}

\item{b.confMax}{upper 95\\% confidence limits of \emph{b}}
}
\value{
\code{binValsTibble} with extra columns \code{estBiomassMLE} and
\code{estBiomassNormMLE} for the estimated biomass and normalised biomass for
\code{b.MLE}, extra columns \code{estBiomassConfMin} and \code{estBiomassNormConfMin} for
the same but using \code{b.confMin}, \code{estBiomassConfMax} and
\code{estBiomassNormConfMax} for \code{b.confMax}.
}
\description{
Only takes a tibble as the input data (unlike \code{pBiomassBins()} and
\code{pBiomass()}. Currently needs \code{wmin}, \code{wmax} and \code{Number} as columns.
}
\author{
Andrew Edwards
}
