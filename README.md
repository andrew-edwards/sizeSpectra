# sizeSpectra
R package for fitting size spectra to ecological data (including binned data)

 <!-- badges: start -->
 [![R-CMD-check](https://github.com/andrew-edwards/sizeSpectra/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/andrew-edwards/sizeSpectra/actions/workflows/R-CMD-check.yaml)
 <!-- [![Codecov test coverage](https://codecov.io/gh/andrew-edwards/sizeSpectra/branch/master/graph/badge.svg)](https://codecov.io/gh/andrew-edwards/sizeSpectra?branch=master) -->
  <!-- badges: end -->

## Description
This R package contains functions for fitting size spectra to ecological data. In particular, it contains functionalised code to reproduce all the results in [1] and [2], and for users to apply the methods to their own data. See the [vignettes](http://htmlpreview.github.io/?https://github.com/andrew-edwards/sizeSpectra/blob/master/doc/vignettes_overview.html) to see what the package can do.

[1] **Testing and recommending methods for fitting size spectra to data** by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum and Julia L. Blanchard. ***Methods in Ecology and Evolution*** (2017, 8:57-67). Freely available at <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12641/full>.

[2] **Accounting for the bin structure of data removes bias when fitting size spectra** by Andrew M. Edwards, James P. W. Robinson, Julia L. Blanchard, Julia K. Baum and Michael J. Plank. ***Marine Ecology Progress Series*** (2020, 636:19-33). Freely available at <https://www.int-res.com/abstracts/meps/v636/p19-33/>.

The size spectrum of an ecological community characterizes how a property, such as abundance or biomass, varies with body size. Size spectra are often used as ecosystem indicators of marine systems. Past applications have included groundfish trawl surveys, visual surveys of fish in kelp forests and coral reefs, sediment samples of benthic invertebrates and satellite remote sensing of chlorophyll, as well as terrestrial systems. Various methods have been used to fit size spectra over the past decades, and in [1] we tested eight of them and recommend the use of maximum likelihood. In [2] we extended the likelihood method to properly account for the bin structure of data.

Here is a movie showing fits to 30 years of International Bottom Trawl Survey data fit and plotted using our new methods, where uncertainty due to the bin structure is properly accounted for (see [2] for full explanations):

![IBTS 30 years of data](vignettes/IBTS_movie.gif)

## Vignettes

See the [overview vignette](http://htmlpreview.github.io/?https://github.com/andrew-edwards/sizeSpectra/blob/master/doc/vignettes_overview.html) for a summary of all vignettes. They are available on this GitHub site, and also within the package in the usual way (see below).

The vignettes explain how to use the functions in the package to reproduce all results in both papers, and to analyse new data sets using our functions. The vignettes are descriptions of how to use the code to implement the methods. The two papers should be consulted first to understand the methods (I have tried to avoid repeating text from the papers in the vignettes). In the vignettes, [1] is referred to as the 'MEE paper' and [2] as the 'MEPS paper'. 

**In particular** see the new `MLEbin_recommend.Rmd` which starts with the core functions needed to fit and plot data using the MLEbin method. Use this as a template if you want to fit your own data. 

## Main updates since version 1.0.0.0 (released Dec 2019)

New vignette (`MLEbin_recommend.Rmd`), linked in the [overview vignette](http://htmlpreview.github.io/?https://github.com/andrew-edwards/sizeSpectra/blob/master/doc/vignettes_overview.html), that gives the core steps for using the MLEbin method, plus also explores some new plotting approaches for binned data using the new functions:
 - `ISD_bin_plot_nonoverlapping()` -- Recommended plotting for binned data with non-overlapping bins, which is the usual case.
 - `LBN_bin_plot()` -- Biomass size spectrum plot for binned data demonstrating uncertainties, showing the bin widths explicitly and the normalised biomass in each bin (with resulting uncertainties). So extending MEE Fig. 6 for already binned data, using a new approach motivated by MEPS Fig. 7.
 - `plot_binned_fitted()` -- Add horizontal bars and shaded rectangles to `LBN_bin_plot()`, to show estimated normalised biomasses in each bin.


Some other new functions have been added to analyse data in segments, but are not yet
part of vignettes. These include: 
 - `pBiomass()` -- function for the biomass distribution function from equations A.4 and A.8 of the MEE paper.
 - `pBiomassBins()` -- to give total and normalised biomass in each bin for a fitted distribution and given bin breaks.
 - `pBiomassBinsConfs()` -- to call `pBiomassBins()` for three values of b (MLE and confidence interval values).
 - `calcLikeSegments()` -- Calculate MLEs and 95% confidence intervals of `b` for analysing separate body-mass segments of a data set


## Install instructions

To install the latest version of `sizeSpectra` directly from GitHub you need the package `remotes`, so if you don't have it install it (once):

```
install.packages("remotes")    # If you do not already have the "remotes" package
remotes::install_github("andrew-edwards/sizeSpectra")
```

If you get an errors, then try:
```
remotes::install_github("andrew-edwards/sizeSpectra", build_vignettes = FALSE)
```
as there's something in the vignettes that is now not working to do with recent versions of R or of packages (as of 9th May 2023); you can still see the built vignettes as described above. Note that I've had to comment out some of the code to build one of the figures in the MEPS_IBTS_recommend vignette -- something has changed between R versions 4.2.1 and 4.3.0 that stops it building. However, the figure (Figure S.4 of MEPS paper) is likely not of major interest for users to produce. As mentioned, the rendered vignette above has been left as the earlier working version. I build the the original vignette locally, just not when doing `install(dependencies = FALSE, build_vignettes = TRUE)`. Will look into it further (should make an Issue).

Then:
```
library(sizeSpectra)
browseVignettes("sizeSpectra")
```
to list the vignettes (and links to their Rmarkdown, R and html versions). `Vignettes_overview` gives an overview of all vignettes, and is available either from the previous command or by just running
```
vignette("Vignettes_overview")
```
[Note that if you are using Rstudio there is a [known issue](https://github.com/rstudio/rstudio/issues/2253) that equations don't render properly in the viewer; they are fine in a usual html viewer. Also `browseVignettes("sizeSpectra")` didn't seem to work in Rstudio either. Running from a different R console (e.g. the default non-Rstudio version) should work fine.]

## Issues, problems

Please report any problems as a [GitHub Issue](https://github.com/andrew-edwards/sizeSpectra/issues), using a minimal working example if possible (and please check the closed issues first).

Note that this code was written over several years, and then converted into a proper R package. As such, I have not used consistent naming conventions (e.g. some column names are camelCase while others are not) like I would now if starting a package from scratch -- these may be partly corrected if time permits.