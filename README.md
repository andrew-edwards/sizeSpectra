# sizeSpectra
R package for fitting size spectra to ecological data (including binned data)

 <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/andrew-edwards/sizeSpectra.svg?branch=master)](https://travis-ci.org/andrew-edwards/sizeSpectra)
  [![Codecov test coverage](https://codecov.io/gh/andrew-edwards/sizeSpectra/branch/master/graph/badge.svg)](https://codecov.io/gh/andrew-edwards/sizeSpectra?branch=master)
  <!-- badges: end -->

##In development - not usable yet

This R package contains functions for fitting size spectra to ecological data. In particular, it contains functionalised code to reproduce all the results in [1] and [2], and for users to apply the methods to their own data.

[1] **Testing and recommending methods for fitting size spectra to data** by Andrew M. Edwards, James P. W. Robinson, Michael J. Plank, Julia K. Baum and Julia L. Blanchard. ***Methods in Ecology and Evolution*** (2017, 8:57-67). Freely available at <http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12641/full>

[2] ....

The size spectrum of an ecological community characterizes how a property, such as abundance or biomass, varies with body size. Size spectra are often used as ecosystem indicators of marine systems. Past applications have included groundfish trawl surveys, visual surveys of fish in kelp forests and coral reefs, sediment samples of benthic invertebrates and satellite remote sensing of chlorophyll, as well as terrestrial systems. Various methods have been used to fit size spectra over the past decades, and in [1] we tested eight of them and recommend the use of maximum likelihood. In [2] we extended the likelihood method to properly account for the bin structure of data.

TODO: include some figures

The vignettes (TODO include links here, do an overview one) explain how to use the functions in the package to reproduce all results and to analyse new data sets. The vignettes are for implementing the methods (and the two papers should be consulted for full descriptions of all methods - I have tried to avoid repeat the papers in the vignettes; note that [1] is referred to as the 'MEE paper' and [2] as the 'MEPS paper').

TODO: Install instructions

TODO: Vignette instructions

TODO: LInk to Issues
