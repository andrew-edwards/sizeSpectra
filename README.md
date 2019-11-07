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




## Temp, keep for now, from **fitting-size-spectra**:

This is a **branch** of the **fitting-size-spectra** repository that contains the new code for our second manuscript (concerning dealing with data that are only available in binned form) in the directory **codeBinStructure/**. The branch will be merged into **master** at some point, I just wanted to keep it separate for now. If you don't understand that GitHub terminology do not worry, just see the file 
**codeBinStructure/readMeBinStructure.txt** to get started. 

Below is the original code repository for the ***Methods in Ecology and Evolution*** paper. The code for the new manuscript is in the directory **codeBinStructure/**, but relies on some of the functions from the first paper, which is why I have put the new code into the same repository.

To download the repository from the GitHub site just click the 'Clone or Download' button (near the top on the right) and select 'Download ZIP' (which is how I will create the **.zip** file to submit with the manuscript).

[Note to self: for submitting manuscript, download the zipped version, then delete all the **code/** files except for **code/PLBfunctions.r** that are referred to from the new code. And replace this **readme.md** with a simple **readme.txt** that mentions **code/PLBfunctions.r** and points the user to 
**codeBinStructure/readMeBinStructure.txt**.

---


The aim of sharing this code is so that others can repeat (and extend) our simulation study, and also analyse their own data. I have endeavoured to properly document the code.

To download the code from the GitHub site just click the 'Clone or Download' button (near the top on the right) and select 'Download ZIP'. Make a note of the 'Latest commit' number in case you have any questions for me. If you use GitHub then feel free to fork and adapt the code.

If you want to run new simulations or apply the code to your own data then just use the latest version that will be automatically displayed on the GitHub site.

To *exactly* reproduce the results in the paper
you should download release version 1.0.0 (click on 'release' tab in the GitHub site) and see the notes in **code/readReRun.txt** concerning the seed and how **require(dplyr)** generates a random number. Later updates of the code will have some generalisation of functions that should not affect the older code (but I won't re-test it all for back compatibility). 

There are functions (in **code/PLBfunctions.r**) that may be of more general use, such as **logTicks()** for adding tick marks to a log-log plot, and **legJust()** for right-justifying a legend (based on an example in **?legend**). 

If you have problems with the code then please contact me. Some of it has been independently used by co-author James Robinson, who got it working fine for his own data.

Thanks,

Andrew Edwards. 

<http://www.chebucto.ns.ca/~english>

Andrew.Edwards@dfo-mpo.gc.ca - this may change at some point, and if it does not work then try andrew.edwards.dfo@gmail.com which will always automatically forward to my correct email.

# Repository Contents

**README.md** - this file (can be read in a Markdown viewer or simply any text editor).

**.gitignore** - ignore this if you don't know what it is.

**code/** - directory containing all the R code.

**code/readMeCode.txt** - readme file for the code directory. Contains instructions and details of the files and further subdirectories.

The subdirectories of **code/** are summarised below, but see **readMeCode.txt** for full details.

**code/single/** - testing the eight methods on a single data set.

**code/multiple/** - testing the eight methods on multiple (10,000) simulated data sets. Contains subdirectories for the sensitivity tests.

**code/MLEbin/** - MLEbin method for likelihood when the data are only available in binned form.

**code/recommend/** - recommended likelihood calculations and resulting plots of data and fitted size spectrum (Figure 6).

**codeBinStructure/** - directory containing code for second manuscript, see **readMeBinStructure.txt** for details.

# Updates for Methods in Ecology and Evolution code
 
15th May 2017 - corrected a minor error regarding the likelihood function for the MLEbin method when b=-1. In practice this should not have any effect. See Issue #7.

12th December 2017 - fully corrected the above point. Phil Wallhead correctly pointed out that I had not completely corrected the code for b=-1. Again, in practice this should not have any effect (though it would for simulating data with b=-1). See Issue #7.


