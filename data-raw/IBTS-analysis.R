# The main analyses of the IBTS data in the MEPS paper.

load("nSeaFungAnalysis.RData")    # Output from nSeaFungAnalysis.Snw, TODO
                                  # should then become output from vignette.
usethis::use_data(fullResults, overwrite = TRUE)

usethis::use_data(trendResults, overwrite = TRUE)

usethis::use_data(fullResults.MLE, overwrite = TRUE) # TODO seems to be contained in
                                        # fullResults, so may not need
