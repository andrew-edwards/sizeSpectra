# Figures 4 and 5 for MEPS paper TODO... that take too long to produce when building a vignette. Code here ensures reproducibility.
# Change filename of data when decided on a format for multiple data sets.


Have done these two:

MLEbin.MEPS.default <- MLEbin.simulate()
usethis::use_data(MLEbin.MEPS.default, overwrite = TRUE)

MLEbin.MEPS.xmin16 <- MLEbin.simulate(xmin.known = 16)
usethis::use_data(MLEbin.MEPS.xmin16, overwrite = TRUE)

MLEbin.MEPS.cutoff16 <- MLEbin.simulate(cut.off = 16)
usethis::use_data(MLEbin.MEPS.cutoff16, overwrite = TRUE)
