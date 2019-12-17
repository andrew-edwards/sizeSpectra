# Simulated data for Figures 4 and 5 for MEPS paper that take too long to produce when building a vignette.
#  See ?sizeSpectra::MLEbin.MEPS.default  etc. for details

MLEbin.MEPS.default <- MLEbin.simulate()
usethis::use_data(MLEbin.MEPS.default, overwrite = TRUE)

MLEbin.MEPS.xmin16 <- MLEbin.simulate(xmin.known = 16)
usethis::use_data(MLEbin.MEPS.xmin16, overwrite = TRUE)

MLEbin.MEPS.cutoff16 <- MLEbin.simulate(cut.off = 16)
usethis::use_data(MLEbin.MEPS.cutoff16, overwrite = TRUE)
