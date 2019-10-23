# Figures 4 and 5 for MEPS paper TODO... that take too long to produce when building a vignette. Code here ensures reproducibility.
# Change filename of data when decided on a format for multiple data sets.



MLEbin.MEPS.default <- MLEbin.simulate()
usethis::use_data(MLEbin.MEPS.default, overwrite = TRUE)


# Do default value then do:
expect_equal(MLEbin.MEPS.default$MLE.array, MLE.array)
expect_equal(MLEbin.MEPS.default$MLEconf.array, MLEconf.array)
expect_equal(MLEbin.MEPS.default$MLE.array.parameters, MLE.array.known)  # may
                                        # not work with name change?

# Can delete these then:
usethis::use_data(MLE.array, overwrite = TRUE)
usethis::use_data(MLEconf.array, overwrite = TRUE)
usethis::use_data(MLE.array.known, overwrite = TRUE)


# Then expect to do these:

MLEbin.MEPS.xmin16 <- MLEbin.simulate(xmin.known = 16)
usethis::use_data(MLEbin.MEPS.xmin16, overwrite = TRUE)

MLEbin.MEPS.cutoff16 <- MLEbin.simulate(xmin.known = 16, ***)
usethis::use_data(MLEbin.MEPS.cutoff16, overwrite = TRUE)
