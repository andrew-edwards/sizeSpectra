# The MLEbins analyses of the IBTS data in the MEPS paper.

# dataBin is created (and saved once) in vignettes/MEPS_IBTS_MLEbins.Rmd

load("nSeaMLEbins.RData")         # Output from nSeaMLEbins.Snw. TODO
                                  # should then become output from vignette maybe





# TODO some results also, so someone can compare them
# Need MLEbins.res, but it's done in vignette. Well, it is saved in

# From the end of MEPS_IBTS_MLEbins.Rmd (since this was originally in
# nSeaMLEbins-recommend.Snw):

# TODO: double check .Rmd code didn't change since fcb88a3
MLEbins.res = MLEbins.nSeaFung.new
MLEbins.res = dplyr::mutate(MLEbins.res,
                            C = (b != -1 ) * (b+1) / ( xmax^(b+1) - xmin^(b+1) ) +
                                (b == -1) * 1 / ( log(xmax) - log(xmin) )
                           )
MLEbins.res = dplyr::select(MLEbins.res, -stdErr)

usethis::use_data(MLEbins.res, overwrite = TRUE)
