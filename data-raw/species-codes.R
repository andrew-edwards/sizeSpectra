# Species codes from ICES DATRAS

specCodeNames = dplyr::tbl_df(read.csv("speccodes.csv"))

usethis::use_data(specCodeNames, overwrite = TRUE)
