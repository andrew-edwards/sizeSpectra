# Species codes from ICES DATRAS

specCodeNames = tibble::as_tibble(read.csv("speccodes.csv"))

usethis::use_data(specCodeNames, overwrite = TRUE)
