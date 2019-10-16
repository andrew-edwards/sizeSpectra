# Species codes from ICES DATRAS

specCodeNames = read.csv("speccodes.csv")

usethis::use_data(specCodeNames, overwrite = TRUE)
