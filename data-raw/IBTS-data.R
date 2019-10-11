# The full IBTS data set used in the MEPS paper.

load("ibtsQ1cpuelength.RData")    # Original file as downloaded by Julia
                                  # Blanchard. Contains data frame q1.
dataOrig = dplyr::tbl_df(q1)

usethis::use_data(dataOrig, overwrite = TRUE)


load("nSeaFungImport.RData") # TODO prob change me

data = dplyr::ungroup(data)   # Else some groups get kept
IBTS.data = data

usethis::use_data(IBTS.data, overwrite = TRUE)
