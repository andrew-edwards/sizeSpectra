# The full IBTS data set used in the MEPS paper.

# Needs running once to save tibbles as tibbles (it puts tibble in Imports)
usethis::use_tibble()

load("ibtsQ1cpuelength.RData")    # Original file as downloaded by Julia
                                  # Blanchard. Contains data frame q1.
dataOrig = tibble::as_tibble(q1)  # TODO - not run yet with as_tibble

usethis::use_data(dataOrig, overwrite = TRUE)

# See TODO vignette
load("nSeaFungImport.RData") # TODO prob change me

IBTS_data = tibble::as_tibble(data)

IBTS_data = dplyr::ungroup(IBTS_data)   # Else some groups get kept

usethis::use_data(IBTS_data, overwrite = TRUE)

## use_tibble() said to document a returned tibble like so: [maybe add when
#doing documentation, TODO]
#' @return a [tibble][tibble::tibble-package]
#'
