# The full IBTS data set used in the MEPS paper.

# Needs running once to save tibbles as tibbles (it puts tibble in Imports)
usethis::use_tibble()

load("ibtsQ1cpuelength.RData")    # Original file as downloaded by Julia
                                  # Blanchard. Contains data frame q1.
dataOrig = tibble::as_tibble(q1)  # TODO - not run yet with as_tibble

usethis::use_data(dataOrig, overwrite = TRUE)

# The steps to convert the data into IBTS_data are now all in the vignette MEPS_IBTS_1.

## use_tibble() said to document a returned tibble like so: [maybe add when
#doing documentation, TODO]
#' @return a [tibble][tibble::tibble-package]
#'
