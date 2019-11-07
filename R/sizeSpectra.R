## As instructed from using usethis::use_tibble() [else saved data tibbles don't
##  get used properly]

## usethis namespace: start
#' @importFrom tibble tibble
## usethis namespace: end
NULL

# As for gfiphc repo, need these to avoid warnings related to dplyr
# commands. These aren't actually needed but some might.
# if (getRversion() >= "2.15.1") utils::globalVariables(c("."))
# if (getRversion() >= "2.15.1") {
#   utils::globalVariables(c(
#           "xsmallticks",
#           "xbigticks",
#           "vertCol"
#  ))
#}
