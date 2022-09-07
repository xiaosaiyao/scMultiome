
#' determine dataset being queried
#'
#' Obtain name of dataset from the (grand)parent call.
#'
#' This function reads the name of the function two steps up the call stack and
#' returns is as a string. It is called by \code{loadExp} and \code{getMetadata}
#' in order to determine the multiome dataset being queried.
#'
#' @return A character string reflecting the name of the dataset.
#'
#' @keywords internal
#'
determineDataset <- function() {
    deparse(sys.call(-2L)[[1]])
}
