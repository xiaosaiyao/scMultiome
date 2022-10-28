
#' check if a file is a valid HDF5 file
#'
#' Check if a file path argument points to an non-corrupt HDF5 file.
#'
#' Compares the first 8 bytes of a file to those of the standard HDF5 file header.
#'
#' @param path path to file to test
#'
#' @return
#' Returns invisible \code{path} if check is successful, otherwise signals an error.
#'
#' @examples
#' fileName1 <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(fileName1)
#' rhdf5::h5write(mtcars, fileName1, "mtcars")
#' rhdf5::h5closeAll()
#' fileName2 <-  tempfile(fileext = ".csv")
#' write.csv(mtcars, fileName2)
#'
#' assertHDF5(fileName1)   # passes
#' \donttest{
#' assertHDF5(fileName2)   # fails
#' }
#'
#' @section Further development:
#' The HDF5 file header contains 8 bytes, which hold specific meanings.
#' Currently the function only tests that the header of the file
#' specified by \code{path} is identical to a healthy HDF5 file and signals
#' a general error if that is not the case.
#' Reporting specific types of corruption can be implemented.
#'
#' @references
#' <http://web.ics.purdue.edu/~aai/HDF5/html/H5.format.html#BootBlock>
#'
#' @export
#' @author Aleksander Chlebowski
assertHDF5 <- function(path) {

  fileHead <- readBin(path, "raw", n = 8L)
  hdf5Head <- as.raw(c(137L, 72L, 68L, 70L, 13L, 10L, 26L, 10L))

  ident <- identical(fileHead, hdf5Head)

  if (ident) {
    return(invisible(path))
  } else {
    stop("Assertion on \"path\" failed: Must be a hdf5 file.", call. = FALSE)
  }
}
