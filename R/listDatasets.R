
#' list all available data sets
#'
#' Summary information for all data sets available in the package.
#'
#' @return
#' A \code{DataFrame} listing all available data sets,
#' with one data set per row and the following columns:
#' \itemize{
#'   \item{Call: function call used to access the data set directly}
#'   \item{Author: original data set author}
#'   \item{Title: data set name}
#'   \item{Species: species name}
#'   \item{Lineage: sample lineage}
#'   \item{CellNumber: number of cells in the data set}
#'   \item{Multiome: paired or unpaired}
#'   \item{DiskSize: size of the dataset in storage (also size of the download)}
#'   \item{Version: data set version number or upload date}
#' }
#'
#' @examples
#' listDatasets()
#'
#' @export
#'
listDatasets <- function () {
    path <- system.file("extdata", "manifest.csv", package = "scMultiome")
    ans <- S4Vectors::DataFrame(utils::read.csv(path, stringsAsFactors = FALSE))

    return(ans)
}
