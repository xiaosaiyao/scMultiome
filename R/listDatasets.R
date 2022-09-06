
#' list all available datasets
#'
#' Summary information for all datasets available in the package.
#'
#' @return
#' A \code{DataFrame} listing all available datasets,
#' with one datasest per row and the following columns:
#' \itemize{
#'   \item{Title: data set name}
#'   \item{Species: species name}
#'   \item{Type: sample type, e.g. cell culture, tissue}
#'   \item{Multiome: paired or unpaired}
#'   \item{DiskSize: size of the dataset in storage (also size of the download)}
#'   \item{MemorySize: size of the dataset in memory}
#'   \item{Accessor: function used to access the dataset directly}
#' }
#'
#' @export
#'
listDatasets <- function () {
    path <- system.file("extdata", "manifest.csv", package = "scMultiome")
    ans <- S4Vectors::DataFrame(utils::read.csv(path, stringsAsFactors = FALSE))

    return(ans)
}
