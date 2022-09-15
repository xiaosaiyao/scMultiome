
#' file utilities
#'
#' Helper functions for new data sets.
#'
#' \code{testFile} can be used to test whether a data set loads correctly from a local file.
#' It calls \code{loadMAE} and extracts all experiment verbosely.
#'
#' \code{uploadFile} will load a single file to the Bioconductor staging directory.
#'
#' Both functions are internal as they are meant to be run in developer mode.
#'
#' @name fileUtils
#'
#' @param file path to hdf5 file containing a data set
#' @param sasToken access token to \code{endpoint}
#' @param endpoint Bioconductor's data bucket endpoint url
#'
#' @return
#' \code{testFile} returns the \code{MultiAssayExperiment} stored in \code{file} in its entirety.
#' \code{uploadFile} returns TRUE invisibly.
#'
#' @rdname fileUtils
#'
#' @keywords internal
#'
testFile <- function(file) {
    checkmate::assertFileExists(file, access = "r", extension = "h5")

    # list file contents
    fc <- rhdf5::h5ls(file)
    # list experiments (groups in root)
    experiments <- fc[fc$group == "/", "name"]

    cat("testing loading MAE from file:\t", file, "\n")
    cat("loading all experiments:\t", paste(experiments, collapse = ", "), "\n")
    ans <- loadMAE(file, experiments, verbose = TRUE)
    return(ans)
}



#' @rdname fileUtils
#'
#' @keywords internal
#'
uploadFile <- function(file, sasToken, endpoint = "https://bioconductorhubs.blob.core.windows.net") {
    # checkmate::assertFileExists(file, access = "r", extension = "h5")
    checkmate::assertString(sasToken)
    checkmate::assertString(endpoint)

    if (!requireNamespace("AzureStor")) {
        stop("uploading files requires package AzureStor, which is not insatlled")
    }

    message("establishing connection")
    # create endpoint object
    ep <- AzureStor::storage_endpoint(endpoint, sas = sasToken)
    # create container
    container <- AzureStor::storage_container(ep, "staginghub")

    message("files present in staging directory")
    # list files
    print(AzureStor::list_storage_files(container))

    message("commencing upload")
    # upload file
    AzureStor::storage_upload(
        container,
        src = normalizePath(file),
        dest = file.path("scMultiome", basename(file)))

    message("upload complete")
    # list files
    print(AzureStor::list_storage_files(container))

    return(invisible(TRUE))
}
