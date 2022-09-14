
#' retrieve data set or its metadata
#'
#' Retrieve and return the data or metadata for the currently queried dataset.
#'
#' This is a generic accessor function that is used by user-level accessor functions
#' to access datasets or their metadata.
#' Dataset name is derived from the name of the accessor function that called this one.
#'
#' @param metadata logical flag specifying whether to return the resource or only its metadata
#' @param experiments character string specifying which experiments to extract
#' @param verbose logical flag specifying loading verbosity
#'
#' @return
#' If \code{metadata = FALSE}, a \code{MultiAssayExperiment},
#' otherwise an \code{ExperimentHub} object.
#'
#' @import AnnotationHub
#' @import ExperimentHub
#'
#' @keywords internal
#'
retrieve <- function(metadata, experiments, verbose = FALSE) {
    checkmate::assertFlag(metadata)
    checkmate::assertCharacter(experiments)
    checkmate::assertFlag(verbose)

    # determine data set
    dataset <- determineDataset()
    # query database for the data set
    eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(), c("scMultiome", dataset))
    # backstop for missing data set
    if (length(eh) == 0L) stop("resource ", dataset, " not found in ExperimentHub")
    # retrieve meatadata or data (i.e. path to file in bucket)
    ans <- if (metadata) {
        eh[1]
    } else {
        loadMAE(file = eh[[1]], experiments, verbose)
    }

    return(ans)
}
