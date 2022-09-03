
#' retrieve data set metadata
#'
#' Retrieve and return metadata for the currently queried dataset in ExperimentHub format.
#'
#' This is an internal function that is used by accessor functions when metadata is requested.
#' Dataset name is derived from the name of the accessor function which calls this one.
#'
#' @return An \code{ExperimentHub} object.
#'
getMetadata <- function() {
    if (!isNamespaceLoaded("ExperimentHub")) attachNamespace("ExperimentHub") # TODO - test if this can be dropped

    # determine dataset
    dataset <- determineDataset()
    ans <- AnnotationHub::query(ExperimentHub::ExperimentHub(), c("scMultiome", dataset))
    return(ans)
}

