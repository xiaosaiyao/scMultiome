#'
#' Transcription factor motifs
#'
#' Transcription factor motifs sets from the `chromVARmotifs` R package
#'
#' This data set stores transcription factor motifs for human and
#' mouse genome, which can be used with the package epiregulon to compute scores of
#' transcription factor-regulatory element links.
#'
#' @inheritParams prostateENZ
#' @param species character string specifying the species of interest
#'
#' @return A list of position weight matrices, one for each transcription factor.
#'
#' @format
#' \code{PWMatrixList} object containing information on
#' transcription factor motifs.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{human_pwms_v2}: PWMatrixList object of length 1558}
#'   \item{\strong{mouse_pwms_v2}: PWMatrixList object of length 1558}
#' }
#'
#' @references
#' Schep, A., Wu, B., Buenrostro, J. &  Greenleaf, W. J. (2017) chromVAR:
#' inferring transcription-factor-associated accessibility from single-cell
#' epigenomic data.
#' \emph{Nature Methods} 14, 975–978.
#' \href{http://dx.doi.org/0.1038/nmeth.4401}{doi:0.1038/nmeth.4401}
#'
#' @section Data storage and access:
#' The transcription factor motifs are stored separately for each species in .rds files encoding \code{PWMatrixList}.
#' Data for both species are be accessed with the same function \code{tfMotifs}.
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-tfMotifs.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' # check metada of dataset
#' tfMotifs("mouse", metadata = TRUE)
#' # download data
#' \dontrun{
#' tfMotifs("human")
#' }
#'
#' @export
#'
tfMotifs <- function(species = c("human", "mouse"),
                      metadata = FALSE) {
    checkmate::assertFlag(metadata)
    species <- match.arg(species, several.ok = FALSE)
    species <- c(human="Homo sapiens", mouse="Mus musculus")[species]
    eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(),
                               pattern = c("scMultiome", "TF motifs", species))
    eh_ID <- eh$ah_id
    ans <-
        if (metadata) {
            eh[eh_ID]
        } else {
            readRDS(eh[[eh_ID]])
        }

    return(ans)
}
