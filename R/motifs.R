
#' Transcription Factor Motif Info
#'
#' Transcription factor binding motif information.
#'
#' This is a special data set that stores transcription factor binding motifs
#' (for several genome builds), which can be used with multiome data to compute epiregulons.
#'
#' @param metadata logical flag specifying whether to retrieve the data (FALSE) or metadata only (TRUE)
#' @param genome character string specifying the genome build
#'
#' @return
#' A \code{GRangeList} object containing binding site information of transcription factor
#' ChIP-seq combined from Cistrome database and ENCODE.
#'
#' @section Data storage and access:
#' All motif info objects (\code{GRangeList}) are stored in one hdf5 file.
#' The function will return the one specified in the call.
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-motifs.Rmd", package = "scMultiome")}
#' ```
#'
#' @export
#'
#'
motifs <-
    function(metadata = FALSE,
             genome = c("hg38", "hg19", "mm10")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(genome, several.ok = FALSE)


        eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(), c("scMultiome", "motifs"))
        ans <- if (metadata) {
            eh[genome]
        } else {
            fileName <- eh[[genome]]
            fc <- rhdf5::h5ls(fileName)
            TFs <- fc[fc[["group"]] == sprintf("/%s", genome), "name"]
            motifInfo <- lapply(TFs, function(x) rhdf5::h5read(file = fileName, name = sprintf("/%s/%s", genome, x)))
            motifInfo <- lapply(motifInfo, restoreGR)
            names(motifInfo) <- TFs
            GenomicRanges::GRangesList(motifInfo, compress = FALSE)
        }

        return(ans)
    }
