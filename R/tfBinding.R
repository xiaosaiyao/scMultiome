#'
#' TF Binding Info
#'
#' Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE.
#' This is a special data set that stores transcription factor binding sites for human and
#' mouse genomic builds, which can be used with multiome data to compute epiregulons.
#'
#' @inheritParams prostateENZ
#' @param genome String indicating the genomic build
#'
#' @return Granges List of TF binding sites
#'
#' @format
#' \code{GRangeList} object containing binding site information of transcription factor
#' ChIP-seq merged from ChIPAtlas database and ENCODE.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{hg38}: GRangesList object of length 1558}
#'   \item{\strong{hg19}: GRangesList object of length 1558}
#'   \item{\strong{mm10}: GRangesList object of length 768}
#' }
#'
#' @references
#' ChIP-Atlas 2021 update: a data-mining suite for exploring epigenomic landscapes by
#' fully integrating ChIP-seq, ATAC-seq and Bisulfite-seq data.
#' Zou Z, Ohta T, Miura F, Oki S.
#' \emph{Nucleic Acids Research. Oxford University Press (OUP);} 2022.
#' \href{http://dx.doi.org/10.1093/nar/gkac199}{doi:10.1093/nar/gkac199}
#'
#' ChIP‐Atlas: a data‐mining suite powered by full integration of public ChIP‐seq data.
#' Oki S, Ohta T, Shioi G, Hatanaka H, Ogasawara O, Okuda Y, Kawaji H, Nakaki R, Sese J, Meno C.
#' \emph{EMBO}; Vol. 19, EMBO reports. 2018.
#' \href{http://dx.doi.org/10.15252/embr.201846255}{doi: 10.15252/embr.201846255}
#'
#' ENCODE {https://www.encodeproject.org/}
#'
#'
#' @section Data storage and access:
#' Each genomic build is a separate \code{GRangesList} object. Each \code{GRangesList}
#' is split into individual \code{GRanges} objects, converted into data frames,
#' and then stored in a single hdf5 file. Data can be accessed with a special
#' function that extracts the requested genomic build and converts the data frame back to
#' \code{GRangesList}.
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-tfBinding.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' tfBinding("hg38")
#'
#' @export
#'
#'
tfBinding <- function(genome = c("hg38", "hg19", "mm10"),
                                       metadata = FALSE) {
    checkmate::assertFlag(metadata)
    genome <- match.arg(genome, several.ok = FALSE)


    eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(), c("scMultiome", "tfBinding"))
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

