#'
#' TF Binding Info
#'
#' Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE
#'
#' This is a special data set that stores transcription factor binding sites for human and
#' mouse genomic builds, which can be used with the package epiregulon to compute regulons.
#'
#' @inheritParams prostateENZ
#' @param genome character string specifying the genomic build
#' @param source character string specifying the ChIP-seq data source
#'
#' @return A list of TF binding sites as a \code{GrangesList} object.
#'
#' @format
#' \code{GRangesList} object containing binding site information
#' of transcription factor ChIP-seq.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{hg38_atlas}: GRangesList object of length 1558}
#'   \item{\strong{hg19_atlas}: GRangesList object of length 1558}
#'   \item{\strong{mm10_atlas}: GRangesList object of length 768}
#'   \item{\strong{hg38_encode.sample}: List object of length 171}
#'   \item{\strong{hg19_encode.sample}: List object of length 171}
#'   \item{\strong{mm10_encode.sample}: List object of length 31}
#'   \item{\strong{hg38_atlas.sample}: List object of length 1112}
#'   \item{\strong{hg19_atlas.sample}: List object of length 1112}
#'   \item{\strong{mm10_atlas.sample}: List object of length 517}
#'   \item{\strong{hg38_atlas.tissue}: List object of length 22}
#'   \item{\strong{hg19_atlas.tissue}: List object of length 22}
#'   \item{\strong{mm10_atlas.tissue}: List object of length 23}
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
#' \href{http://dx.doi.org/10.15252/embr.201846255}{doi:10.15252/embr.201846255}
#'
#' ENCODE: {https://www.encodeproject.org/}
#'
#' @section Data storage and access:
#' Each genomic build is a separate \code{GRangesList} object, stored in a separate RDS file.
#' All genomic builds can be accessed with the same function \code{tfBinding}.
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-tfBinding.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' # check metada of dataset
#' tfBinding("mm10", metadata = TRUE)
#' # download data
#' \dontrun{
#' tfBinding("mm10", "atlas")
#' }
#'
#' @export
#'
tfBinding <- function(genome = c("hg38", "hg19", "mm10"),
                      source = c("atlas", "encode.sample", "atlas.sample","atlas.tissue"),
                      metadata = FALSE) {
    checkmate::assertFlag(metadata)
    genome <- match.arg(genome, several.ok = FALSE)
    source <- match.arg(source, several.ok = FALSE)

    eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(),
                               pattern = c("scMultiome", "tfBinding", source, genome))

    if (source %in% c("atlas")) {
        eh_ID <- sort(eh$ah_id)[1]
    } else {
        eh_ID <- eh$ah_id
    }


    ans <-
        if (metadata) {
            eh[eh_ID]
        } else {
            readRDS(eh[[eh_ID]])
        }

    return(ans)
}
