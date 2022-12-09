#'
#' TF Binding Info
#'
#' Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE or
#' from CistromeDB and ENCODE.
#' This is a special data set that stores transcription factor binding sites for human and
#' mouse genomic builds, which can be used with the package epiregulon to compute regulons.
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
#'   \item{\strong{hg38_cistrome}: GRangesList object of length 1269}
#'   \item{\strong{hg19_cistrome}: GRangesList object of length 1271}
#'   \item{\strong{mm10_cistrome}: GRangesList object of length 544}
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
#' Cistrome Data Browser: expanded datasets and new tools for gene regulatory analysis.
#' Zheng R, Wan C, Mei S, Qin Q, Wu Q, Sun H, Chen CH, Brown M, Zhang X, Meyer CA, Liu XS
#' \emph{Nucleic Acids Res}, 2018 Nov 20.
#' \href{https://academic.oup.com/nar/article/47/D1/D729/5193328}{doi:10.1093/nar/gky1094}
#'
#' Cistrome data browser: a data portal for ChIP-Seq and chromatin accessibility data in human and mouse.
#' Mei S, Qin Q, Wu Q, Sun H, Zheng R, Zang C, Zhu M, Wu J, Shi X, Taing L, Liu T, Brown M, Meyer CA, Liu XS
#' \emph{Nucleic Acids Res}, 2017 Jan 4;45(D1):D658-D662.
#' \href{https://academic.oup.com/nar/article/45/D1/D658/2333932}{doi:10.1093/nar/gkw983}
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
#' tfBinding("mm10")
#'
#' @export
#'
tfBinding <- function(genome = c("hg38", "hg19", "mm10", "hg38_cistrome", "hg19_cistrome", "mm10_cistrome"),
                      metadata = FALSE) {
    checkmate::assertFlag(metadata)
    genome <- match.arg(genome, several.ok = FALSE)

    genomeName <- sprintf("tfBinding_%s", genome)
    eh <- AnnotationHub::query(ExperimentHub::ExperimentHub(), c("scMultiome", genomeName))

    ans <-
        if (metadata) {
            eh[genomeName]
        } else {
            readRDS(eh[[genomeName]])
        }

    return(ans)
}
