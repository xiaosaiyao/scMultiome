#'
#' Single-cell analysis of samples from healthy human colon
#'
#' ATACseq and RNAseq data obtained by the colon tissues analysis. Samples were
#' collected from adult human donors.
#'
#' @inheritParams prostateENZ
#'
#' @inherit prostateENZ return
#'
#' @format
#' \code{MultiAssayExperiment} obtained from an \code{ArchR} project. Annotated with the Hg38 genome build.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{TileMatrix500}: SingleCellExperiment with 6062095 rows and 59231 columns}
#'   \item{\strong{GeneIntegrationMatrix}: SingleCellExperiment with 19020 rows and 59231 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 24919 rows and 59231 columns}
#'   \item{\strong{MotifMatrix}: SingleCellExperiment with 870 rows and 59231 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 406946 rows and 59231 columns}
#' }
#'
#' @references
#' 1. Zhang K, Hocker JD, Miller M, Hou X, Chiou J, Poirion OB, Qiu Y, Li YE,
#' Gaulton KJ, Wang A, Preissl S, Ren B. A single-cell atlas of
#' chromatin accessibility in the human genome.
#' Cell. 2021 Nov 24;184(24):5985-6001.e19. doi: 10.1016/j.cell.2021.10.024.
#' Epub 2021 Nov 12. PMID: 34774128; PMCID: PMC8664161.
#'
#' 2. Becker, W.R., Nevins, S.A., Chen, D.C. et al. Single-cell analyses
#' define a continuum of cell state and composition changes in the malignant
#' transformation of polyps to colorectal cancer. Nat Genet 54, 985–995 (2022).
#'  https://doi.org/10.1038/s41588-022-01088-x
#'
#' @inheritSection prostateENZ Data storage and access
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-colonHealthy.Rmd", package = "scMultiome")}
#' ```
#'
#' @export
#'
#'
colonHealthy <-
    function(metadata = FALSE,
             experiments = c("TileMatrix500",
                             "GeneIntegrationMatrix",
                             "GeneScoreMatrix",
                             "MotifMatrix",
                             "PeakMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

        retrieve("colonHealthy", metadata, experiments, verbose = FALSE)
    }
