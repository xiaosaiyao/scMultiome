
#' create dummy data sets
#'
#' Create dummy SCE and MAE objects.
#'
#' @param features character string specifying which (optional) features to create in the SCE
#' @param experiments named list of character vectors
#'                    specifying experiments to create and their features
#'
#' @return
#' \code{dummySCE} returns a \code{SingleCellExperiment}.
#' \code{dummyMAE} returns a \code{MultiAssayExperiment}.
#'
#' @name dummies
#'
#' @examples
#' scMultiome:::dummySCE()
#' scMultiome:::dummySCE("rowData")
#' scMultiome:::dummySCE("rowRanges")
#' scMultiome:::dummySCE("reducedDims")
#' scMultiome:::dummySCE("altExps")
#'
#' scMultiome:::dummyMAE(list("dummyExperiment" = NULL))
#'



#' @import SingleCellExperiment
#'
#' @keywords internal
#' @rdname dummies
#'
dummySCE <- function(features = c("rowData", "rowRanges", "reducedDims", "altExps", "none")) {
    features <- match.arg(features, several.ok = TRUE)

    ncells <- 10
    nfeats <- 20

    u <- matrix(stats::rpois(nfeats * ncells, 5), ncol = ncells)
    v <- log2(u + 1)

    u <- methods::as(u, "dgCMatrix")
    v <- methods::as(v, "dgCMatrix")

    pca <- matrix(stats::runif(ncells * 5), ncells)
    tsne <- matrix(stats::rnorm(ncells * 2), ncells)

    ans <- SingleCellExperiment(assays = list(counts = u, logcounts = v))

    rownames(ans) <- paste("feat", seq_len(nfeats), sep = "_")
    colnames(ans) <- paste("cell", seq_len(ncells), sep = "_")
    colData(ans) <- S4Vectors::DataFrame(
        "Cell" = colnames(ans),
        "one" = sample(letters, ncells),
        "two" = stats::rnorm(ncells),
        row.names = colnames(ans)
    )
    if (is.element("rowRanges", features)) {
        rowRanges(ans) <- GenomicRanges::GRanges(
            data.frame(
                "seqnames" = paste("chr", seq_along(rownames(ans))),
                "start" = seq_along(rownames(ans)),
                "end" = seq_along(rownames(ans)) * 5,
                "width" = seq_along(rownames(ans)) * 5,
                "strand" = rep_len(c("+", "-"), length(rownames(ans))),
                row.names = rownames(ans)
            )
        )
    }
    if (is.element("rowData", features)) {
        rowData(ans) <- S4Vectors::DataFrame(
            "Feat" = paste("chr", seq_along(rownames(ans))),
            "idx" = seq_along(rownames(ans)),
            row.names = rownames(ans)
        )
    }
    if (is.element("reducedDims", features)) {
        reducedDims(ans) <- S4Vectors::SimpleList(PCA = pca, tSNE = tsne)
    }
    if (is.element("altExps", features)) {
        alt.mat.1 <- matrix(stats::rpois(nfeats * ncells, 5), ncol = ncells)
        colnames(alt.mat.1) <- colnames(ans)

        alt.mat.2 <- matrix(stats::rpois(nfeats * ncells, 5), ncol = ncells)
        colnames(alt.mat.2) <- colnames(ans)

        altExp(ans, "Hash") <- SingleCellExperiment(list(counts = alt.mat.1))
        altExp(ans, "Protein") <- SingleCellExperiment(list(counts = alt.mat.2))
    }

    return(ans)
}



#' @import MultiAssayExperiment
#'
#' @keywords internal
#' @rdname dummies
#'
dummyMAE <- function(experiments = list("EXP1" = NULL, "EXP2" = NULL)) {
    checkmate::assertList(experiments)
    lapply(experiments, checkmate::assertCharacter, null.ok = TRUE)

    experiments <- lapply(experiments, function(x) {
        if (is.null(x)) c("rowData", "rowRanges", "reducedDims", "altExps") else x
    })

    expList <- lapply(experiments, dummySCE)
    mae <- MultiAssayExperiment(experiments = expList)
    return(mae)
}
