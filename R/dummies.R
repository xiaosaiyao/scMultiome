
#' create dummy data sets
#'
#' Create summy SCE and MAE objects
#'
#' @param experiments character string of experiment names
#'
#' @return
#' \code{dummySCE} returns a \code{SingleCellExperiment}.
#' \code{dummyMAE} returns a \code{MultiAssayExperiment}.
#'
#' @name dummies
#'
#' @examples
#' scMultiome:::dummySCE()
#' scMultiome:::dummyMAE("dummyExperiment")
#'



#' @import SingleCellExperiment
#'
#' @keywords internal
#' @rdname dummies
#'
dummySCE <- function(...) {
    ncells <- 10
    nfeats <- 20

    u <- matrix(stats::rpois(nfeats * ncells, 5), ncol = ncells)
    v <- log2(u + 1)

    u <- methods::as(u, "dgCMatrix")
    v <- methods::as(v, "dgCMatrix")

    pca <- matrix(stats::runif(ncells * 5), ncells)
    tsne <- matrix(stats::rnorm(ncells * 2), ncells)

    ans <- SingleCellExperiment(
        assays = list(counts = u, logcounts = v),
        reducedDims = S4Vectors::SimpleList(PCA = pca, tSNE = tsne))
    rownames(ans) <- paste("feat", seq_len(nfeats), sep = "_")
    colnames(ans) <- paste("cell", seq_len(ncells), sep = "_")

    colData(ans) <- S4Vectors::DataFrame(
        "Cell" = colnames(ans),
        "one" = sample(letters, ncells),
        "two" = stats::rnorm(ncells),
        row.names = colnames(ans)
    )
    rowData(ans) <- S4Vectors::DataFrame(
        "Feat" = rownames(ans),
        "idx" = seq_along(rownames(ans)),
        row.names = rownames(ans)
    )
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

    alt.mat.1 <- matrix(stats::rpois(nfeats * ncells, 5), ncol = ncells)
    colnames(alt.mat.1) <- colnames(ans)

    alt.mat.2 <- matrix(stats::rpois(nfeats * ncells, 5), ncol = ncells)
    colnames(alt.mat.2) <- colnames(ans)

    altExp(ans, "Hash") <- SingleCellExperiment(list(counts = alt.mat.1))

    altExp(ans, "Protein") <- SingleCellExperiment(list(counts = alt.mat.2))

    return(ans)
}


#' @import MultiAssayExperiment
#'
#' @keywords internal
#' @rdname dummies
#'
dummyMAE <- function(experiments = c("EXP1", "EXP2")) {
    checkmate::assertCharacter(experiments)

    expList <- lapply(experiments, dummySCE)
    names(expList) <- experiments
    mae <- MultiAssayExperiment(experiments = expList)
    return(mae)
}
