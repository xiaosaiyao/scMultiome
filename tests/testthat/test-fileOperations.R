
library(S4Vectors)
library(GenomicRanges)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

f <- function() {
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
        reducedDims = SimpleList(PCA = pca, tSNE = tsne))
    rownames(ans) <- paste("feat", seq_len(nfeats), sep = "_")
    colnames(ans) <- paste("cell", seq_len(ncells), sep = "_")

    colData(ans) <- DataFrame(
        "Cell" = colnames(ans),
        "one" = sample(letters, ncells),
        "two" = rnorm(ncells),
        row.names = colnames(ans)
    )
    return(ans)
}

sce1 <- f()
rowData(sce1) <- DataFrame(
    "Feat" = rownames(sce1),
    "idx" = seq_along(rownames(sce1)),
    row.names = rownames(sce1)
)

sce2 <- f()
rowRanges(sce2) <- GRanges(
    data.frame(
        "seqnames" = paste("chr", seq_along(rownames(sce2))),
        "start" = seq_along(rownames(sce2)),
        "end" = seq_along(rownames(sce2)) * 5,
        "width" = seq_along(rownames(sce2)) * 5,
        "strand" = rep_len(c("+", "-"), length(rownames(sce2))),
        row.names = rownames(sce2)
    )
)


mae <- MultiAssayExperiment(experiments = list("EXP1" = sce1, "EXP2" = sce2))


fileName <- tempfile(fileext = ".h5")
fileName2 <- tempfile(fileext = ".csv")
write.csv(head(iris), fileName2)

# saving
suppressMessages({
    ansSave <- saveMAE(mae, fileName)
    ansSave2 <- saveMAE(mae, fileName, overwrite = TRUE)
})
test_that("saving mae", {
    expect_true(is.list(ansSave))
    expect_true(all(unlist(ansSave)))
    expect_true(is.list(ansSave2))
    expect_true(all(unlist(ansSave2)))
    expect_true(file.exists(fileName))
})


# loading
suppressMessages(
    ansLoad <- loadMAE(fileName, experiments = c("EXP1", "EXP2"), verbose = TRUE)
)
test_that("loading mae", {
    expect_s4_class(ansLoad, "MultiAssayExperiment")
    expect_identical(length(experiments(ansLoad)), 2L)
})
test_that("assay idenitty pre- and post", {
    expect_identical(assay(mae[[1]], 1), methods::as(assay(ansLoad[[1]], 1), "dgCMatrix"))
    expect_identical(assay(mae[[1]], 2), methods::as(assay(ansLoad[[1]], 2), "dgCMatrix"))
    expect_identical(assay(mae[[2]], 1), methods::as(assay(ansLoad[[2]], 1), "dgCMatrix"))
    expect_identical(assay(mae[[2]], 2), methods::as(assay(ansLoad[[2]], 2), "dgCMatrix"))
})


# loading with wrapper
suppressMessages(
    ansTest <- testFile(fileName)
)

test_that("loading mae with wrapper", {
    expect_identical(ansLoad, ansTest)
})



# failures
test_that("bad arguments", {
    expect_error(loadMAE(fileName2))
})

test_that("missing arguments", {
    expect_error(saveMAE(mae))
    expect_error(saveMAE(file = fileName))
    expect_error(loadMAE(fileName))
})

test_that("conflicts", {
    expect_error(saveMAE(mae, fileName, overwrite = FALSE))
})

unlink(fileName, force = TRUE)
unlink(fileName2, force = TRUE)

test_that("missing arguments, continued", {
    expect_error(loadMAE(fileName))
})
