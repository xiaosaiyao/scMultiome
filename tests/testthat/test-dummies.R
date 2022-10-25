

# helper function for comparing SCE slots
f <- function(sce) {
    ans <- logical(4L)
    names(ans) <- c("rowData", "rowRanges", "reducedDims", "altExps")
    ans["rowData"] <- prod(dim(SummarizedExperiment::rowData(sce))) > 0L
    ans["rowRanges"] <- inherits(SummarizedExperiment::rowRanges(sce), "GRanges")
    ans["reducedDims"] <- length(SingleCellExperiment::reducedDims(sce)) > 0L
    ans["altExps"] <- length(SingleCellExperiment::altExps(sce)) > 0L
    return(names(ans[ans]))
}


test_that("dummy SCE is created with all features", {
    sce1 <- dummySCE()
    sce2 <- dummySCE(features = c("rowData", "rowRanges", "reducedDims", "altExps"))
    expect_s4_class(sce1, "SingleCellExperiment")
    expect_s4_class(sce2, "SingleCellExperiment")
    expect_identical(f(sce1), f(sce2))
})

test_that("dummy SCE can be created with some features", {
    sce1 <- dummySCE(features = c("rowData"))
    sce2 <- dummySCE(features = c("rowData", "rowRanges"))
    sce3 <- dummySCE(features = c("rowData", "rowRanges", "reducedDims"))
    sce4 <- dummySCE(features = c("rowRanges", "altExps"))
    expect_identical(f(sce1), c("rowData"))
    expect_identical(f(sce2), c("rowData", "rowRanges"))
    expect_identical(f(sce3), c("rowData", "rowRanges", "reducedDims"))
    expect_identical(f(sce4), c("rowRanges", "altExps"))
})




test_that("dummy MAE is created", {
  mae <- dummyMAE()
  expect_s4_class(mae, "MultiAssayExperiment")
  expect_equal(length(mae), 2L)
})

test_that("dummy MAE is created containing SCEs with features", {
  mae <- dummyMAE()
  expect_s4_class(mae, "MultiAssayExperiment")
  expect_length(mae, 2L)

  mae <- dummyMAE(experiments = list("E1" = c("rowRanges", "reducedDims"),
                                     "E2" = c("reducedDims", "altExps"),
                                     "E3" = c("rowData", "altExps", "rowRanges")))
  features <- lapply(MultiAssayExperiment::experiments(mae), f)
  expect_identical(vapply(features, length, integer(1L), USE.NAMES = FALSE), c(2L, 2L, 3L))
  expect_identical(names(mae), c("E1", "E2", "E3"))
})




test_that("bad arguments", {
    expect_error(
        dummySCE("all")
    )
    expect_error(
        dummyMAE(experiments = c("E1", "E2", "E3"))
    )
})


