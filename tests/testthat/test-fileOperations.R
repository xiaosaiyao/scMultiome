
# sce0 <- dummySCE("none")
# sce1 <- dummySCE("rowData")
# sce2 <- dummySCE("rowRanges")
# sce3 <- dummySCE("reducedDims")
# sce4 <- dummySCE("altExps")
# sce5 <- dummySCE()

mae <- MultiAssayExperiment(experiments = list(
    "EXP0" = dummySCE("none"),
    "EXP1" = dummySCE("rowData"),
    "EXP2" = dummySCE("rowRanges"),
    "EXP3" = dummySCE("reducedDims"),
    "EXP4" = dummySCE("altExps"),
    "EXP5" = dummySCE()
))


fileName <- tempfile(fileext = ".h5")
fileName2 <- tempfile(fileext = ".csv")
write.csv(head(iris), fileName2)

# saving
suppressMessages({
    ansSave <- saveMAE(mae, fileName, verbose = TRUE)
    ansSave2 <- saveMAE(mae, fileName, verbose = FALSE, overwrite = TRUE)
})
test_that("saving mae", {
    expect_true(is.list(ansSave))
    expect_true(all(vapply(ansSave, isTRUE, logical(1L))))
    expect_true(is.list(ansSave2))
    expect_true(all(vapply(ansSave2, isTRUE, logical(1L))))
    expect_true(file.exists(fileName))
})


# loading
ansLoad <- loadMAE(fileName, experiments = sprintf("EXP%i", 0:5), verbose = FALSE)

test_that("loading mae", {
    expect_s4_class(ansLoad, "MultiAssayExperiment")
    expect_identical(length(experiments(ansLoad)), 6L)
})
test_that("assay identity pre- and post", {
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
