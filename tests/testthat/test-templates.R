
datasetName <- "fictitiousDataset"

ansMD <- makeMakeData(datasetName)
ansMM <- makeMakeMetadata(datasetName)
ansMR <- makeR(datasetName)

path.Rmd <- file.path(system.file("scripts", package = "scMultiome"), paste0("make-data-", datasetName, ".Rmd"))
path.metaR <- file.path(system.file("scripts", package = "scMultiome"), paste0("make-metadata-", datasetName, ".R"))
path.R <- file.path(system.file("R", package = "scMultiome"), paste0(datasetName, ".R"))

test_that("file ceration", {
    expect_true(ansMD)
    expect_true(ansMM)
    expect_true(ansMR)
})

test_that("files exist", {
    expect_true(file.exists(path.Rmd))
    expect_true(file.exists(path.metaR))
    expect_true(file.exists(path.R))
})


unlink(c(path.Rmd, path.metaR, path.R))
