
fileName <- file.path(system.file("scripts", package = "scMultiome"), "datasetList.Rmd")

test_that("file is created", {
    metadata <- file.path(system.file("extdata", "metadata.csv", package = "scMultiome"))
    metadata <- read.csv(metadata)
    expect_true(makeDataSetList(metadata))
    expect_true(file.exists(fileName))
})
