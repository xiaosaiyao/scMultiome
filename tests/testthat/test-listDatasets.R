
fileName <- system.file("extdata", "manifest.csv", package = "scMultiome")

test_that("locating manifest file", {
    expect_true(file.exists(fileName))
})

test_that("DataFrame is returned", {
    ans <- listDatasets()
    expect_s4_class(ans, "DataFrame")
})
