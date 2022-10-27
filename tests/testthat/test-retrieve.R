
test_that("resource is returned", {
    skip("skipping until data is uploaded") # TODO - remove when ready
    testthat::expect_s4_class(retrieve("tfBinding", metadata = TRUE, experiments = "mm10"),
                              "ExperimentHub")
    testthat::expect_s4_class(retrieve("tfBinding", metadata = FALSE, experiments = "mm10"),
                              "MultiAssayExperiment")
})

test_that("missing resources signal errors", {
    expect_error(retrieve("nonexistentdataset",
                          metadata = FALSE,
                          experiments = c("experiment1", "experiment2")))
    expect_error(retrieve("nonexistentdataset",
                          metadata = TRUE,
                          experiments = c("experiment1", "experiment2")))
})
