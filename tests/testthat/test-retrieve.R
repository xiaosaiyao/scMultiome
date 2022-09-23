
test_that("missing resources", {
    expect_error(retrieve("nonexistentdataset", metadata = FALSE, experiments = c("experiment1", "experiment2")))
    expect_error(retrieve("nonexistentdataset", metadata = TRUE, experiments = c("experiment1", "experiment2")))
})
