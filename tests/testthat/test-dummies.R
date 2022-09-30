

test_that("dummy SCE is created", {
  expect_s4_class(dummySCE(), "SingleCellExperiment")
})

test_that("dummy MAE is created", {
  expect_s4_class(dummyMAE(), "MultiAssayExperiment")
})
