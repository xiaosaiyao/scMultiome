
# create hdf5 file
fileName1 <- tempfile(fileext = ".h5")
rhdf5::h5createFile(fileName1)
rhdf5::h5write(mtcars, fileName1, "mtcars")
rhdf5::h5closeAll()

# create a different file
fileName2 <-  tempfile(fileext = ".csv")
write.csv(mtcars, fileName2)


test_that("passing assertion", {
    expect_invisible(assertHDF5(fileName1))
    expect_type(assertHDF5(fileName1), "character")
})
test_that("failing assertion", {
    expect_error(assertHDF5(fileName2))
})
