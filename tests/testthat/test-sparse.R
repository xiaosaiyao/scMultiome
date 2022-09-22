# This tests the sparse readers/writers.
# library(testthat); library(artificer.matrix); source("test-sparse.R")

library(Matrix)

test_that("writing to a sparse matrix works as expected", {
    for (i in 1:3) {
        x <- rsparsematrix(100, 20, 0.5)
        if (i == 2) {
            x <- DelayedArray(x) * 1 # force use of the block method.
        } else if (i == 3) {
            x <- BiocGenerics::cbind(DelayedArray(x), DelayedArray(rsparsematrix(100, 10, 0.5)))
        }

        tmp <- tempfile(fileext=".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        writeSparseMatrix(x, tmp, "csr_matrix", column = FALSE)
        writeSparseMatrix(x, tmp, "tenx_matrix", tenx = TRUE)

        library(HDF5Array)
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(x)))
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csr_matrix"))), unname(as.matrix(x)))
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "tenx_matrix"))), unname(as.matrix(x)))
    }
})

test_that("writing to a sparse matrix works with tiny chunks", {
    # Forcing use of the block method.
    x <- DelayedArray(rsparsematrix(100, 20, 0.5)) * 2

    tmp <- tempfile(fileext=".h5")
    setAutoBlockSize(max(dim(x))*8)
    writeSparseMatrix(x, tmp, "csc_matrix")
    writeSparseMatrix(x, tmp, "csr_matrix", column = FALSE)

    library(HDF5Array)
    expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(x)))
    expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csr_matrix"))), unname(as.matrix(x)))

    setAutoBlockSize()
})

get_type <- function(tmp, path) {
    fhandle <- H5Fopen(tmp)
    on.exit(H5Fclose(fhandle))
    dhandle <- H5Dopen(fhandle, path)
    on.exit(H5Dclose(dhandle), add = TRUE, after = FALSE)
    .Call("_getDatatypeName", H5Dget_type(dhandle), PACKAGE = "rhdf5")
}

test_that("writing to a sparse matrix works with guessed type", {
    set.seed(1000)
    for (i in 1:3) {
        core <- function() round(rsparsematrix(100, 20, 0.5))
        if (i == 1) {
            FUN <- function(f) f(core())
        } else if (i == 2) {
            FUN <- function(f) DelayedArray(f(core())) * 1 # force use of the block method.
        } else if (i == 3) {
            FUN <- function(f) BiocGenerics::cbind(DelayedArray(f(core())), DelayedArray(f(round(rsparsematrix(100, 10, 0.5)))))
        }

        # Signed:
        tmp <- tempfile(fileext = ".h5")
        y <- FUN(identity)
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "I32")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))

        # Unsigned short
        tmp <- tempfile(fileext = ".h5")
        y <- FUN(function(x) abs(x * 1000))
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "U16")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))

        # Unsigned word
        tmp <- tempfile(fileext = ".h5")
        y <- FUN(function(x) abs(x * 1000000))
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "I32")
        expect_equal(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))

        # Unsigned long
        tmp <- tempfile(fileext = ".h5")
        y <- FUN(function(x) abs(x * 10000000000))
        writeSparseMatrix(y, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/data"), "F64")
        expect_identical(unname(as.matrix(H5SparseMatrix(tmp, "csc_matrix"))), unname(as.matrix(y)))
    }
})

test_that("writing to a sparse matrix works with guessed index type for block method", {
    set.seed(1000)
    for (i in 1:2) {
        x <- round(rsparsematrix(100, 20, 0.5))
        if (i != 1) {
            x <- DelayedArray(x) * 1 # force use of the block method.
        }

        # 16-bit:
        tmp <- tempfile(fileext = ".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/indices"), "U16")

        # 32-bit:
        x <- round(rsparsematrix(100000, 20, 0.001))
        tmp <- tempfile(fileext = ".h5")
        writeSparseMatrix(x, tmp, "csc_matrix")
        expect_match(get_type(tmp, "csc_matrix/indices"), "U32")
    }
})

test_that("autoloader recognizes the sparse matrix", {
    x <- rsparsematrix(100, 20, 0.2)
    tmp <- tempfile(fileext = ".h5")

    library(rhdf5)
    h5createFile(tmp)
    h5createGroup(tmp, "samp_data")
    writeSparseMatrix(x, tmp, "samp_data/data")

    out <- .createRawMatrixSeed(list(`$schema` = "sparse_matrix/v2.json",
                                     sparse_matrix = list(dimensions = dim(x))), tmp)
    expect_s4_class(out, "H5SparseMatrixSeed")
})
