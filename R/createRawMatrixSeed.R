#' Matrix loading utilities
#'
#' Utilities for loading a matrix saved by \code{stageObject}.
#'
#' @param info A named list of metadata for this matrix.
#' @param path String containing the path to the file containing said matrix.
#' @param names Logical scalar indicating whether the seed should be annotated with dimnames (if available).
#'
#' @return \code{.createRawMatrixSeed} returns a seed that can be used in the \code{DelayedArray} constructor.
#'
#' \code{.extractMatrixDimnames} returns a list of two character vectors or \code{NULL}, containing the dimnames.
#'
#' @details
#' For \code{.extractDimnames}, \code{path} is expected to be a HDF5 file
#' with row and column names in \code{samp_data/features} and \code{samp_data/samples}, respectively.
#'
#' @author Aaron Lun
#'
#' @section Note:
#' This function has been kindly contributed by Aaron Lun.
#'
#' @name createRawMatrixSeed
#'
#' @keywords internal
#'
.createRawMatrixSeed <- function(info, path, names = TRUE) {
    if ("dense_matrix" %in% names(info)) {
        out <- HDF5Array::HDF5ArraySeed(filepath = path, name = "samp_data/data")
        .matrix_namer(out, path, names)
    } else if ("sparse_matrix" %in% names(info)) {
        out <- HDF5Array::H5SparseMatrixSeed(filepath = path, group = "samp_data/data")
        .matrix_namer(out, path, names)
    } else if ("delayed_array" %in% names(info)) {
        chihaya::loadDelayed(filepath = path, "data")
    } else {
        # Just assume it's one of the old objects for now.
        HDF5Array::HDF5ArraySeed(filepath = path, name = "samp_data/data")
        # stop("unsupported type '", info$`_extra`$type, "'")
    }
}

#' @keywords internal
#'
.matrix_namer <- function(seed, path, names) {
    if (names) {
        dnames <- .extractMatrixDimnames(path)
        if (!is.null(dnames)) {
            mat <- DelayedArray::DelayedArray(seed)
            dimnames(mat) <- dnames
            seed <- mat@seed
        }
    }
    seed
}

#' @rdname createRawMatrixSeed
#'
#' @keywords internal
#'
.extractMatrixDimnames <- function(path) {
    info <- rhdf5::h5ls(path)
    available <- file.path(info$group, info$name)

    dimnames <- vector("list", 2)

    ftag <- "/samp_data/__features__"
    if (ftag %in% available) {
        dimnames[[1]] <- as.character(rhdf5::h5read(path, ftag))
    }

    stag <- "/samp_data/__samples__"
    if (stag %in% available) {
        dimnames[[2]] <- as.character(rhdf5::h5read(path, stag))
    }

    dimnames
}
