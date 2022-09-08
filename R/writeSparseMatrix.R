
#' Write a sparse matrix
#'
#' Writes a sparse matrix to file in a compressed sparse format.
#'
#' @param x A sparse matrix of some sort.
#' This includes sparse \linkS4class{DelayedMatrix} objects.
#' @param file String containing a path to the HDF5 file.
#' The file is created if it is not already present.
#' @param name String containing the name of the group to store \code{x}.
#' @param chunk Integer scalar specifying the chunk size for the indices and values.
#' @param column Logical scalar indicating whether to store as compressed sparse column format.
#' @param tenx Logical scalar indicating whether to use the 10X compressed sparse column format.
#' @param guess.integer Logical scalar specifying whether to guess an appropriate integer type from \code{x}.
#'
#' @details
#' This writes a sparse matrix to file in various formats:
#' \itemize{
#' \item \code{column = TRUE} and \code{tenx = FALSE} uses H5AD's \code{csr_matrix} format.
#' \item \code{column = FALSE} and \code{tenx = FALSE} uses H5AD's \code{csc_matrix} format.
#' \item \code{tenx = TRUE} uses 10X Genomics' HDF5 matrix format.
#' }
#' For the first two formats, the apparent transposition is deliberate,
#' because columns in R are interpreted as rows in H5AD.
#' This allows us to retain consistency the interpretation of samples
#' (columns in R, rows in H5AD) and features (vice versa).
#' Constructors for classes like \linkS4class{H5SparseMatrix}
#' will automatically transpose so no extra work is required.
#'
#' If \code{guess.integer = TRUE}, we attempt to save \code{x}'s values
#' into the smallest type that will accommodate all of its values.
#' If \code{x} only contains unsigned integers, we will attempt
#' to save either 8-, 16- or 32-bit unsigned integers.
#' If \code{x} contains signed integers, we will fall back to 32-bit signed integers.
#' For all other values, we will fall back to double-precision floating point values.
#'
#' We attempt to save \code{x}'s indices to unsigned 16-bit integers
#' if the relevant dimension of \code{x} is small enough.
#' Otherwise we will save it as an unsigned 32-bit integer.
#'
#' @return
#' A \code{NULL} invisibly.
#' The contents of \code{x} are written to \code{name} in \code{file}.
#'
#' @author Aaron Lun
#'
#' @section Note:
#' This function has been kindly contributed by Aaron Lun.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' x <- rsparsematrix(100, 20, 0.5)
#' tmp <- tempfile(fileext = ".h5")
#' writeSparseMatrix(x, tmp, "csc_matrix")
#' writeSparseMatrix(x, tmp, "csr_matrix", column = FALSE)
#' writeSparseMatrix(x, tmp, "tenx_matrix", tenx = TRUE)
#'
#' rhdf5::h5ls(tmp)
#' library(HDF5Array)
#' H5SparseMatrix(tmp, "csc_matrix")
#' H5SparseMatrix(tmp, "csr_matrix")
#' H5SparseMatrix(tmp, "tenx_matrix")
#' }
#'
#' @export
writeSparseMatrix <- function(x,
                              file,
                              name,
                              chunk = 10000,
                              column = TRUE,
                              tenx = FALSE,
                              guess.integer = TRUE) {
    if (!file.exists(file)) {
        rhdf5::h5createFile(file)
    }
    .write_CS_matrix(file,
                     name,
                     x,
                     chunk_dim = chunk,
                     by_column = column,
                     use_tenx = tenx,
                     guess_type = guess.integer)
    invisible(NULL)
}

.write_CS_matrix <- function(file,
                             name,
                             mat,
                             chunk_dim = 10000,
                             by_column = TRUE,
                             use_tenx = FALSE,
                             guess_type = TRUE) {
    handle <- rhdf5::H5Fopen(file)
    on.exit(rhdf5::H5Fclose(handle))
    rhdf5::h5createGroup(handle, name)

    if (use_tenx) {
        by_column <- TRUE # must be column major.
        rhdf5::h5writeDataset(dim(mat), handle, file.path(name, "shape"))
    } else {
        if (by_column) {
            format <- "csr_matrix"
        } else {
            format <- "csc_matrix"
        }

        ghandle <- rhdf5::H5Gopen(handle, name)
        on.exit(rhdf5::H5Gclose(ghandle), add = TRUE, after = FALSE)

        rhdf5::h5writeAttribute(
            format,
            ghandle,
            "encoding-type",
            variableLengthString = TRUE,
            encoding = "UTF8",
            asScalar = TRUE)
        rhdf5::h5writeAttribute(
            "0.1.0",
            ghandle,
            "encoding-version",
            variableLengthString = TRUE,
            encoding = "UTF8",
            asScalar = TRUE)
        rhdf5::h5writeAttribute(
            rev(dim(mat)),
            ghandle,
            "shape")
    }

    # Scan through and guess the best types.
    type <- .guess_integer_type(mat)
    htype <- paste0("H5T_NATIVE_", type)

    max.index <- if (by_column) {
        nrow(mat)
    } else {
        ncol(mat)
    }

    itype <- if (max.index < 2^16) {
        "H5T_NATIVE_UINT16"
    } else {
        "H5T_NATIVE_UINT32"
    }

    out <- NULL
    if (methods::is(mat, "dgCMatrix") && by_column) {
        # Directly save when possible.
        rhdf5::h5createDataset(
            handle,
            file.path(name, "data"),
            dims = length(mat@x),
            H5type = htype,
            chunk = min(length(mat@x), chunk_dim)
        )
        rhdf5::h5writeDataset(
            mat@x,
            handle,
            file.path(name, "data")
        )
        rhdf5::h5createDataset(
            handle,
            file.path(name, "indices"),
            dims = length(mat@i),
            H5type = itype,
            chunk = min(length(mat@i), chunk_dim)
        )
        rhdf5::h5writeDataset(
            mat@i,
            handle,
            file.path(name, "indices")
        )

        out <- mat@p
    } else {
        # Otherwise writing by block to save memory.
        rhdf5::h5createDataset(
            handle,
            file.path(name, "data"),
            dims = 0,
            maxdims = rhdf5::H5Sunlimited(),
            H5type = htype,
            chunk = chunk_dim
        )
        rhdf5::h5createDataset(
            handle,
            file.path(name, "indices"),
            dims = 0,
            maxdims = rhdf5::H5Sunlimited(),
            H5type = itype,
            chunk = chunk_dim
        )

        last <- 0 # make this a double to avoid overflow.

        seeds <- .is_delayed_cbind_dgCMatrix(mat)
        if (!is.null(seeds) && by_column) {
            out <- vector("list", length(seeds))
            for (i in seq_along(seeds)) {
                cout <- .cbind_dgCMatrix_sparse_writer(
                    seeds[[i]],
                    last,
                    file = handle,
                    name = name)
                out[[i]] <- cout$number
                last <- cout$last
            }
        } else {
            if (by_column) {
                grid <- DelayedArray::colAutoGrid(mat)
            } else {
                grid <- DelayedArray::rowAutoGrid(mat)
            }
            out <- vector("list", length(grid))
            for (i in seq_along(grid)) {
                block <- DelayedArray::read_sparse_block(mat, grid[[i]])
                cout <- .blockwise_sparse_writer(
                    block,
                    last,
                    file = handle,
                    name = name,
                    by_column = by_column)
                out[[i]] <- cout$number
                last <- cout$last
            }
        }

        out <- as.double(unlist(out))
        out <- c(0, cumsum(out))
    }

    iname <- file.path(name, "indptr")
    rhdf5::h5createDataset(
        handle,
        iname,
        dims   = length(out),
        H5type = "H5T_NATIVE_UINT64"
    )

    rhdf5::h5writeDataset(out, handle, iname)
}

.guess_integer_type_dsparseMatrix <- function(mat) {
    any.nonint <- FALSE
    i <- 0
    while (i <= length(mat@x)) {
        end <- i + 10000000L
        current <- mat@x[i:min(end, length(mat@x))]
        if (any(current %% 1 != 0)) {
            any.nonint <- TRUE
            break
        }
        i <- end
    }

    limits <- range(mat@x)
    list(
        negative = limits[1] < 0,
        extreme = max(abs(limits)),
        non.integer = any.nonint
    )
}

.is_delayed_cbind_dgCMatrix <- function(mat) {
    if (methods::is(mat, "DelayedMatrix")) {
        seed <- mat@seed
        if (methods::is(seed, "DelayedSetDimnames")) {
            seed <- seed@seed
        }
        if (methods::is(seed, "DelayedAbind") && seed@along == 2L) {
            if (all(vapply(seed@seeds, methods::is, class = "dgCMatrix", TRUE))) {
                return(seed@seeds)
            }
        }
    }
    return(NULL)
}

.check_worst_type <- function(sparse) {
    vals <- DelayedArray::nzdata(sparse)
    any.negative <- any(vals < 0)
    any.nonint <- any(vals != round(vals))
    extreme.val <- max(abs(vals))
    list(negative = any.negative, non.integer = any.nonint, extreme = extreme.val)
}

.guess_integer_type <- function(mat) {
    if (methods::is(mat, "dsparseMatrix")) {
        details <- .guess_integer_type_dsparseMatrix(mat)
        any.neg <- details$negative
        any.nonint <- details$non.integer
        extreme <- details$extreme
    } else {
        seeds <- .is_delayed_cbind_dgCMatrix(mat)
        if (!is.null(seeds)) {
            out <- lapply(seeds, .guess_integer_type_dsparseMatrix)
        } else {
            out <- DelayedArray::blockApply(
                mat,
                FUN = .check_worst_type,
                as.sparse = TRUE,
                grid = DelayedArray::colAutoGrid(mat))
        }
        any.neg <- any(vapply(out, function(y) y$negative, TRUE))
        any.nonint <- any(vapply(out, function(y) y$non.integer, TRUE))
        extreme <- max(unlist(lapply(out, function(y) y$extreme)))
    }

    if (any.nonint) {
        type <- "DOUBLE"
    } else if (any.neg) {
        type <- "INT32"
    } else {
        # We don't try to save it as a byte, because that gets interpreted by
        # rhdf5 as a raw vector... not helpful.
        if (extreme < 2^16) {
            type <- "UINT16"
        } else if (extreme < 2^31) {
            # Don't attempt to use unsigned 32-bit ints; values between 2^31
            # and 2^32 will fail to be converted to R's signed integers.
            type <- "INT32"
        } else {
            type <- "DOUBLE"
        }
    }

    type
}

.blockwise_sparse_writer <- function(block, last, file, name, by_column = FALSE) {
    nzdex <- DelayedArray::nzindex(block)
    if (by_column) {
        primary <- nzdex[, 2]
        secondary <- nzdex[, 1]
        ndim <- ncol(block)
    } else {
        primary <- nzdex[, 1]
        secondary <- nzdex[, 2]
        ndim <- nrow(block)
    }
    v <- DelayedArray::nzdata(block)

    o <- order(primary, secondary)
    primary <- primary[o]
    secondary <- secondary[o]
    v <- v[o]

    newlast <- last + length(primary)
    index <- list(last + seq_along(primary))

    iname <- file.path(name, "indices")
    rhdf5::h5set_extent(file, iname, newlast)
    rhdf5::h5writeDataset(
        secondary - 1L,
        file,
        iname,
        index = index)

    vname <- file.path(name, "data")
    rhdf5::h5set_extent(
        file,
        vname,
        newlast)
    rhdf5::h5writeDataset(
        v,
        file,
        vname,
        index = index)

    list(number = tabulate(primary, ndim), last = newlast)
}

.cbind_dgCMatrix_sparse_writer <- function(mat, last, file, name) {
    newlast <- last + length(mat@x)
    index <- list(last + seq_along(mat@i))

    iname <- file.path(name, "indices")
    rhdf5::h5set_extent(file, iname, newlast)
    rhdf5::h5writeDataset(
        mat@i,
        file,
        iname,
        index = index)

    vname <- file.path(name, "data")
    rhdf5::h5set_extent(file, vname, newlast)
    rhdf5::h5writeDataset(
        mat@x,
        file,
        vname,
        index = index)

    list(number = diff(mat@p), last = newlast)
}
