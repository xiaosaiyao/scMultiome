% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeSparseMatrix.R
\name{writeSparseMatrix}
\alias{writeSparseMatrix}
\title{Write a sparse matrix}
\usage{
writeSparseMatrix(
  x,
  file,
  name,
  chunk = 10000,
  column = TRUE,
  tenx = FALSE,
  guess.integer = TRUE
)
}
\arguments{
\item{x}{A sparse matrix of some sort.
This includes sparse \linkS4class{DelayedMatrix} objects.}

\item{file}{String containing a path to the HDF5 file.
The file is created if it is not already present.}

\item{name}{String containing the name of the group to store \code{x}.}

\item{chunk}{Integer scalar specifying the chunk size for the indices and values.}

\item{column}{Logical scalar indicating whether to store as compressed sparse column format.}

\item{tenx}{Logical scalar indicating whether to use the 10X compressed sparse column format.}

\item{guess.integer}{Logical scalar specifying whether to guess an appropriate integer type from \code{x}.}
}
\value{
A \code{NULL} invisibly.
The contents of \code{x} are written to \code{name} in \code{file}.
}
\description{
Writes a sparse matrix to file in a compressed sparse format.
}
\details{
This writes a sparse matrix to file in various formats:
\itemize{
\item \code{column = TRUE} and \code{tenx = FALSE} uses H5AD's \code{csr_matrix} format.
\item \code{column = FALSE} and \code{tenx = FALSE} uses H5AD's \code{csc_matrix} format.
\item \code{tenx = TRUE} uses 10X Genomics' HDF5 matrix format.
}
For the first two formats, the apparent transposition is deliberate,
because columns in R are interpreted as rows in H5AD.
This allows us to retain consistency the interpretation of samples
(columns in R, rows in H5AD) and features (vice versa).
Constructors for classes like \linkS4class{H5SparseMatrix}
will automatically transpose so no extra work is required.

If \code{guess.integer = TRUE}, we attempt to save \code{x}'s values
into the smallest type that will accommodate all of its values.
If \code{x} only contains unsigned integers, we will attempt
to save either 8-, 16- or 32-bit unsigned integers.
If \code{x} contains signed integers, we will fall back to 32-bit signed integers.
For all other values, we will fall back to double-precision floating point values.

We attempt to save \code{x}'s indices to unsigned 16-bit integers
if the relevant dimension of \code{x} is small enough.
Otherwise we will save it as an unsigned 32-bit integer.
}
\section{Note}{

This function has been kindly contributed by Aaron Lun.
}

\examples{
library(Matrix)
x <- rsparsematrix(100, 20, 0.5)
tmp <- tempfile(fileext = ".h5")
scMultiome:::writeSparseMatrix(x, tmp, "csc_matrix")
scMultiome:::writeSparseMatrix(x, tmp, "csr_matrix", column = FALSE)
scMultiome:::writeSparseMatrix(x, tmp, "tenx_matrix", tenx = TRUE)

rhdf5::h5ls(tmp)
library(HDF5Array)
H5SparseMatrix(tmp, "csc_matrix")
H5SparseMatrix(tmp, "csr_matrix")
H5SparseMatrix(tmp, "tenx_matrix")

}
\author{
Aaron Lun
}
