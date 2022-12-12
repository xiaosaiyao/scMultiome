
#' save and load data sets
#'
#' Functions to disassemble and save, and load and reassemble MultiAssayExperiment data sets.
#'
#' These are utilities for developers to add new data sets to the package.
#' Most will usually be called internally.
#'
#' \code{saveMAE} is used to save a \code{MultiAssayExperiment} to a hdf5 file.
#' It creates the \code{file} and passes individual \code{experiments} to \code{saveExp}.
#'
#' \code{saveExp} is called by \code{saveMAE} to disassemble experiment \code{exp}
#' and save its elements in \code{file}. A group hierarchy is created with
#' the top level group called \code{expName}, e.g. "GeneScoreMatrix",
#' and lower level groups to store specific elements:
#' \itemize{
#'   \item{\code{experiment} class is saved in group "class"}
#'   \item{assays are saved as sparse matrices in subgroup "assays", using \code{writeSparseMatrix}}
#'   \item{\code{colData} and \code{colnames} are saved in subgroup "properties"}
#'   \item{if \code{experiment} has \code{rownames}, the are saved in "properties"}
#'   \item{\code{rowData} is saved in "properties", unless it is an empty \code{DataFrame}}
#'   \item{if \code{experiment} has a \code{rowRanges} component,
#'         it is converted to a data frame and saved in "properties"}
#'   \item{if \code{experiment} has a \code{metadata} component,
#'         it is deparsed to a string and saved in "properties"}
#'   \item{if \code{experiment} has a \code{reducedDims} component,
#'         they are saved in subgroup "reducedDims"}
#'   \item{if \code{experiment} has a \code{altExps} component,
#'         they are saved in subgroup "altExps" by recursively calling \code{saveExp}}
#' }
#' \code{DataFrames} (e.g. \code{colData}, \code{rowData}, embeddings)
#' are converted to \code{data.frames} and saved as compound type.
#'
#' \code{loadMAE} is called by accessor functions to retrieve data. It locates the hdf5 file
#' in which the data set is stored and uses \code{loadExp} to extract the specified \code{experiments}.
#'
#' \code{loadExp} first checks which property elements are stored for the experiment in question,
#' loads all elements of the experiment and builds a \code{SummarizedExperiment} object.
#' If the experiment was originally a \code{SingleCellExperiment} or a subclass thereof,
#' that class as well as possible additional slots are restored.
#'
#' \code{testFile} can be used to test whether a data set loads correctly from a local file.
#' It calls \code{loadMAE} and extracts all experiment verbosely.
#'
#' \code{uploadFile} will upload a single file to Bioconductor's staging directory.
#'
#' @name fileOperations
#'
#' @param mae object of class \code{MultiAssayExperiment}
#' @param file path to a hdf5 file
#' @param experiments character string specifying which experiments to save/load
#' @param verbose logical flag specifying operation verbosity
#' @param overwrite logical flag specifying whether to allow overwriting the hdf5 file
#' @param exp an experiment object that inherits from class \code{SummarizedExperiment},
#'            usually a \code{SingleCellExperiment}
#' @param expName name of the experiment, i.e. name of the \code{ArchR} Matrix
#' @param sasToken access token to \code{endpoint}
#' @param endpoint Bioconductor's data bucket endpoint url
#'
#' @return
#' \code{saveExp} returns \code{TRUE} invisibly if the save was successful.
#' \code{saveMAE} returns a named list of \code{TRUE} values.
#' \code{loadExp} returns a \code{SingleCellExperiment} or an object of a subclass.
#' \code{loadMAE} returns a \code{MultiAssayExperiment}.
#' \code{testFile} returns the \code{MultiAssayExperiment} stored in \code{file} in its entirety.
#' \code{uploadFile} returns TRUE invisibly.
#'
#' @seealso
#' Vignette "rhdf5 - HDF5 interface for R" (\code{vignette} or \code{browseVignettes})
#' for details of hdf5 file construction.\cr
#' \code{writeSparseMatrix} for details of saving sparse matrices.
#'
#' @examples
#'
#' # create dummy MultiAssayExperiment
#' mae <- scMultiome:::dummyMAE()
#'
#' fileName <- tempfile(fileext = ".h5")
#' saveMAE(mae, fileName)                                # save MAE
#' remae <- loadMAE(fileName, c("EXP1", "EXP2"), TRUE)   # load MAE (internal)
#' remae_exp1 <- loadMAE(fileName, "EXP1", TRUE)         # load MAE with one experiment (internal)
#'
#' # create dummy SingleCellExperiment
#' sce <- scMultiome:::dummySCE()
#'
#' saveExp(sce, "EXP3", fileName, TRUE)                  # save one experiment (internal)
#' resce <- loadExp(fileName, "EXP3", TRUE)              # load one experiment (internal)
#'
#' testFile(fileName)                                    # load whole MAE



#' @import MultiAssayExperiment
#'
#' @rdname fileOperations
#'
#' @export
#' @author Aleksander Chlebowski and Xiaosai
#'
saveMAE <- function(mae, file, experiments = NULL, verbose = TRUE, overwrite = FALSE) {
    # TODO:
    #     consider saving colData only once in root rather once for each experiment
    #     would save storage but interfere with listing available experiments
    checkmate::assertClass(mae, "MultiAssayExperiment")
    lapply(experiments(mae), checkmate::assertClass, classes = "SummarizedExperiment")
    checkmate::assertCharacter(experiments, null.ok = TRUE)
    checkmate::assertChoice(experiments, names(mae), null.ok = TRUE)
    checkmate::assertDirectoryExists(dirname(file), access = "w")
    checkmate::assertFlag(verbose)
    checkmate::assertFlag(overwrite)

    if (file.exists(file)) {
        if (overwrite) {
            if (verbose) message("file ", file, " exists and will be overwritten")
            unlink(file, force = TRUE)
        } else {
            stop("file ", file, " already exists")
        }
    }

    if (is.null(experiments)) experiments <- names(mae)

    # create file
    if(verbose) message("creating h5 file ")
    rhdf5::h5createFile(file)

    # save and track progress
    ans <- lapply(experiments, function(x) {
        if(verbose) message("saving experiment ", x)
        ans <- saveExp(mae[[x]], x, file, verbose)
        if (verbose) message(" ... done", "\n")
        return(ans)
    })
    names(ans) <- experiments

    return(ans)
}



#' @import MultiAssayExperiment
#'
#' @rdname fileOperations
#'
#' @export
#'
loadMAE <- function(file, experiments, verbose) {
    checkmate::assertFileExists(file, access = "r")
    assertHDF5(file)
    checkmate::assertCharacter(experiments)
    checkmate::assertFlag(verbose)

    if (verbose) {
        dataset <- dynGet("dataset",
                          ifnotfound = paste("stored in file",
                                             dynGet("file",
                                                    ifnotfound = stop("resource missing"))))
        message("loading dataset ", dataset)
    }

    # list experiments in file
    fileContents <- rhdf5::h5ls(file)
    expStored <- fileContents[fileContents[["group"]] == "/", "name"]
    lapply(experiments, checkmate::assertChoice, choices = expStored)

    # load files
    if (verbose) message("loading experiments")
    expList <- lapply(experiments, loadExp, file = file, verbose = verbose)
    names(expList) <- experiments

    # build mae
    if (verbose) message("building MultiAssayExperiment")
    mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = expList)

    return(mae)
}



#' @rdname fileOperations
#'
#' @export
#'
methods::setGeneric(
    name = "saveExp",
    def = function(exp, expName, file, verbose) standardGeneric("saveExp"))



#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
#' @rdname fileOperations
#'
#' @export
#'
methods::setMethod(
    f = "saveExp",
    signature = "SummarizedExperiment",
    definition = function(exp, expName, file, verbose) {
        checkmate::assertString(expName)
        checkmate::assertFileExists(file, access = "w")
        assertHDF5(file)
        checkmate::assertFlag(verbose)

        # PART ONE: disassemble experiment
        if (verbose) message("\t dissassembling experiment")
        assays <- assays(exp)
        colData <- SummarizedExperiment::colData(exp)
        rowData <- SummarizedExperiment::rowData(exp)
        colnames <- colnames(exp)
        rownames <- rownames(exp)
        rowRanges <- SummarizedExperiment::rowRanges(exp)
        metadata <- S4Vectors::metadata(exp)

        # PART TWO: save experiment elements
        # create master group for this experiment
        rhdf5::h5createGroup(file, expName)

        # write object class
        if (verbose) message("\t writing class")
        rhdf5::h5write(
            obj = class(exp),
            file = file,
            name = sprintf("%s/class", expName))

        # write assays
        if (verbose) message("\t writing assays")
        rhdf5::h5createGroup(file, sprintf("%s/assays", expName))
        for (ass in names(assays)) {
            if (verbose) message("\t ... ", ass)
            writeSparseMatrix(
                x = methods::as(assays[[ass]],"dgCMatrix"),
                file = file,
                name = sprintf("%s/assays/%s", expName, ass))
        }

        # write properties
        if (verbose) message("\t writing properties")
        rhdf5::h5createGroup(file, sprintf("%s/properties", expName))

        # colData is only saved if it has columns
        if (ncol(colData) != 0L) {
            if (verbose) message("\t ...colData")
            rhdf5::h5write(
                obj = as.data.frame(colData),
                file = file,
                name = sprintf("%s/properties/colData", expName))
        }

        # write colnames
        if (verbose) message("\t ...colnames")
        rhdf5::h5write(
            obj = colnames,
            file = file,
            name = sprintf("%s/properties/colnames", expName))

        # rowData is only saved if it has columns
        # column-less rowData will be built from rownames
        if (ncol(rowData) != 0L) {
            if (verbose) message("\t ...rowData")
            rhdf5::h5write(
                obj = as.data.frame(rowData),
                file = file,
                name = sprintf("%s/properties/rowData", expName))
        }

        # rownames are only saved if there are any
        if (!is.null(rownames)) {
            if (verbose) message("\t ...rownames")
            rhdf5::h5write(
                obj = rownames,
                file = file,
                name = sprintf("%s/properties/rownames", expName))
        }

        # write rowRanges
        if (inherits(rowRanges, "GRanges")) {
            if (verbose) message("\t ... rowRanges")
            rhdf5::h5write(
                obj = storeGR(rowRanges),
                file = file,
                name = sprintf("%s/properties/rowRanges", expName)
            )
        } else if (inherits(rowRanges, "DataFrame")) {
            if (verbose) message("\t ... rowRanges")
            rhdf5::h5write(
                obj = as.data.frame(rowRanges),
                file = file,
                name = sprintf("%s/properties/rowRanges", expName)
            )
         }

        # write metadata, if any
        if (!is.null(metadata)) {
            if (verbose) message("\t ... metadata")
            rhdf5::h5write(
                obj = deparse(metadata),
                file = file,
                name = sprintf("%s/properties/metadata", expName)
            )
        }

        return(invisible(TRUE))
    })



#' @import SingleCellExperiment
#'
#' @rdname fileOperations
#'
#' @export
#'
methods::setMethod(
    f = "saveExp",
    signature = "SingleCellExperiment",
    definition = function(exp, expName, file, verbose) {
        checkmate::assertString(expName)
        checkmate::assertFileExists(file, access = "w")
        assertHDF5(file)
        checkmate::assertFlag(verbose)

        methods::callNextMethod()

        # write reduced dimensions, if any
        reducedDims <- SingleCellExperiment::reducedDims(exp)

        if (length(reducedDims) > 0L) {
            if (verbose) message("\t writing reducedDims")
            rhdf5::h5createGroup(file, sprintf("%s/reducedDims", expName))
            for (red in names(reducedDims)) {
                if (verbose) message("\t ... ", red)
                rhdf5::h5write(
                    obj = reducedDims[[red]],
                    file = file,
                    name = sprintf("%s/reducedDims/%s", expName, red))
            }
        }

        # write alternative experiments, if any
        altExps <- SingleCellExperiment::altExps(exp)

        if (length(altExps) > 0L) {
            if (verbose) message("\t writing altExps")
            rhdf5::h5createGroup(file, sprintf("%s/altExps", expName))
            for (alt in names(altExps)) {
                if (verbose) message("\t ... ", alt)
                saveExp(exp = altExps[[alt]],
                        expName = sprintf("%s/altExps/%s", expName, alt),
                        file = file,
                        verbose = verbose)
            }
        }
        return(invisible(TRUE))
    })



#' @import SingleCellExperiment
#'
#' @rdname fileOperations
#'
#' @export
#'
loadExp <- function(file, expName, verbose) {
    checkmate::assertFileExists(file, access = "r")
    assertHDF5(file)
    checkmate::assertString(expName)
    checkmate::assertFlag(verbose)

    # load class
    expClass <- rhdf5::h5read(file, sprintf("%s/class", expName))

    if (verbose) message("loading ", expName)

    # list file contents
    fileContents <- rhdf5::h5ls(file)
    colnames.present <- is.element("colnames",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    colData.present <- is.element("colData",
                                  fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    rownames.present <- is.element("rownames",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    rowData.present <- is.element("rowData",
                                  fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    rowRanges.present <- is.element("rowRanges",
                                    fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    metadata.present <- is.element("metadata",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/properties", expName), "name"])
    reducedDims.present <- is.element("reducedDims",
                                      fileContents[fileContents[["group"]] == sprintf("/%s", expName), "name"])
    altExps.present <- is.element("altExps",
                                  fileContents[fileContents[["group"]] == sprintf("/%s", expName), "name"])


    # load everything in group specified by experiment

    # load properties
    if (verbose) message("\t loading properties")

    # cell information
    if (colnames.present) {
        if (verbose) message("\t ... colnames")
        colnames <- rhdf5::h5read(file, sprintf("%s/properties/colnames", expName))
    }
    if (colData.present) {
        if (verbose) message("\t ... colData")
        colData <- rhdf5::h5read(file, sprintf("%s/properties/colData", expName))
    }

    # feature information
    if (rownames.present) {
        if (verbose) message("\t ... rownames")
        rownames <- rhdf5::h5read(file, sprintf("%s/properties/rownames", expName))
    }
    if (rowData.present) {
        if (verbose) message("\t ... rowData")
        rowData <- rhdf5::h5read(file, sprintf("%s/properties/rowData", expName))
    }

    # load row ranges if saved
    if (rowRanges.present) {
        if (verbose) message("\t ... rowRanges")
        rowRanges <- rhdf5::h5read(file = file, name = sprintf("%s/properties/rowRanges", expName))
        rowRanges <- restoreGR(rowRanges)
    }

    # load metadata if saved
    if (metadata.present) {
        if (verbose) message("\t ... metadata")
        metadata <- rhdf5::h5read(file = file, name = sprintf("%s/properties/metadata", expName))
        metadata <- eval(parse(text = metadata))
    }

    # load assays
    if (verbose) message("\t loading assays")
    assaysStored <- fileContents[fileContents[["group"]] == sprintf("/%s/assays", expName), "name"]
    assays <- lapply(assaysStored, function(ass) {
        if (verbose) message("\t ... ", ass)
        HDF5Array::H5SparseMatrix(filepath = file, group = sprintf("/%s/assays/%s", expName, ass))
    })
    names(assays) <- assaysStored

    # load reduced dimensions if saved
    if (reducedDims.present) {
        if (verbose) message("\t loading reducedDims")
        reducedDimsStored <- fileContents[fileContents[["group"]] == sprintf("/%s/reducedDims", expName), "name"]
        reducedDims <- lapply(reducedDimsStored, function(red) {
            if (verbose) message("\t ... ", red)
            rhdf5::h5read(file = file, name = sprintf("%s/reducedDims/%s", expName, red))
        })
        names(reducedDims) <- reducedDimsStored
    }

    # load alternative experiment
    if (altExps.present) {
        if (verbose) message("\t loading altExps")
        altExpsStored <- fileContents[fileContents[["group"]] == sprintf("/%s/altExps", expName), "name"]
        altExps <- lapply(altExpsStored, function(alt) {
            if (verbose) message("\t ... ", alt)
            loadExp(file = file,
                    expName = sprintf("%s/altExps/%s", expName, alt),
                    verbose = verbose)
        })
        names(altExps) <- altExpsStored
    }


    # rebuild experiment
    if (verbose) message("\t rebuilding experiment")
    ans <- SummarizedExperiment::SummarizedExperiment(assays = assays)
    if (colData.present) {
        if (verbose) message("\t ...colData")
        SummarizedExperiment::colData(ans) <- S4Vectors::DataFrame(colData)
    }
    if (colnames.present) {
        if (verbose) message("\t ...colnames")
        colnames(ans) <- colnames
    }
    if (rowData.present) {
        if (verbose) message("\t ...rowData")
        SummarizedExperiment::rowData(ans) <- S4Vectors::DataFrame(rowData)
    }
    if (rowRanges.present) {
        if (verbose) message("\t ...rowRanges")
        SummarizedExperiment::rowRanges(ans) <- rowRanges
    }
    if (rownames.present) {
        if (verbose) message("\t ...rownames")
        rownames(ans) <- rownames
    }
    if (metadata.present) {
        if (verbose) message("\t ...metadata")
        S4Vectors::metadata(ans) <- metadata
    }

    # if inherits from SCE, convert and add more slots
    if (methods::extends(expClass, "SingleCellExperiment")) {
        # direct conversion from (Ranged)SummarizedExperiment fails
        # due to different requirements for the int_elementMetadata slot
        ans <- methods::as(ans, "SingleCellExperiment")
        if (!methods::is(ans, expClass)) {
            ans <- methods::as(ans, expClass)
        }

        if (reducedDims.present) {
            if (verbose) message("\t ...reducedDims")
            SingleCellExperiment::reducedDims(ans, withDimnames = FALSE) <- reducedDims
        }
        if (altExps.present) {
            if (verbose) message("\t ...altExps")
            SingleCellExperiment::altExps(ans, withDimnames = FALSE) <- altExps
        }
    }

    return(ans)
}



#' @rdname fileOperations
#'
#' @export
#'
testFile <- function(file) {
    checkmate::assertFileExists(file, access = "r")
    assertHDF5(file)

    # list file contents
    fc <- rhdf5::h5ls(file)
    # list experiments (groups in root)
    experiments <- fc[fc$group == "/", "name"]

    message("TESTING loading MAE from file:\t", file)
    message("LOADING ALL EXPERIMENTS:\t", paste(experiments, collapse = ", "), "\n")
    ans <- loadMAE(file, experiments, verbose = TRUE)
    return(ans)
}



#' @rdname fileOperations
#'
#' @export
#'
uploadFile <- function(file, sasToken, endpoint = "https://bioconductorhubs.blob.core.windows.net") {
    checkmate::assertFileExists(file, access = "r")
    checkmate::assertString(sasToken)
    checkmate::assertString(endpoint)

    if (!requireNamespace("AzureStor")) {
        stop("uploading files requires package AzureStor, which is not insatlled")
    }

    message("establishing connection")
    # create endpoint object
    ep <- AzureStor::storage_endpoint(endpoint, sas = sasToken)
    # create container
    container <- AzureStor::storage_container(ep, "staginghub")

    message("files present in staging directory")
    # list files
    show(AzureStor::list_storage_files(container))

    message("commencing upload")
    # upload file
    AzureStor::storage_upload(
        container,
        src = normalizePath(file),
        dest = file.path("scMultiome", basename(file)))

    message("upload complete")
    # list files
    show(AzureStor::list_storage_files(container))

    return(invisible(TRUE))
}
