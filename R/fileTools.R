
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import MultiAssayExperiment

saveMAE <- function(mae, file, experiments = NULL, verbose = TRUE, overwrite = FALSE) {
    # TODO: save colData once in root - but that will interfere with listing available experiments!
    checkmate::assertClass(mae, "MultiAssayExperiment")
    checkmate::assertCharacter(experiments, null.ok = TRUE)
    checkmate::assertChoice(experiments, names(mae), null.ok = TRUE)
    checkmate::assertDirectoryExists(dirname(file), access = "w")
    if (tools::file_ext(file) != "h5") stop("\"file\" must be a .h5 file")
    checkmate::assertFlag(verbose)
    checkmate::assertFlag(overwrite)

    if (file.exists(file)) {
        if (overwrite) {
            if (verbose) message("file ", file, " exists and will be overwritten")
            unlink(file, force = TRUE)
        } else {
            stop("saveMAE: file ", file, " already exists")
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
        if (verbose) message("\t ... done")
        return(ans)
    })
    names(ans) <- experiments

    return(ans)
}



loadMAE <- function(experiments, verbose = FALSE) {
    checkmate::assertCharacter(experiments)
    checkmate::assertFlag(verbose)

    # determine dataset
    dataset <- determineDataset()

    if (verbose) message("loading dataset ", dataset)

    # check for file
    if (verbose) message("\t locating file")
    file <- file.path(system.file("extdata", package = "scMultiome"), sprintf("%s.h5", dataset))
    checkmate::assertFileExists(file, access = "r", extension = "h5")

    # list experiments in file
    if (verbose) message("\t checking experiments")
    fileContents <- rhdf5::h5ls(file)
    expStored <- fileContents[fileContents[["group"]] == "/", "name"]
    lapply(experiments, checkmate::assertChoice, choices = expStored)

    # load files
    if (verbose) message("\t loading experiments")
    expList <- sapply(experiments, loadExp, file = file, verbose = verbose, simplify = FALSE)

    # build mae
    if (verbose) message("\t building MultiAssayExperiment")
    mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = expList)

    return(mae)
}



saveExp <- function(exp, expName, file, verbose) {
    checkmate::assertClass(exp, "SummarizedExperiment")
    checkmate::assertString(expName)
    checkmate::assertFileExists(file, access = "w", extension = "h5")

    # PART ONE: disassemble experiment
    if (verbose) message("\t dissassembling experiment")
    assays <- assays(exp)
    colData <- SummarizedExperiment::colData(exp)
    rowData <- SummarizedExperiment::rowData(exp)
    colnames <- colnames(exp)
    rownames <- rownames(exp)
    reducedDims <- SingleCellExperiment::reducedDims(exp)

    # PART TWO: save experiment elements

    # create groups for this experiment within file
    if (verbose) message("\t creating groups")
    rhdf5::h5createGroup(file, expName)
    rhdf5::h5createGroup(file, sprintf("%s/assays", expName))
    rhdf5::h5createGroup(file, sprintf("%s/metadata", expName))
    # reducedDims group is only created if there any reducedDims objects
    if (length(reducedDims) > 0L) {
        rhdf5::h5createGroup(file, sprintf("%s/reducedDims", expName))
    }

    # write object class
    if (verbose) message("\t writing class")
    rhdf5::h5write(obj = class(exp), file = file, name = sprintf("%s/class", expName))

    # write package that defines object class
    if (verbose) message("\t writing parent package")
    rhdf5::h5write(obj = attr(class(exp), "package"), file = file, name = sprintf("%s/package", expName))

    # write assays
    if (verbose) message("\t writing assays")
    for (ass in names(assays)) {
        if (verbose) message("\t\t ", ass)
        writeSparseMatrix(assays[[ass]], file = file, name = sprintf("%s/assays/%s", expName, ass))
    }

    # write metadata
    if (verbose) message("\t writing metadata")
    rhdf5::h5write(obj = as.data.frame(colData), file = file, name = sprintf("%s/metadata/colData", expName))
    # rowData is only saved if it has columns
    # column-less rowData will be built from rownames
    if (ncol(rowData) != 0L) {
        rhdf5::h5write(obj = as.data.frame(rowData), file = file, name = sprintf("%s/metadata/rowData", expName))
    }
    rhdf5::h5write(obj = colnames, file = file, name = sprintf("%s/metadata/colnames", expName))
    # rownames are only saved if there are any
    if (!is.null(rownames)) {
        rhdf5::h5write(obj = rownames, file = file, name = sprintf("%s/metadata/rownames", expName))
    }
    # write reduced dimensions, if any
    if (length(reducedDims) > 0L) {
        if (verbose) message("\t writing reduced dimensions")
        for (red in names(reducedDims)) {
            if (verbose) message("\t\t ", red)
            rhdf5::h5write(obj = reducedDims[[red]], file = file,
                           name = sprintf("%s/reducedDims/%s", expName, red))
        }
    }

    return(invisible(TRUE))
}



loadExp <- function(file, expName, verbose) {
    checkmate::assertFileExists(file, access = "r", extension = "h5")
    checkmate::assertString(expName)
    checkmate::assertFlag(verbose)

    if (verbose) message("loading ", expName)

    # list file contents
    fileContents <- rhdf5::h5ls(file)
    rownames.present <- is.element("rownames",
                                   fileContents[fileContents[["group"]] == sprintf("/%s/metadata", expName), "name"])
    rowData.present <- is.element("rowData",
                                  fileContents[fileContents[["group"]] == sprintf("/%s/metadata", expName), "name"])
    reducedDims.present <- is.element("reducedDims",
                                      fileContents[fileContents[["group"]] == sprintf("/%s", expName), "name"])

    # load everything in group specified by experiment

    # load class
    expClass <- rhdf5::h5read(file, sprintf("%s/class", expName))

    # load parent package of class
    expClassPkg <- rhdf5::h5read(file, sprintf("%s/package", expName))

    # load assays
    if (verbose) message("\t loading assays")
    assaysStored <- fileContents[fileContents[["group"]] == sprintf("/%s/assays", expName), "name"]
    assays <- lapply(assaysStored, function(x) {
        if (verbose) message("\t\t ", x)
        HDF5Array::H5SparseMatrix(filepath = file, group = sprintf("/%s/assays/%s", expName, x))
    })
    names(assays) <- assaysStored

    # load metadata
    if (verbose) message("\t loading metadata")
    colnames <- rhdf5::h5read(file, sprintf("%s/metadata/colnames", expName))
    colData <- rhdf5::h5read(file, sprintf("%s/metadata/colData", expName))
    # load rownames if saved
    if (rownames.present) {
        rownames <- rhdf5::h5read(file, sprintf("%s/metadata/rownames", expName))
    }
    # load rowData if saved saved, otherwise construct column-less data frame from rownames
    rowData <- if (rowData.present) {
        rhdf5::h5read(file, sprintf("%s/metadata/rowData", expName))
    } else {
        data.frame(row.names = rownames)
    }

    # load reduced dimensions if saved
    if (reducedDims.present) {
        if (verbose) message("\t loading reduced dimensions")
        reducedDimsStored <- fileContents[fileContents[["group"]] == sprintf("/%s/reducedDims", expName), "name"]
        reducedDims <- lapply(reducedDimsStored, function(x) {
            if (verbose) message("\t\t ", x)
            rhdf5::h5read(file = file, name = sprintf("%s/reducedDims/%s", expName, x))
        })
        names(reducedDims) <- reducedDimsStored
    }

    # rebuild experiment
    if (verbose) message("\t rebuilding experiment")
    ans <- SingleCellExperiment::SingleCellExperiment(assays = assays,
                                                      rowData = rowData,
                                                      colData = colData)
    colnames(ans) <- colnames
    if (rownames.present) {
        rownames(ans) <- rownames
    }
    if (reducedDims.present) {
        SingleCellExperiment::reducedDims(ans, withDimnames = FALSE) <- reducedDims
    }

    # attempt to restore class of experiment if different to SingleCellExperiment
    # only if required package is available
    if (expClass != "SingleCellExperiment") {
        if (requireNamespace(expClassPkg, quietly = TRUE)) {
            if (verbose) message("\t converting to ", expClass)
            ans <- switch(expClass,
                          # class-specific conversions

                          "SingleCellAccessibilityExperiment" = {
                              tileSize <- as.integer(sub("([A-Za-z]+)(\\d+)", "\\2", "TileMatrix500"))
                              dsdb.plus::SingleCellAccessibilityExperiment(ans, tile.size = tileSize)
                          },

                          # generic conversion (default)
                          methods::as(ans, expClass))
        } else {
            if (verbose) {
                message("\t cannot convert to class ", expClass, " as package ", expClassPkg, " is unavailable")
                message("\t ... returning as SingleCellExperiment")
            }
        }
    }

    return(ans)
}



testFile <- function(file) {
    checkmate::assertFileExists(file, access = "r", extension = "h5")

    # list file contents
    fc <- rhdf5::h5ls(file)
    # list experiments (groups in root)
    experiments <- fc[fc$group == "/", "name"]

    # define accessor function
    fooDef <- function() {
        cat("testing accessor function:\t", deparse(match.call()[[1]]), "\n")
        cat("loading all experiments:\t", paste(experiments, collapse = ", "), "\n")
        ans <- loadMAE(experiments, verbose = TRUE)
        return(ans)
    }
    # obtain function name
    fooName <- tools::file_path_sans_ext(basename(file))

    # assign accessor with proper name
    assign(fooName, fooDef)

    # call accessor
    eval(call(fooName))
}



uploadFile <- function(file, sasToken, endpoint = "https://bioconductorhubs.blob.core.windows.net") {
    # checkmate::assertFileExists(file, access = "r", extension = "h5")
    checkmate::assertString(sasToken)
    checkmate::assertString(endpoint, pattern = "https://\\w+\\.\\w+\\.\\w+\\.\\w+\\.[a-z]{3}")

    message("establishing connection")
    # create endpoint object
    ep <- AzureStor::storage_endpoint(endpoint, sas = sasToken)
    # create container
    container <- AzureStor::storage_container(ep, "staginghub")

    message("files present in staging directory")
    # list files
    print(AzureStor::list_storage_files(container))

    message("commencing upload")
    # upload file
    AzureStor::storage_upload(container, src = normalizePath(file), dest = file.path("scMultiome", basename(file)))

    message("upload complete")
    # list files
    print(AzureStor::list_storage_files(container))

    return(invisible(TRUE))
}
