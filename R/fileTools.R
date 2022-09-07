
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
    file <- system.file("data", sprintf("%s.h5", dataset), package = "scMultiome")
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
    mae <- MultiAssayExperiment(experiments = expList)

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

    # write assays
    if (verbose) message("\t writing assays")
    for (ass in names(assays)) {
        if (verbose) message("\t\t ", ass)
        artificer.matrix::writeSparseMatrix(assays[[ass]], file = file,
                                            name = sprintf("%s/assays/%s", expName, ass))
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
    expClass <- rhdf5::h5read(file, sprintf("%s/class", expName))

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
    if (expClass == "SingleCellAccessibilityExperiment") {
        tileSize <- as.integer(sub("([A-Za-z]+)(\\d+)", "\\2", "TileMatrix500"))
        ans <- dsdb.plus::SingleCellAccessibilityExperiment(ans, tile.size = tileSize)
    }

    return(ans)
}
