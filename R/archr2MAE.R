#' Create a MultiAssayExperiment from an ArchR project
#'
#' Import the matrices from a ArchR project to into a MultiAssayExperiment where each ArchR matrix is a SingleCellExperiment. Metadata is not added to the MultiAssayExperiment or Experiments.
#'
#' @param archrDir String specifying the ArchR project directory. The directory should contain a Save-ArchR-Project.rds.
#' @param numThreads Integer specifying threads to be used for parallel computing during ArchR project setup
#'
#' @return A \linkS4class{MultiAssayExperiment}
#'
#' @author Natalie Fox
#' @export
#' @importFrom ArchR getArchRThreads
#' @importFrom SummarizedExperiment colData
#' @importFrom MultiAssayExperiment colData colData<- MultiAssayExperiment ExperimentList listToMap
#' @importFrom S4Vectors DataFrame

archr2MAE <- function(archrDir, numThreads = 1) {

    # Get list of Single Cell Experiments
    all.exp <- create.exp.list.from.archr(archrDir = archrDir,
                                          numThreads = numThreads)

    # reorder experiments so tile matrices are first in MultiAssayExperiment
    tile.matrix.names <- names(all.exp)[grep('TileMatrix', names(all.exp))]
    se.names <- c(tile.matrix.names[order(as.numeric(sub('TileMatrix','', tile.matrix.names)))],
                  setdiff(names(all.exp), tile.matrix.names))

    all.exp <- all.exp[se.names]

    # create the sample Map for the MultiAssayExperiment
    el <- ExperimentList(all.exp)
    maplist <- lapply(all.exp, function(se) {
        data.frame(primary = se$Sample, colname = colnames(se), stringsAsFactors = FALSE)})

    sampMap <- listToMap(maplist)

    # create the MultiAssayExperiment
    mae <- MultiAssayExperiment(el,
                                sampleMap = sampMap,
                                colData = DataFrame(row.names=unique(sampMap$primary)));

    return(mae)
}

#' @importFrom ArchR getAvailableMatrices loadArchRProject getArchRThreads getArchRLogging addArchRLogging
#' @importFrom BiocGenerics start
#' @importFrom methods as
#'
#' @keywords internal
create.exp.list.from.archr <- function(archrDir,
                                       numThreads = getArchRThreads()) {

    archr.logging <- getArchRLogging()
    addArchRLogging(useLogs = FALSE)

    # Check that there is a saved ArchR project in the directory
    if(! 'Save-ArchR-Project.rds' %in% list.files(archrDir)) {
        stop('Are you sure this is an ArchR project directory?
             No Save-ArchR-Project.rds file found in the ArchR project directory')
    }

    archrProj <- suppressMessages(loadArchRProject(archrDir))

    # map which Embeddings should be saved to which SingleCellExperiment
    embedding.map <- create.embedding.map(archrProj)

    # import the existing matrices into Experiment objects
    all.exp <- list()
    available.matrices <- getAvailableMatrices(archrProj)

    for(matrix.type in available.matrices) {
        sce <- load.matrix(matrix.type, archrProj, embedding.map, numThreads)
        if(matrix.type == 'TileMatrix') {
            tile.size <- start(rowRanges(sce))[2] - start(rowRanges(sce))[1]
            tile.size.check <- start(rowRanges(sce))[3] - start(rowRanges(sce))[2]
            if(tile.size != tile.size.check) {
                stop('not sure of the tile size')
            }
            matrix.type <- paste0(matrix.type, tile.size)
        }

        all.exp[[matrix.type]] <- sce
    }

    addArchRLogging(useLogs = archr.logging)

    return(all.exp)
}

#' @importFrom ArchR getEmbedding
#' @keywords internal
create.embedding.map <- function(archrProj) {
    reduced.dim.names <- names(archrProj@reducedDims)
    reduced.dim.matrix.type <- unlist(lapply(reduced.dim.names,
                                             function(x) {archrProj@reducedDims[[x]]$useMatrix}))

    embedding.colnames <- unlist(lapply(names(archrProj@embeddings),
                                 function(x) {colnames(getEmbedding(archrProj, x))}))
    embedding.reduced.dim.mapping <- unique(sub('#.*','', embedding.colnames))
    embedding.map <- list()

    for(reduced.dim.name in names(reduced.dim.matrix.type)) {
        matrix.name <- reduced.dim.matrix.type[reduced.dim.name]

        if(matrix.name == 'TileMatrix') {
            matrix.name <- paste0(matrix.name,archrProj@reducedDims[[reduced.dim.name]]$tileSize);
        }
        embedding.map[[matrix.name]] <- c(
            embedding.map[[matrix.name]],
            colnames(embedding.reduced.dim.mapping)[which(embedding.reduced.dim.mapping == reduced.dim.name)]
        )
    }

    # match Embeddings without a specified matrix to TileMatrix500
    unspecified.embedding.names <- setdiff(colnames(embedding.reduced.dim.mapping), unlist(embedding.map))
    if(length(unspecified.embedding.names) > 0) {
        embedding.map[['TileMatrix500']] <- c(embedding.map[['TileMatrix500']], unspecified.embedding.names)
    }

    return(embedding.map)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start end end<-
#' @importFrom SummarizedExperiment rowData
#' @importFrom GenomeInfoDb seqlevels sortSeqlevels
#' @keywords internal
assign.tile.rowranges <- function(se, chrSizes) {
    # Create the GRanges for the tile matrix
    se.row.data <- rowData(se)
    tile.size <- se.row.data$start[2] - se.row.data$start[1];
    tile.size.check <- se.row.data$start[3] - se.row.data$start[2];
    if(tile.size != tile.size.check) {
        stop('not sure of the tile size')
    }
    se.tile.size <- tile.size
    if(length(tile.size) == 1) {
        rowRanges(se) <- GRanges(
            seqnames = rowData(se)$seqnames,
            ranges = IRanges(
                start = rowData(se)$start+1,
                end = rowData(se)$start+tile.size
            )
        )
        # adjust the tiles on the end of the chromosomes to match the chromosome ends
        for(chr in seqlevels(rowRanges(se))){
            chr.end <- end(chrSizes)[as.logical(seqnames(chrSizes) == chr)]
            end(rowRanges(se))[as.logical(seqnames(rowRanges(se)) == chr) &
                                   end(rowRanges(se)) > chr.end] <- chr.end
        }
        rowRanges(se) <- sortSeqlevels(rowRanges(se))
        se <- se[order(rowRanges(se)),]
        rownames(se) <- paste0(seqnames(rowRanges(se)),':',start(rowRanges(se)),'-',end(rowRanges(se)))
    }
    return(se)
}



#' @importFrom SingleCellExperiment SingleCellExperiment reducedDims<-
#' @importFrom ArchR getMatrixFromProject getEmbedding getArrowFiles getGenomeAnnotation getCellColData getReducedDims
#' @importFrom SummarizedExperiment rowRanges assayNames assayNames<-
#' @importFrom GenomicRanges width
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom S4Vectors SimpleList
#' @importFrom methods as
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom rhdf5 h5read
#' @keywords internal
load.matrix <- function(matrix.type, archrProj, embedding.map = NULL, numThreads = getArchRThreads()) {
    message(matrix.type)

    # check if the matrix was created as binary or not
    binarize <- (h5read(getArrowFiles(archrProj)[1], paste0(matrix.type,'/Info/Class')) == 'Sparse.Binary.Matrix')

    # Get the SummarizedExperimets with one column per cell and then convert to SingleCellExperiments
    se <- suppressMessages(getMatrixFromProject(archrProj, matrix.type, binarize = binarize, threads = numThreads))

    # check that colData colnames match the getCellColData colnames and replace if not matching
    if('projColData[, colnames(projColData) %ni% colnames(colData)]' %in% colnames(colData(se))) {
        cell.col.data <- getCellColData(archrProj)
        missing.colname <- setdiff(colnames(cell.col.data),colnames(colData(se)))
        if(length(missing.colname) == 1 &&
           length(which(colnames(colData(se)) == 'projColData[, colnames(projColData) %ni% colnames(colData)]')) == 1) {
            colnames(colData(se))[which(colnames(colData(se)) ==
                                            'projColData[, colnames(projColData) %ni% colnames(colData)]')] <- missing.colname;
        }
    }


    # Create the rowRanges for the experiment
    if(matrix.type == 'TileMatrix') {
        se <- assign.tile.rowranges(se,chrSizes = getGenomeAnnotation(archrProj)$chromSizes)
        tile.size <- width(rowRanges(se))[1]
        matrix.type <- paste0(matrix.type,tile.size)
    } else if(is.null(rowRanges(se)) & all(c('start','end','strand') %in% names(rowData(se)))) {
        row.data <- rowData(se)
        row.data$start[row.data$strand == 2] <- rowData(se)$end[row.data$strand == 2]
        row.data$end[row.data$strand == 2] <- rowData(se)$start[row.data$strand == 2]
        rowRanges(se) <- sortSeqlevels(GRanges(
            seqnames = row.data$seqnames,
            ranges = IRanges(
                start = row.data$start,
                end = row.data$end),
            strand = c('+','-','*')[as.numeric(row.data$strand)]
        ))
        # columns 1 through 4 are the start, end, seqnames and strand so remove them from the rowData
        rowData(se) <- row.data[,-seq(1,4)]
    }

    sce <- as(se, 'SingleCellExperiment')

    # change the assay Name
    if(any(grep('Matrix',assayNames(sce)))) {
        assayNames(sce) <- 'counts'
    }

    # add reduced dimensions
    reduced.dims.list <- SimpleList()
    reduced.dim.matrix.type <- unlist(sapply(
        names(archrProj@reducedDims),
        function(x) {
            id <- archrProj@reducedDims[[x]]$useMatrix;
            if(is.null(id)) {
                id <- 'TileMatrix'
                x <- which(unlist(sapply(names(archrProj@reducedDims),function(x) {archrProj@reducedDims[[x]]$useMatrix == 'TileMatrix'})));
            }
            ifelse(id == 'TileMatrix',paste0(id,archrProj@reducedDims[[x]]$tileSize),id)
        }))
    for(rd.name in names(reduced.dim.matrix.type)[which(reduced.dim.matrix.type == matrix.type)]) {
        reduced.dims.list[[rd.name]] <- getReducedDims(archrProj, rd.name)
    }
    # add reduced dimensions - Embeddings
    if(matrix.type %in% names(embedding.map)) {
        embedding.list <- SimpleList()
        for(embedding.name in embedding.map[[matrix.type]]) {
            embedding.list[[embedding.name]] <- getEmbedding(archrProj, embedding.name)
        }
        reduced.dims.list <- append(reduced.dims.list, embedding.list)
    }
    if(length(reduced.dims.list) > 0) {
        for(rd.name in names(reduced.dims.list)) {
            if(all(rownames(reduced.dims.list[[rd.name]]) %in% colnames(sce))) {
                reduced.dims.list[[rd.name]] <- reduced.dims.list[[rd.name]][colnames(sce),]
            }
        }

        reducedDims(sce) <- reduced.dims.list
    }

    # remove <Rle> from colData
    colData(sce)$Sample <- as.character(colData(sce)$Sample)

    return(sce)
}
