#'
#' create make-metadata data file
#'
#' Creates a template for data set metadata file.
#'
#' @param dataset name of data set as character string
#'
#' @return Invisible TRUE.
#' Creates \code{inst/scripts/make-metadata-<dataset>.R} file and navigates to it, if possible.
#'
makeMakeMetadata <- function(dataset) {
    checkmate::assertString(dataset)

    fileName <- file.path(system.file("scripts", package = "scMultiome"),
                          paste0("make-metadata-", dataset, ".R"))
    file.create(fileName)

    writeLines(con = fileName, text = c(
        "",
        "# ExperimentHub metadata",
        "# see ?ExperimentHubData::makeExperimentHubMetadata for details",
        paste0("metadata.", dataset, " <- data.frame("),
        "    Title = \"<YOUR DATA TITLE>\",",
        "    Description = \"<YOUR DATA TITLE>\",",
        "    BiocVersion = \"<CURRENT BIOCONDUCTOR VERSION>\",",
        "    Genome = \"<GENOME VERSION USED IN ANALYSIS>\",",
        "    SourceType = \"<SOURCE FILE TYPE>\", # refers to raw data",
        "    SourceUrl = \"<SOURCE FILE LOCATION>\", # refers to raw data; optional",
        "    SourceVersion = \"<VERSION NUMEBR OR DATE>\", # no commas!",
        "    Species = \"<SPECIES NAME>\", # e.g. \"Homo sapiens\"",
        "    TaxonomyId = \"<SPECIES ID>\", # e.g. \"9606\"",
        "    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE",
        "    DataProvider = , # refers to raw data",
        "    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer",
        "    RDataClass = \"MultiAssayExperiment\", # class that is returned from hub",
        "    DispatchClass = \"FilePath\", # format saved on disk; FilePath only returns file location",
        "    # Location_Prefix = \"\", # SKIP if data stored in the Bioconductor AWS S3",
        paste0("    RDataPath = \"scMultiome/", dataset, "\""),
        ")",
        "",
        "# dataset manifest metadata",
        "# see ?listDatasets for details",
        paste0("manifest.", dataset, " <- data.frame("),
        "    Title = \"<YOUR DATA TITLE>\",",
        "    Species = \"<SPECIES NAME>\", # e.g. \"Homo sapiens\"",
        "    Type = \"<CELL CULTURE OR TISSUE>\",",
        "    Multiome = \"<PAIRED OR UNPAIRED>\",",
        "    DiskSize = \"<H5 FILE SIZE ON DISK>\",",
        "    MemorySize = \"<MAE SIZE IN MEMORY>\",",
        "    Accessor = \"<NAME OF ACCESSOR FUNCTION>\",",
        paste0("    Version = \"", Sys.Date(), "\""),
        ")"
    ))

    if ("tools:rstudio" %in% search() && requireNamespace("rstudioapi", quietly = TRUE)) {
        rstudioapi::navigateToFile(fileName)
    } else {
        message(fileName, " created")
    }

    return(invisible(TRUE))
}
