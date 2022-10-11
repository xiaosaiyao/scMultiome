
#' documentation templates
#'
#' Create templates for data set documentation.
#'
#' Functions to facilitate documenting a new data set.
#' Each function creates a file and attempts to open it in RStudio for editing.
#'
#' \code{makeMakeData} creates an Rmarkdown report called \code{inst/scripts/make-data-<dataset>.Rmd}.\cr
#' \code{makeMakeMetadata} creates an R script called \code{inst/scripts/make-metadata-<dataset>.R}.\cr
#' \code{makeR} creates an R file called \code{R/<dataset>.Rmd}.\cr
#'
#' @param dataset name of data set as character string
#'
#' @return
#' All functions return TRUE invisibly.
#'
#' @name templates
#'
#' @rdname templates
#'
makeMakeData <- function(dataset) {
    checkmate::assertString(dataset)

    fileName <- file.path(system.file("scripts", package = "scMultiome"),
                          paste0("make-data-", dataset, ".Rmd"))
    file.create(fileName)
    writeLines(con = fileName, text = c(
        "---",
        "title: \"Your Data Set Title\"",
        "author: \"Your Name\"",
        paste("date:", Sys.Date()),
        "output:",
        "  BiocStyle::html_document:",
        "    titlecaps: false",
        "    toc_float: true",
        "---",
        "",
        "```{r setup, include = FALSE}",
        "knitr::opts_chunk$set(echo = TRUE, eval = FALSE)",
        "```",
        ""))

    if ("tools:rstudio" %in% search() && requireNamespace("rstudioapi", quietly = TRUE)) {
        rstudioapi::navigateToFile(fileName)
    } else {
        message(fileName, " created")
    }

    return(invisible(TRUE))
}



#' @rdname templates
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
        paste0("metadata.", dataset, " <-"),
        "  data.frame(",
        "    Title = \"<YOUR DATA TITLE>\",",
        "    Description = \"<YOUR DATA TITLE>\",",
        paste0("    BiocVersion = \"", as.character(utils::packageVersion("BiocVersion")[,seq_len(2)]), "\","),
        "    Genome = \"<GENOME VERSION USED IN ANALYSIS>\",",
        "    SourceType = \"<SOURCE FILE TYPE>\", # refers to raw data",
        "    SourceUrl = \"<SOURCE FILE LOCATION>\", # refers to raw data",
        "    SourceVersion = \"<VERSION NUMEBR OR DATE>\", # no commas!",
        "    Species = \"<SPECIES NAME>\", # e.g. \"Homo sapiens\"",
        "    TaxonomyId = \"<SPECIES ID>\", # e.g. \"9606\"",
        "    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE",
        "    DataProvider = \"<SOURCE DATA PROVIDER (INSTITUTION)>\", # refers to raw data",
        "    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer",
        "    RDataClass = \"MultiAssayExperiment\", # class that is returned from hub",
        "    DispatchClass = \"FilePath\", # format saved on disk; FilePath only returns file location",
        "    # Location_Prefix = \"\", # SKIP if data stored in the Bioconductor AWS S3",
        paste0("    RDataPath = \"scMultiome/", dataset, "\""),
        ")",
        "",
        "# dataset manifest metadata",
        "# see ?listDatasets for details",
        paste0("manifest.", dataset, " <-"),
        paste0("    Call = \"", dataset, "()\","),
        "    Author = \"<DATA SET AUTHOR>\",",
        "    Title = \"<YOUR DATA TITLE>\",",
        "    Species = \"<SPECIES NAME>\", # e.g. \"Homo sapiens\"",
        "    Lineage = \"<TISSUE OR ORGAN>\",",
        "    Cell_Num = \"<NUMBER OF CELLS>\",",
        "    Multiome = \"<PAIRED OR UNPAIRED>\",",
        "    DiskSize = \"<H5 FILE SIZE ON DISK>\",",
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



#' @rdname templates
#'
makeR <- function(dataset) {
    checkmate::assertString(dataset)

    fileName <- file.path(system.file("R", package = "scMultiome"),
                          paste0(dataset, ".R"))
    file.create(fileName)

    writeLines(con = fileName, text = c(
        "#'",
        "#' <YOUR DATA SET TITLE>",
        "#'",
        "#' <YOUR DATA SET DESCRIPTION",
        "#'",
        "#' @inheritParams prostateENZ",
        "#'",
        "#' @inherit prostateENZ return",
        "#'",
        "#' @format",
        "#' <YOUR DATA FORMAT>",
        "#' \\code{MultiAssayExperiment} obtained from an \\code{ArchR} project.",
        "#' Annotated with the <YOUR GENOME> genome build.",
        "#' Contains the following experiments:",
        "#' \\itemize{",
        "#'   \\item{}",
        "#' }",
        "#'",
        "#' @references",
        "#' <YOUR DATA REFERENCE>",
        "#'",
        "#' @inheritSection prostateENZ Data storage and access",
        "#'",
        "#' @section Data preparation:",
        paste0("#' ```{r child = system.file(\"scripts\", \"make-data-", dataset, ".Rmd\", package = \"scMultiome\")}"),
        "#' ```",
        "#'",
        "#' @examples",
        paste0("#' ", dataset, "()"),
        "#'",
        "#' @export",
        "#'",
        "#'",
        paste(dataset, "<-"),
        "    function(metadata = FALSE,",
        "             experiments = c(\"TileMatrix500\",",
        "                             \"GeneScoreMatrix\",",
        "                             \"GeneIntegrationMatrix\",",
        "                             \"PeakMatrix\",",
        "                             \"MotifMatrix\")) {",
        "        checkmate::assertFlag(metadata)",
        "        experiments <- match.arg(experiments, several.ok = TRUE)",
        "",
        paste0("retrieve(\"", dataset, "\", metadata, experiments, verbose = FALSE)"),
        "    }",
        "",
        "# place methods to convert your custom class to SingleCellExperiment and vice versa here",
        "# skip if conversion can be done with `as`",
        "# see ?convertSCE and the Adding Data Sets vignette (Custom Classes) for details"
    ))

    if ("tools:rstudio" %in% search() && requireNamespace("rstudioapi", quietly = TRUE)) {
        rstudioapi::navigateToFile(fileName)
    } else {
        message(fileName, " created")
    }

    return(invisible(TRUE))
}
