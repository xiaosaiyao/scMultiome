#'
#' create accessor function file
#'
#' Creates a template for accessor function and data set help.
#'
#' @param dataset name of data set as character string
#'
#' @return Invisible TRUE.
#' Creates \code{R/<dataset>.Rmd} file and navigates to it, if possible.
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
        "        retrieve(metadata, experiments, verbose = FALSE)",
        "    }"
        ))

    if ("tools:rstudio" %in% search() && requireNamespace("rstudioapi", quietly = TRUE)) {
        rstudioapi::navigateToFile(fileName)
    } else {
        message(fileName, " created")
    }

    return(invisible(TRUE))
}
