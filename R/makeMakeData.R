#'
#' create make-data file
#'
#' Creates a template for data set history file.
#'
#' @param dataset name of data set as character string
#'
#' @return Invisible TRUE.
#' Creates \code{inst/scripts/make-data-<dataset>.Rmd} file and navigates to it, if possible.
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
