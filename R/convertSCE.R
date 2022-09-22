
#' convert SingleCellExperiment to custom class
#'
#' Converts an SCE object to a custom class, if a method is defined.
#'
#' This function attempts to convert a SingleCellExperiment object to another class.
#' A class-specific method for \code{newClass} is looked for in this package's namespace.
#' The method name \strong{must} follow the format: "convertTo\code{<newClass>}".
#' If no such method is found, the default is to use `as`.
#' If that fails, the original \code{SingleCellExperiment} will be returned.
#'
#' @param sce a \code{SingleCellExperiment}
#' @param newClass character string specifying the target class
#'
#' @return Object of class \code{newClass}. If conversion fails, the original \code{sce}.
#'
#' @seealso \code{loadExp}
#'
convertSCE <- function(sce, newClass) {
    checkmate::assertClass("SingleCellExperiment")
    checkmate::assertString(newClass)

    # find contents of scMultiome namespace
    cont <- ls(pattern = "^convertTo.+", envir = getNamespace("scMultiome"))
    # find method for newClass
    meth <- grep(newClass, cont, value = TRUE)

    if (length(meth) > 1L) {
        stop("multiple conversion methods found for class ", newClass)
    } else if (length(meth) == 1L) {
        meth <- get(meth, envir = getNamespace("scMultiome"))
        ans <- meth(sce)
    } else {
        # generic method (default)
        ans <- tryCatch(methods::as(sce, newClass),
                        error = function(e) {
                            warning("conversion to class ", newClass,
                                    " failed, returning SingleCellExperiment")
                            return(sce)
                        })
    }

    return(ans)
}
