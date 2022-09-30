
#' convert custom experiments classes to and from SingleCellExperiment
#'
#' Converts a custom class experiment object to SCE and \emph{vice versa},
#' if methods are defined.
#'
#' Multiome data sets are saved as MultiAssayExperiments composed of SingleCellExperiments.
#' If a data set contains other classes, they must be converted to SCE before saving
#' and reconverted to those classes after loading. Some packages provide handlers
#' for such conversions, either utilizing the \code{methods::as} generic,
#' or defining their own methods.
#' These functions provide a catch-all mechanism to do the conversions automatically
#' during saving and loading experiments.
#' Authors adding new data sets may choose to do it manually.
#'
#' \code{convergeSCE} is called by \code{saveExp}.
#' It attempts to convert an experiment of another class to SingleCellExperiment.
#' A class-specific method is searched for in this package's namespace.
#' This method \strong{must} follow the naming format: "convertFrom\code{<customClass>}".
#' If no such method is found, the default is to use `as`.
#' If that fails, an error is signaled.
#'
#' \code{convertSCE} is called by \code{loadExp}.
#' It attempts to convert a SingleCellExperiment to another class.
#' A class-specific method is search for in this package's namespace.
#' This method \strong{must} follow the naming format: "convertTo\code{<customClass>}".
#' If no such method is found, the default is to use `as`.
#' If that fails, the original \code{SingleCellExperiment} will be returned.
#'
#' @param exp an experiment other than a \code{SingleCellExperiment}
#' @param sce a \code{SingleCellExperiment}
#' @param customClass character string specifying the custom class
#'
#' @return
#' \code{convergeSCE} returns an object of class \code{SingleCellExperiment}
#' or an error if conversion fails.
#' \code{convertSCE} returns an object of class \code{customClass}
#' or the original \code{sce} if conversion fails.
#'
#' @seealso \code{saveExp}, \code{loadExp}
#'
#' @name experimentClass
#'
#' @rdname experimentClass
#' @keywords internal
#'
convergeSCE <- function(exp, customClass) {
    checkmate::assertString(customClass)

    # find contents of scMultiome namespace
    cont <- ls(pattern = "^convertFrom.+", envir = getNamespace("scMultiome"))
    # find method for customClass
    meth <- grep(customClass, cont, value = TRUE)

    if (length(meth) > 1L) {
        stop("multiple conversion methods found for class ", customClass)
    } else if (length(meth) == 1L) {
        meth <- get(meth, envir = getNamespace("scMultiome"))
        ans <- meth(exp)
    } else {
        # generic method (default)
        ans <- tryCatch(methods::as(exp, customClass),
                        error = function(e) {
                            stop("Conversion from class ", customClass,
                                 " to SingleCellExperiment failed. ",
                                 "Try converting manually.",
                                 call. = FALSE)
                        })
    }
    return(ans)
}



#' @rdname experimentClass
#' @keywords internal
#'
convertSCE <- function(sce, customClass) {
    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertString(customClass)

    # find contents of scMultiome namespace
    cont <- ls(pattern = "^convertTo.+", envir = getNamespace("scMultiome"))
    # find method for customClass
    meth <- grep(customClass, cont, value = TRUE)

    if (length(meth) > 1L) {
        stop("multiple conversion methods found for class ", customClass)
    } else if (length(meth) == 1L) {
        meth <- get(meth, envir = getNamespace("scMultiome"))
        ans <- meth(sce)
    } else {
        # generic method (default)
        ans <- tryCatch(methods::as(sce, customClass),
                        error = function(e) {
                            warning("Conversion to class ", customClass, " failed. ",
                                    "Returning SingleCellExperiment")
                            return(sce)
                        })
    }

    return(ans)
}
