
#' manipulate GeneRanges
#'
#' Prepare and recover \code{GRanges} objects for and after storing.
#'
#' \code{GRanges} objects, which can be encountered in \code{rowRanges} slots of
#' \code{SingleCellExperiment}s, are stored as data frames (of type compound).
#'
#' \code{storeGR} converts \code{GRanges} to a data frames converts factors to characters.
#' \code{restoreGR} resets data types in the basic columns and re-instantiates \code{GRanges}.
#'
#' @name RGtools
#'
#' @param x object of class \code{GRanges}
#' @param df a \code{data.frame}
#'
#' @return
#' \code{storeGR} returns a data frame, \code{restoreGR} returns a \code{GRanges} object.
#'

#' @rdname GRtools
#'
#' @keywords internal
#'
storeGR <- function(x) {
    checkmate::assertClass(x, "GRanges")
    ans <- as.data.frame(x)
    rownames(ans) <- names(x)
    ans[vapply(ans, is.factor, logical(1L))] <- lapply(Filter(is.factor, ans), as.character)

    return(ans)
}

#' @rdname GRtools
#'
#' @keywords internal
#'
restoreGR <- function(df) {
    checkmate::assertDataFrame(df)
    basic <- c("seqnames", "start", "end", "width", "strand")
    df[basic] <- lapply(df[basic], function(x) methods::as(x, typeof(x)))
    ans <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)

    return(ans)
}

