
#' pad numbers in strings with zeros
#'
#' Insert "0" characters to unify number of characters in character vector.
#'
#' Given a character vector composed of letters and digits,
#' add zeros at the beginning of the numeric part so that all strings are of the same length.
#'
#' @param x a character vector
#'
#' @return A character where all elements have the same number of characters.
#'
#' @examples
#' padZeros(paste0("C", 1:5))
#' padZeros(paste0("C", 1:10))
#'
#' @export
#'
padZeros <- function(x) {
    checkmate::assertCharacter(x)

    ans1 <- sub("(\\D+)(\\d+)", "\\1", x)
    ans2 <- sub("(\\D+)(\\d+)", "\\2", x)
    n <- max(nchar(ans2))
    fmt <- paste0("%s%0", n, "i")
    ans <- sprintf(fmt, ans1, as.integer(ans2))

    return(ans)
}
