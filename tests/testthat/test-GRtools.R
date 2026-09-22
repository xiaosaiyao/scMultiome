
# # NOTE:
# In the wild the reconverted GRanges object may not be identical to the original
# as factor levels will be arranged alphabetically.

gr0 <- GenomicRanges::GRanges(
    "seqnames" = S4Vectors::Rle(c("chr1", "chr1", "chr2", "chr3"), c(1, 3, 2, 4)),
    "ranges" = IRanges::IRanges(1:10, width = 10:1, names = head(letters, 10)),
    "strand" = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    "score" = 1:10,
    "GC" = seq(1, 0, length = 10))


# prepare for storage
gr0df <- storeGR(gr0)

# restore from storage
GR0 <- restoreGR(gr0df)

test_that("conversion of GR to data.frame", {
    expect_s3_class(gr0df, "data.frame")
    expect_false(any(vapply(gr0df, is.factor, logical(1L))))
})

test_that("reconversion of data.frame to GR", {
    expect_s4_class(GR0, "GRanges")
    expect_identical(GR0, gr0)
})
