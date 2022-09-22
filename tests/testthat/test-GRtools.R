
# # NOTE:
# In the wild the reconverted GRanges object may not be identical to the original
# as factor levels will be arranged alphabetically.

gr0 <- GenomicRanges::GRanges(
    Rle(c("chr1", "chr1", "chr2", "chr3"), c(1, 3, 2, 4)),
    IRanges(1:10, width = 10:1)
)
strand(gr0) <- Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2))


# prepare for storage
gr0df <- storeGR(gr0)

# restore from storage
GR0 <- restoreGR(gr0df)

test_that("conversion of GR to data.frame", {
    expect_s3_class(gr0df, "data.frame")
    expect_false(any(sapply(gr0df, is.factor)))
})

test_that("reconversion of data.frame to GR", {
    expect_s4_class(GR0, "GRanges")
    expect_identical(GR0, gr0)
})
