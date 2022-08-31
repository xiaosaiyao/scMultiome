
# this script creates the dataset manifest table

### PROCEDURE
# create a 1-row data frame for every dataset and name it "manifest.<3-DIGIT_INDEX>"


### DESCRIBE FIRST DATASET
manifest.001 <- data.frame(
    Species = "Homo sapiens",
    Type = "cell culture",
    Multiome = "unpaired",
    DiskSize = "", # TBD
    MemorySize = "", #TBD
    Accessor = "prostateENZ",
    stringsAsFactors = FALSE
)


### DESCRIBE FURTHER DATASETS

# more objects documented here



### FINAL STEPS

# collate metadata
manifest.all <- ls(envir = .GlobalEnv, pattern = "manifest\\.\\d{3}")
manifest <- do.call(rbind, mget(manifest.all))

# write metadata to file - DO NOT ALTER
utils::write.csv(manifest,
                 file.path(system.file("extdata", package = "scMultiome"), "manifest.csv"),
                 row.names = FALSE)

