
# This script creates package metadata displayed in ExperimentHub.
# It also creates the data set manifest used by the listDatasets function.
# ExperimentHub DEMANDS that the metadata be in a data frame that is written to a csv file.
# See https://bioconductor.org/packages/3.15/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html#overview for details.
# Describe your data set in a file called /inst/scripts/make-metadata-<DATASET_NAME>.R
#   by creating two 1-row data frames called metadata.<DATASET_NAME> and manifest.<DATASET_NAME>.
# See ?ExperimentHubData::makeExperimentHubMetadata and ?listDatasets for field details.

# When your file is ready, run this script.

# discover metadata files
files <- list.files(path = system.file("scripts", package = "scMultiome"),
                    full.names = TRUE,
                    pattern = "make-metadata-.+\\.R")
# create metadata items
lapply(files, source)


# collate manifest
manifest.all <- ls(envir = .GlobalEnv, pattern = "^manifest\\..+")
manifest <- do.call(rbind, mget(manifest.all))

# write manifest to file - DO NOT ALTER
utils::write.csv(manifest,
                 file.path(system.file("extdata", package = "scMultiome"), "manifest.csv"),
                 row.names = FALSE)


# collate metadata
metadata.all <- ls(envir = .GlobalEnv, pattern = "^metadata\\..+")
metadata <- do.call(rbind, mget(metadata.all))

# write metadata to file - DO NOT ALTER
fileLoc <- system.file("extdata", package = "scMultiome")
fileName <- "metadata.csv"
## add current package version to previous metadata file
if (file.exists(file.path(fileLoc, fileName))) {
    oldVersion <- paste0("_v", paste(strsplit(as.character(desc::desc_get("Version")), "\\.")[[1]][1:3], collapse = ""))
    oldFileName <- paste0("metadata", oldVersion, ".csv")
    file.rename(file.path(fileLoc, "metadata.csv"), file.path(fileLoc, oldFileName))
}
## wirte new metadata
utils::write.csv(metadata, file.path(fileLoc, fileName), row.names = FALSE)



# create data set list for package man page
makeDataSetList(metadata)



# clean up to prevent corruption on re-run
rm(list = ls(pattern = "^manifest.?"))
rm(list = ls(pattern = "^metadata.?"))
rm(files, fileLoc, fileName, oldFileName, oldVersion)
