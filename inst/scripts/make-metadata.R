
# This script creates package metadata displayed in ExperimentHub.
# ExperimentHub DEMANDS that the metadata be in a data frame that is written to a csv file.
# See https://bioconductor.org/packages/3.15/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html#overview for details.
# Describe your data set in a file called /inst/scripts/make-metadata-<DATASET_NAME>.R.
# See ?ExperimentHubData::makeExperimentHubMetadata for field details.

# This script also creates the data set manifest used by the listDatasets function.
# See ?listDatasets for details.

# When your file is ready, run this script.

# discover metadata files
files <- list.files(path = system.file("scripts", package = "scMultiome"),
                    full.names = TRUE,
                    pattern = "make-metadata-.+\\.R")
# create metadata items
lapply(files, source)


# collate manifest
manifest.all <- ls(envir = .GlobalEnv, pattern = "manifest\\.\\w+")
manifest <- do.call(rbind, mget(manifest.all))

# write manifest to file - DO NOT ALTER
utils::write.csv(manifest,
                 file.path(system.file("extdata", package = "scMultiome"), "manifest.csv"),
                 row.names = FALSE)


# collate metadata
metadata.all <- ls(envir = .GlobalEnv, pattern = "metadata\\.\\w+")
metadata <- do.call(rbind, mget(metadata.all))

# write metadata to file - DO NOT ALTER
utils::write.csv(metadata,
                 file.path(system.file("extdata", package = "scMultiome"), "metadata.csv"),
                 row.names = FALSE)
