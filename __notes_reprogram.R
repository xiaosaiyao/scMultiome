
### SAVE REPROGRAM DATA ###

# download
mae <- dsassembly::getDataset("DS000013080")

# strip non-public class
## identify offending experiment
tmind <- which(names(mae) == "TileMatrix500")
## build another object with SCEs only
maelim <- MultiAssayExperiment(experiments = c(
    TileMatrix500 = as(experiments(mae)[[tmind]], "SingleCellExperiment"),
    as.list(experiments(mae)[-tmind])
))

# save
saveMAE(maelim, "inst/extdata/reprogram.h5", verbose = TRUE)


# data set sizes to put in manifest file:
# in memory: 660 MB - 2.2 GB
# on disk: 420 MB
