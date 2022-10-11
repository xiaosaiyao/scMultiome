
### SAVEDATA ###

# attach scMultiome
## attaches SingleCellExperiment and MultiAssayExperiment as well
devtools::load_all() # library(scMultiome)

# load data set
mae <- dsassembly::getDataset("DS000013664")

# save and test MAE
saveMAE(mae, "inst/extdata/hematopoiesis.h5")
hm <- testFile("inst/extdata/hematopoiesis.h5")










# materialize whole object in memory to measure size
exps <- as.list(experiments(hm))
exps %>% lapply(assays) %>% sapply(length)

library(magritttr)

assay(exps[[1]], 1) %<>% as("dgCMatrix")
assay(exps[[2]], 1) %<>% as("dgCMatrix")
assay(exps[[3]], 1) %<>% as("dgCMatrix")
assay(exps[[3]], 2) %<>% as("dgCMatrix")
assay(exps[[4]], 1) %<>% as("dgCMatrix")
assay(exps[[5]], 1) %<>% as("dgCMatrix")
assay(exps[[5]], 2) %<>% as("dgCMatrix")
assay(exps[[6]], 1) %<>% as("dgCMatrix")


HM <- MultiAssayExperiment(experiments = exps)

object.size(hm) / 1024 / 1024
object.size(HM) / 1024 / 1024
file.size("inst/extdata/reprogramSeq.h5") / 1024 / 1024

### DATA IS SAVED AND RELOAD WAS TESTED ###




# data set name: reprogram_Xie

# # add author - DONE
# desc::desc_add_author("Shiqi", "Xie", "xie.shiqi@gene.com", "ctb")

# # add metadata and documentation - TODO
# makeMakeMetadata("reprogramSeq")
# makeMakeData("reprogramSeq")
# makeR("reprogramSeq")



# data set sizes to put in manifest file:
# in memory: 660 MB - 2.2 GB
# on disk: 421 MB

# object.size(HM) / 1024 / 1024
# file.size("inst/extdata/reprogramSeq.h5") / 1024 / 1024


