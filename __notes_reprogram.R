
### SAVE REPROGRAM DATA ###

# attach scMultiome
## attaches SingleCellExperiment and MultiAssayExperiment as well
devtools::load_all() # library(scMultiome)

# load data set
reprogram <- dsassembly::getDataset("DS000013080")
# isolate Matrix
GEM <- reprogram[["GeneExpressionMatrix"]]
# add row names
rownames(GEM) <- rowData(GEM)[["name"]]
# rename assay
assay(GEM, "logcounts") <- assay(GEM, "counts")

# load raw counts
seRNA <- readRDS("/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/seRNA.rds")
# subset to fit Matrix
seRNAc <- as(as.matrix(assay(seRNA, "counts"))[rownames(GEM), colnames(GEM)], "dgCMatrix")

# add assay to Matrix
assay(GEM, "counts") <- seRNAc

# rebuild MAE
expList <- experiments(reprogram)
expList[["GeneExpressionMatrix"]] <- GEM
mae <- MultiAssayExperiment(experiments = expList)

# save and test MAE
saveMAE(mae, "inst/extdata/reprogram_Xie.h5")
hm <- testFile("inst/extdata/reprogram_Xie.h5")










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
file.size("inst/extdata/reprogram_Xie.h5") / 1024 / 1024

### DATA IS SAVED AND RELOAD WAS TESTED ###




# data set name: reprogram_Xie

# # add author - DONE
# desc::desc_add_author("Shiqi", "Xie", "xie.shiqi@gene.com", "ctb")

# # add metadata and documentation - TODO
# makeMakeMetadata("reprogram_Xie")
# makeMakeData("reprogram_Xie")
# makeR("reprogram_Xie")



# data set sizes to put in manifest file:
# in memory: 660 MB - 2.2 GB
# on disk: 421 MB

# object.size(HM) / 1024 / 1024
# file.size("inst/extdata/reprogram_Xie.h5") / 1024 / 1024


