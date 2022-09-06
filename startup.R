

### CLEAN UP

unlink(list.files("R", full.names = TRUE))
unlink(list.files("man", full.names = TRUE))
unlink("NAMESPACE")

### ADD FILES
dir.create("inst/scripts", recursive = TRUE)
dir.create("inst/extdata", recursive = TRUE)
file.create("./inst/scripts/make-data.R")
file.create("./inst/scripts/make-metadata.R")

# MODIFY DESCRIPTION:
desc::desc_del(c("Author", "Maintainer"))
desc::desc_set_version("0.0.0.9000")
usethis::use_package("SummarizedExperiment", type = "Depends")
usethis::use_package("SingleCellExperiment", type = "Depends")
usethis::use_package("MultiAssayExperiment", type = "Depends")
usethis::use_package("dsassembly", type = "Depends")
usethis::use_package("dsdb.plus", type = "Depends")
usethis::use_package("ExperimentHub", type = "Depends")
rstudioapi::navigateToFile("DESCRIPTION")
# edit fields
desc::desc_set("Package", "scMultiome")
desc::desc_set("Title", "Collection of Public Single-Cell Multiome (scATAC + scRNAseq) Datasets")
desc::desc_set("Description", "Single cell multiome data, containing chromatin accessibility (scATAC-seq) and gene expression (scRNA-seq) information analyzed with the ArchR package and presented as MultiAssayExperiment objects.")
# add authors with desc::desc_add_author (one author and one role at a time)
desc::desc_add_author("Xiaosai", "Yao", "xyao19@gene.com", "cre")
desc::desc_add_author("Aleksander", "Chlebowski", "aleksander.chlebowski@contractors.roche.com", "aut")
# add biocViews field (required by Bioconductor)
desc::desc_set("biocViews", "ExperimentHub, SingleCellData, ExpressionData, Homo_sapiens_Data, CellCulture, Tissue, GEO")
# enable Rmarkdown support for help files (allows inserting unevaluated code chunks with ```R)
desc::desc_set("Roxygen", "list(markdown = TRUE)")

desc::desc_normalize()

### ADD DOCUMENTATION
usethis::use_package_doc()
# edit package documentation

# wrap up
devtools::document()




### CREATE VIGNETTE
usethis::use_vignette(name = "multiome", title = "Analyzing Epiregulon Data")


### NOTES

# # change output in vignette(s)
# output:
#     rmarkdown::html_vignette:
#     toc: true
# number_section: true
# self_contained: true
# titlecaps: true



### ADDING NEW DATASETS
# 0. open the package in developer mode and run devtools::load_all()
#   this is necessary in order for the scripts below to write to files properly
# 1. add metadata for that object to inst/scripts/make-metadata.R; run the script
#   (OR use HubPuB::hubMetadata and HubPuB::add_resource)
# 2. run ExperimentHubData::makeExperimentHubMetadata(file.path("/gstore/home", Sys.getenv("USER"), "scMultiome"))
#    and ensure it passes without errors
# 3. add the dataset to the manifest in inst/scripts/make-manifest.R; run the script
# 4. create an R file to document you dataset; copy and adapt the accessor function from R/prostateENZ.R;
#    the default value of the "experiments" argument must reflect experiment names in your dataset
# 5. add your dataset to the list in R/scMultiome-package.R
# 6. run devtools::document()
# 7. add your dataset to the R/data directory with saveMAE
# 8. test that the accessor function for your dataset works
# 9. run inst/scripts/upload-data.R
#
# 10. update package version with desc::desc_bump_version("minor")
# 11. reorder the Depends field in the DESCRIPTION file so that ExperimentHub is the last package listed there
#
