

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
usethis::use_package("ExperimentHub", type = "Depends")
usethis::use_package("MultiAssayExperiment", type = "Depends")
rstudioapi::navigateToFile("DESCRIPTION")
# edit fields
desc::desc_set("Package", "scMultiome")
desc::desc_set("Title", "Collection of Public Single-Cell Multiome (scATAC + scRNAseq) Datasets")
desc::desc_set("Description", "Single cell multiome data, containing chromatin accessibility (scATAC-seq) and gene expression (scRNA-seq) information analyzed with the ArchR package and presented as MultiAssayExperiment objects.")
# add authors with desc::desc_add_author (one author and one role at a time)
desc::desc_add_author("Aleksander", "Chlebowski", "aleksander.chlebowski@contractors.roche.com", "aut")
desc::desc_add_author("Xiaosai", "Yao", "xyao19@gene.com", "cre")
# add biocViews field
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
# 2. make sure ExperimentHubData::addResources(sprintf("/gstore/home/%s/scMulitome", Sys.getenv("USER"))) runs without errors
# 3. add the dataset to the manifest in inst/scripts/make-manifest.R; run the script
# 4. create an R file to document you dataset
# 5. add your dataset to the list in R/scMultiome-package.R
# 6. run devtools::document() (you may need to comment out the .onLoad function in R/zzz.R for this to pass)
# 7. update package version with desc::desc_bump_version("minor")
#
# reorder the Depends field in the DESCRIPTION file so that ExperimentHub is the last package listed there
#
