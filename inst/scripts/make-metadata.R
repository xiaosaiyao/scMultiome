
# this script creates package metadata
# ExperimentHub DEMANDS that the metadata be in a data frame that is written to a csv file (last line)
# see https://bioconductor.org/packages/3.15/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html#overview
# the FINAL STEPS section should only be altered after careful consideration


### PROCEDURE
# create a 1-row data frame for every dataset and name it "metadata.<3-DIGIT_INDEX>"
# see ?ExperimentHubData::makeExperimentHubMetadata for field details

### DESCRIBE FIRST DATASET
metadata.001 <- data.frame(
    Title = "prostateENZ",
    Description = "Multiome data (scATAC and scRNAseq) for enzalutamide-treated LNCaP cells. Measured (unpaired) with Illumina NextSeq 500 (GPL18573). Results stored in GEO, acc.nos: GSE168667 and GSE168668.",
    BiocVersion = "3.16",
    Genome = "hg38",
    SourceType = "tar.gz",
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168667&format=file, https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168668&format=file", # optional
    SourceVersion = "Sep 08 2021", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = NA, # ?
    DataProvider = "Tampere University",
    Maintainer = "Xiaosai Yao <yaox19@gene.com>",
    RDataClass = "MultiAssayExperiment",
    DispatchClass = "H5File",
    # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
    RDataPath = "", # TBD
    stringsAsFactors = FALSE
)




### DESCRIBE FURTHER DATASETS

# more objects documented here




### FINAL STEPS

# collate metadata
metadata.all <- ls(envir = .GlobalEnv, pattern = "metadata\\.\\d{3}")
metadata <- do.call(rbind, mget(metadata.all))

# write metadata to file - DO NOT ALTER
utils::write.csv(metadata,
                 file.path(system.file("extdata", package = "scMultiome"), "metadata.csv"),
                 row.names = FALSE)
