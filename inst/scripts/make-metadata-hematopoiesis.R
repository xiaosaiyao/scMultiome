
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.hematopoiesis <-
  data.frame(
    Title = "scATAC-seq and unpaired scRNA-seq of hematopoetic cells",
    Description = "Example scATAC-seq data of hematopoietic cells included in ArchR package was integrated with scRNAseq. ScATAC-seq data was obtained from GSE139369 and scRNA-seq obtained from https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds",
    BiocVersion = "3.16",
    Genome = "hg19",
    SourceType = "fragment", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369", # refers to raw data
    SourceVersion = "2019-10-29", # no commas!
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    TaxonomyId = "9606", # e.g. "9606"
    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE
    DataProvider = "Stanford", # refers to raw data
    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer
    RDataClass = "MultiAssayExperiment", # class that is returned from hub
    DispatchClass = "FilePath", # format saved on disk; FilePath only returns file location
    # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
    RDataPath = "scMultiome/hematopoiesis"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.hematopoiesis <-
  data.frame(
    Call = "hematopoiesis()",
    Author = "Granje",
    Title = "scATAC-seq and unpaired scRNA-seq of hematopoetic cells",
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    Lineage = "blood",
    Cell_Num = "10250",
    Multiome = "unpaired",
    DiskSize = "1.0 GB",
    Version = "2019-10-29"
)
