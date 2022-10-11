
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.colonHealthy <- data.frame(
    Title = "Single-cell analysis of samples from healthy human colon",
    Description = "ATACseq and RNAseq data obtained by the colon tissues analysis. Samples were collected from adult human donors.",
    BiocVersion = "3.16",
    Genome = "hg38",
    SourceType = "tar.gz", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE165659&format=file", # refers to raw data; optional
    SourceVersion = "Feb 16 2021", # no commas!
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    TaxonomyId = "9606", # e.g. "9606"
    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE.
    DataProvider = "Stanford University", # refers to raw data
    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer
    RDataClass = "MultiAssayExperiment", # class that is returned from hub
    DispatchClass = "FilePath", # format saved on disk; FilePath only returns file location
    # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
    RDataPath = "scMultiome/colonHealthy"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.colonHealthy <- data.frame(
    Call = "colonHealthy",
    Author = "Zhang",
    Title = "Healthy colon",
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    Lineage = "Colon",
    Cell_Num = "59231",
    Multiome = "unpaired",
    DiskSize = "6.7 GB",
    Version = "2022-09-21"
)
