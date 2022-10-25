
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.reprogramSeq <-
  data.frame(
    Title = "Reprogram-seq of LNCaP cells",
    Description = "scMultiome data of LNCaP infected with FOXA1, NKX2-1, GATA6",
    BiocVersion = "3.16",
    Genome = "hg38",
    SourceType = "HDF5", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/", # refers to raw data
    SourceVersion = "2022-10-06", # no commas!
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    TaxonomyId = "9606", # e.g. "9606"
    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer
    RDataClass = "MultiAssayExperiment", # class that is returned from hub
    DispatchClass = "FilePath", # format saved on disk; FilePath only returns file location
    # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
    RDataPath = "scMultiome/reprogramSeq"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.reprogramSeq <-
    data.frame(
        Call = "reprogramSeq()",
        Author = "Xie",
        Title = "Reprogram-seq of LNCaP cells",
        Species = "Homo sapiens", # e.g. "Homo sapiens"
        Lineage = "Prostate",
        CellNumber = "3903",
        Multiome = "paired",
        DiskSize = "430 MB",
        Version = "2022-10-04"
)

