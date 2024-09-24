
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.TEADi_resistance <-
  data.frame(
    Title = "Resistance of TEAD inhibitor to drug",
    Description = "Cooperation between the Hippo and MAPK pathway activation drives acquired resistance to TEAD inhibition",
    BiocVersion = "3.20",
    Genome = "hg38",
    SourceType = "HDF5", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247442", # refers to raw data
    SourceVersion = "2022-05-19", # no commas!
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    TaxonomyId = "9606", # e.g. "9606"
    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE
    DataProvider = "Genentech, Inc.", # refers to raw data
    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer
    RDataClass = "MultiAssayExperiment", # class that is returned from hub
    DispatchClass = "FilePath", # format saved on disk; FilePath only returns file location
    # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
    RDataPath = "scMultiome/TEADi_resistance"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.TEADi_resistance <-
  data.frame(
    Call = "TEADi_resistance()",
    Author = "Xiaosai Yao",
    Title = "Resistance of TEAD inhibitor to drug",
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    Lineage = "H226 cell line",
    CellNumber = "<NUMBER OF CELLS>",
    Multiome = "paired",
    DiskSize = "403.9 MB",
    Version = "2024-09-17"
)
