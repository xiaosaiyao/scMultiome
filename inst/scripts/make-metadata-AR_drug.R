
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.AR_drug <-
  data.frame(
    Title = "Response of prostate cancer cells to drug treatment",
    Description = "Gene expression and ATACseq data from 6 cell lines (LNCaP, VCaP, MDA, 22Rv1, DU145, H660) after treatment with AR-targeting drugs",
    BiocVersion = "3.20",
    Genome = "hg38",
    SourceType = "HDF5",
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251977",
    SourceVersion = "Sep 09 2024",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech, Inc.",
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "MultiAssayExperiment",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/AR_drug.h5"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.AR_drug <-
  data.frame(
    Call = "AR_drug()",
    Author = "Xiaosai Yao",
    Title = "Response of prostate cancer cells to drug treatment",
    Species = "Homo sapiens",
    Lineage = "Prostate",
    CellNumber = "23118",
    Multiome = "paired",
    DiskSize = "3.3 GB",
    Version = "2024-09-09"
)
