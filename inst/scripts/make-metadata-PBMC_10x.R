
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.dataset <-
  data.frame(
    Title = "PBMC Data Set",
    Description = "PBMC from a Healthy Donor - Granulocytes Removed Through Cell Sorting (10k)",
    BiocVersion = "3.20",
    Genome = "hg38",
    SourceType = "HDF5",
    SourceUrl = "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/",
    SourceVersion = "2021-05-03",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "10X Genomics",
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "MultiAssayExperiment",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/dataset"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.dataset <-
  data.frame(
    Call = "PBMC_10x()",
    Author = "10X Genomics",
    Title = "PBMC Data Set",
    Species = "Homo sapiens",
    Lineage = "Blood",
    CellNumber = "11898",
    Multiome = "paired",
    DiskSize = "1.2 GB",
    Version = "2024-09-17"
)
