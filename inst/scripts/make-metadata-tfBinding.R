
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.tfBinding.hg38 <- data.frame(
    Title = "TF Binding Info Hg38",
    Description = "Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE",
    BiocVersion = "3.16",
    Genome = "hg38",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfbinding"
)
metadata.tfBinding.hg19 <- data.frame(
    Title = "TF Binding Info Hg19",
    Description = "Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE",
    BiocVersion = "3.16",
    Genome = "hg19",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfbinding"
)
metadata.tfBinding.mm10 <- data.frame(
    Title = "TF Binding Info Mm10",
    Description = "Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE",
    BiocVersion = "3.16",
    Genome = "mm10",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfbinding"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.tfBinding.hg38 <- data.frame(
    Call = "tfBinding()",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding Hg38",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "600 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.hg19 <- data.frame(
    Call = "tfBinding()",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding Hg19",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "600 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.mm10 <- data.frame(
    Call = "tfBinding()",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding Mm10",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "600 MB",
    Version = "2022-09-20"
)

