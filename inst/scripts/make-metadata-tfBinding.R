
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.tfBinding.hg38 <- data.frame(
    Title = "TF Binding Info Hg38 (ChIP-Atlas and ENCODE)",
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
    RDataPath = "scMultiome/tfBinding_hg38"
)

metadata.tfBinding.hg38.cistrome <- data.frame(
    Title = "TF Binding Info Hg38 (CistromeDB and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from CistromeDB and ENCODE",
    BiocVersion = "3.16",
    Genome = "hg38",
    SourceType = "BED", # refers to raw data
    SourceUrl = "http://cistrome.org/db/#/bdown",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg38_cistrome"
)
metadata.tfBinding.hg19 <- data.frame(
    Title = "TF Binding Info Hg19 (ChIP-Atlas and ENCODE)",
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
    RDataPath = "scMultiome/tfBinding_hg19"
)
metadata.tfBinding.hg19.cistrome <- data.frame(
    Title = "TF Binding Info Hg19 (CistromeDB and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from CistromeDB and ENCODE",
    BiocVersion = "3.16",
    Genome = "hg19",
    SourceType = "BED", # refers to raw data
    SourceUrl = "http://cistrome.org/db/#/bdown",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg19_cistrome"
)
metadata.tfBinding.mm10 <- data.frame(
    Title = "TF Binding Info Mm10 (ChIP-Atlas and ENCODE)",
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
    RDataPath = "scMultiome/tfBinding_mm10"
)
metadata.tfBinding.mm10.cistrome <- data.frame(
    Title = "TF Binding Info Mm10 (CistromeDB and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from CistromeDB and ENCODE",
    BiocVersion = "3.16",
    Genome = "mm10",
    SourceType = "BED", # refers to raw data
    SourceUrl = "http://cistrome.org/db/#/bdown",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_mm10_cistrome"
)
# dataset manifest metadata
# see ?listDatasets for details
manifest.tfBinding.hg38 <- data.frame(
    Call = "tfBinding(\"hg38\")",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding Hg38 ChIPAtlas+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "280 MB",
    Version = "2022-09-20"
)

manifest.tfBinding.hg38.cistrome <- data.frame(
    Call = "tfBinding(\"hg38_cistrome\")",
    Author = "CistromeDB, ENCODE",
    Title = "TF Binding Hg38 CistromeDB+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "139 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.hg19 <- data.frame(
    Call = "tfBinding(\"hg19\")",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding Hg19 ChIPAtlas+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "275 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.hg19.cistromeDB <- data.frame(
    Call = "tfBinding(\"hg19_cistrome\")",
    Author = "CistromeDB, ENCODE",
    Title = "TF Binding Hg19 CistromeDB+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "138 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.mm10 <- data.frame(
    Call = "tfBinding(\"mm10\")",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding Mm10 ChIPAtlas+ENCODE",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "160 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.mm10.cistrome <- data.frame(
    Call = "tfBinding(\"mm10_cistrome\")",
    Author = "CistromeDB, ENCODE",
    Title = "TF Binding Mm10 CistromeDB+ENCODE",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "41 MB",
    Version = "2022-09-20"
)
