
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.tfBinding.hg38 <- data.frame(
    Title = "TF Binding Info hg38 (ChIP-Atlas and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE",
    BiocVersion = "3.17",
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
    RDataPath = "scMultiome/tfBinding_hg38_atlas.rds"
)

metadata.tfBinding.hg38.cistrome <- data.frame(
    Title = "TF Binding Info hg38 (CistromeDB and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from CistromeDB and ENCODE",
    BiocVersion = "3.17",
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
    RDataPath = "scMultiome/tfBinding_hg38_cistrome.rds"
)
metadata.tfBinding.hg19 <- data.frame(
    Title = "TF Binding Info hg19 (ChIP-Atlas and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE",
    BiocVersion = "3.17",
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
    RDataPath = "scMultiome/tfBinding_hg19_atlas.rds"
)
metadata.tfBinding.hg19.cistrome <- data.frame(
    Title = "TF Binding Info hg19 (CistromeDB and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from CistromeDB and ENCODE",
    BiocVersion = "3.17",
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
    RDataPath = "scMultiome/tfBinding_hg19_cistrome.rds"
)
metadata.tfBinding.mm10 <- data.frame(
    Title = "TF Binding Info mm10 (ChIP-Atlas and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE",
    BiocVersion = "3.17",
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
    RDataPath = "scMultiome/tfBinding_mm10_atlas.rds"
)
metadata.tfBinding.mm10.cistrome <- data.frame(
    Title = "TF Binding Info mm10 (CistromeDB and ENCODE)",
    Description = "Combined transcription factor ChIP-seq data from CistromeDB and ENCODE",
    BiocVersion = "3.17",
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
    RDataPath = "scMultiome/tfBinding_mm10_cistrome.rds"
)
# dataset manifest metadata
# see ?listDatasets for details
manifest.tfBinding.hg38 <- data.frame(
    Call = "tfBinding(\"hg38\", \"atlas\")",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding hg38 ChIPAtlas+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "280 MB",
    Version = "2022-09-20"
)

manifest.tfBinding.hg38.cistrome <- data.frame(
    Call = "tfBinding(\"hg38\", \"cistrome\")",
    Author = "CistromeDB, ENCODE",
    Title = "TF Binding hg38 CistromeDB+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "139 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.hg19 <- data.frame(
    Call = "tfBinding(\"hg19\", \"atlas\")",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding hg19 ChIPAtlas+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "275 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.hg19.cistromeDB <- data.frame(
    Call = "tfBinding(\"hg19\", \"cistrome\")",
    Author = "CistromeDB, ENCODE",
    Title = "TF Binding hg19 CistromeDB+ENCODE",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "138 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.mm10 <- data.frame(
    Call = "tfBinding(\"mm10\", \"atlas\")",
    Author = "ChipAtlas, ENCODE",
    Title = "TF Binding mm10 ChIPAtlas+ENCODE",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "160 MB",
    Version = "2022-09-20"
)
manifest.tfBinding.mm10.cistrome <- data.frame(
    Call = "tfBinding(\"mm10\", \"cistrome\")",
    Author = "CistromeDB, ENCODE",
    Title = "TF Binding mm10 CistromeDB+ENCODE",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "41 MB",
    Version = "2022-09-20"
)
