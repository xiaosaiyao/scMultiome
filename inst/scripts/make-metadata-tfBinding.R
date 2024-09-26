
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
# merged ChIP-seq
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


# sample specific ChIP-Atlas
metadata.tfBinding.hg38.atlas.sample <- data.frame(
    Title = "TF Binding Info hg38 by sample (ChIP-Atlas)",
    Description = "Transcription factor ChIP-seq data from ChIP-Atlas broken down by sample",
    BiocVersion = "3.19",
    Genome = "hg38",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "3.0", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg38_atlas.sample.rds"
)

metadata.tfBinding.hg19.atlas.sample <- data.frame(
    Title = "TF Binding Info hg19 by sample (ChIP-Atlas)",
    Description = "Transcription factor ChIP-seq data from ChIP-Atlas broken down by sample",
    BiocVersion = "3.19",
    Genome = "hg19",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "3.0", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg19_atlas.sample.rds"
)

metadata.tfBinding.mm10.atlas.sample <- data.frame(
    Title = "TF Binding Info mm10 by sample (ChIP-Atlas)",
    Description = "Transcription factor ChIP-seq data from ChIP-Atlas broken down by sample",
    BiocVersion = "3.19",
    Genome = "mm10",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "3.0", # no commas!
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_mm10_atlas.sample.rds"
)

# sample  specific ENCODE
metadata.tfBinding.hg38.encode.sample <- data.frame(
    Title = "TF Binding Info hg38 by sample (ENCODE)",
    Description = "Transcription factor ChIP-seq data from ENCODE broken down by sample",
    BiocVersion = "3.19",
    Genome = "hg38",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://www.encodeproject.org/",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg38_encode.sample.rds"
)

metadata.tfBinding.hg19.encode.sample <- data.frame(
    Title = "TF Binding Info hg19 by sample (ENCODE)",
    Description = "Transcription factor ChIP-seq data from ENCODE broken down by sample",
    BiocVersion = "3.19",
    Genome = "hg19",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://www.encodeproject.org/",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg19_encode.sample.rds"
)

metadata.tfBinding.mm10.encode.sample <- data.frame(
    Title = "TF Binding Info mm10 by sample (ENCODE)",
    Description = "Transcription factor ChIP-seq data from ENCODE broken down by sample",
    BiocVersion = "3.19",
    Genome = "mm10",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://www.encodeproject.org/",
    SourceVersion = "2021-09-21", # no commas!
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_mm10_encode.sample.rds"
)
# tissue specific ChIP-Atlas
metadata.tfBinding.hg38.atlas.tissue <- data.frame(
    Title = "TF Binding Info hg38 by tissue (ChIP-Atlas)",
    Description = "Transcription factor ChIP-seq data from ChIP-Atlas broken down by tissue",
    BiocVersion = "3.19",
    Genome = "hg38",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "3.0", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg38_atlas.tissue.rds"
)

metadata.tfBinding.hg19.atlas.tissue <- data.frame(
    Title = "TF Binding Info hg19 by tissue (ChIP-Atlas)",
    Description = "Transcription factor ChIP-seq data from ChIP-Atlas broken down by tissue",
    BiocVersion = "3.19",
    Genome = "hg19",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "3.0", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_hg19_atlas.tissue.rds"
)

metadata.tfBinding.mm10.atlas.tissue <- data.frame(
    Title = "TF Binding Info mm10 by tissue (ChIP-Atlas)",
    Description = "Transcription factor ChIP-seq data from ChIP-Atlas broken down by tissue",
    BiocVersion = "3.19",
    Genome = "mm10",
    SourceType = "BED", # refers to raw data
    SourceUrl = "https://github.com/inutano/chip-atlas/",
    SourceVersion = "3.0", # no commas!
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "List",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/tfBinding_mm10_atlas.tissue.rds"
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


# sample specific ChIP-Atlas

manifest.tfBinding.hg38.atlas.sample <- data.frame(
    Call = "tfBinding(\"hg38\", \"atlas.sample\")",
    Author = "ChipAtlas",
    Title = "TF Binding hg38 ChIPAtlas by sample",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "660.2 MB",
    Version = "2024-09-26"
)

manifest.tfBinding.hg19.atlas.sample <- data.frame(
    Call = "tfBinding(\"hg19\", \"atlas.sample\")",
    Author = "ChipAtlas",
    Title = "TF Binding hg19 ChIPAtlas by sample",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "637.4 MB",
    Version = "2024-09-26"
)

manifest.tfBinding.mm10.atlas.sample <- data.frame(
    Call = "tfBinding(\"mm10\", \"atlas.sample\")",
    Author = "ChipAtlas",
    Title = "TF Binding mm10 ChIPAtlas by sample",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "311.4 MB",
    Version = "2024-09-26"
)


# sample specific ENCODE
manifest.tfBinding.hg38.encode.sample <- data.frame(
    Call = "tfBinding(\"hg38\", \"encode.tissue\")",
    Author = "ENCODE",
    Title = "TF Binding hg38 ENCODE by sample",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "167.2 MB",
    Version = "2024-09-26"
)

manifest.tfBinding.hg19.encode.sample <- data.frame(
    Call = "tfBinding(\"hg19\", \"encode.sample\")",
    Author = "ENCODE",
    Title = "TF Binding hg19 ENCODE by sample",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "168 MB",
    Version = "2024-09-26"
)

manifest.tfBinding.mm10.encode.sample <- data.frame(
    Call = "tfBinding(\"mm10\", \"encode.sample\")",
    Author = "ENCODE",
    Title = "TF Binding mm10 ENCODE by sample",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "14.7 MB",
    Version = "2024-09-26"
)

# tissue specific ChIP-Atlas

manifest.tfBinding.hg38.atlas.tissue <- data.frame(
    Call = "tfBinding(\"hg38\", \"atlas.tissue\")",
    Author = "ChipAtlas",
    Title = "TF Binding hg38 ChIPAtlas by tissue",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "466 MB",
    Version = "2024-09-26"
)

manifest.tfBinding.hg19.atlas.tissue <- data.frame(
    Call = "tfBinding(\"hg19\", \"atlas.tissue\")",
    Author = "ChipAtlas",
    Title = "TF Binding hg19 ChIPAtlas by tissue",
    Species = "Homo sapiens",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "448.8 MB",
    Version = "2024-09-26"
)

manifest.tfBinding.mm10.atlas.tissue <- data.frame(
    Call = "tfBinding(\"mm10\", \"atlas.tissue\")",
    Author = "ChipAtlas",
    Title = "TF Binding mm10 ChIPAtlas by tissue",
    Species = "Mus musculus",
    Lineage = "All",
    CellNumber = "Bulk",
    Multiome = "n/a",
    DiskSize = "231.6 MB",
    Version = "2024-09-26"
)
