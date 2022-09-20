
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.moitfs.hg38 <- data.frame(
    Title = "TF Motif Info Hg38",
    Description = "Combined transcription factor ChIP-seq data from Cistrone and ENCODE",
    BiocVersion = "3.16",
    Genome = "hg19",
    SourceType = "tar.gz", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/",
    SourceVersion = "2022-09-20", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/moitfs"
)
metadata.moitfs.hg19 <- data.frame(
    Title = "TF Motif Info Hg19",
    Description = "Combined transcription factor ChIP-seq data from Cistrone and ENCODE",
    BiocVersion = "3.16",
    Genome = "hg19",
    SourceType = "tar.gz", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/",
    SourceVersion = "2022-09-20", # no commas!
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/moitfs"
)
metadata.moitfs.mm10 <- data.frame(
    Title = "TF Motif Info Mm10",
    Description = "Combined transcription factor ChIP-seq data from Cistrone and ENCODE",
    BiocVersion = "3.16",
    Genome = "mm10",
    SourceType = "tar.gz", # refers to raw data
    SourceUrl = "https://www.ncbi.nlm.nih.gov/",
    SourceVersion = "2022-09-20", # no commas!
    Species = "Mus musculus",
    TaxonomyId = "10090",
    Coordinate_1_based = TRUE,
    DataProvider = "Genentech", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "GRangesList",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/moitfs"
)

# dataset manifest metadata
# see ?listDatasets for details
manifest.moitfs.hg38 <- data.frame(
    Title = "TF Motifs Hg38",
    Species = "Homo sapiens",
    Type = "cell line",
    Multiome = "n/a",
    DiskSize = "600 MB",
    MemorySize = "250 MB",
    Accessor = "motifs",
    Version = "2022-09-20"
)
manifest.moitfs.hg.19 <- data.frame(
    Title = "TF Motifs Hg19",
    Species = "Homo sapiens",
    Type = "cell line",
    Multiome = "n/a",
    DiskSize = "600 MB",
    MemorySize = "230 MB",
    Accessor = "motifs",
    Version = "2022-09-20"
)
manifest.moitfs.mm10 <- data.frame(
    Title = "TF Motifs Mm10",
    Species = "Mus musculus",
    Type = "cell line",
    Multiome = "n/a",
    DiskSize = "600 MB",
    MemorySize = "60 MB",
    Accessor = "motifs",
    Version = "2022-09-20"
)
