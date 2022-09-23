
# ExperimentHub metadata
metadata.prostateENZ <- data.frame(
    Title = "LNCaP Cells Treated with Enzalutamide",
    Description = "Multiome data (scATAC and scRNAseq) for enzalutamide-treated LNCaP cells. Measured (unpaired) with Illumina NextSeq 500 (GPL18573). Results stored in GEO, acc.nos: GSE168667 and GSE168668.",
    BiocVersion = "3.16",
    Genome = "hg38",
    SourceType = "tar.gz",
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168667&format=file, https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168668&format=file", # refers to raw data; optional
    SourceVersion = "Sep 08 2021",
    Species = "Homo sapiens",
    TaxonomyId = "9606",
    Coordinate_1_based = TRUE, # can stay NA, defaults to TRUE
    DataProvider = "Tampere University", # refers to raw data
    Maintainer = desc::desc_get_maintainer(),
    RDataClass = "MultiAssayExperiment",
    DispatchClass = "FilePath",
    RDataPath = "scMultiome/prostateENZ"
)

# dataset manifest metadata
manifest.prostateENZ <- data.frame(
    Title = "LNCaP Cells",
    Species = "Homo sapiens",
    Type = "cell culture",
    Multiome = "unpaired",
    DiskSize = "2.9 GB",
    MemorySize = "0.6 - 8.5 GB",
    Accessor = "prostateENZ",
    Version = "2022-09-06"
)
