
# ExperimentHub metadata
# see ?ExperimentHubData::makeExperimentHubMetadata for details
metadata.tfMotifs.human <-
  data.frame(
    Title = "TF motifs human",
    Description = "Transcription factor motifs in the human genome from chromVARmotifs package",
    BiocVersion = "3.19",
    Genome = "",
    SourceType = "RDA", # refers to raw data
    SourceUrl = "https://github.com/GreenleafLab/chromVARmotifs/raw/master/data/human_pwms_v2.rda", # refers to raw data
    SourceVersion = "2", # no commas!
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    TaxonomyId = "9606", # e.g. "9606"
    Coordinate_1_based = NA, # can stay NA, defaults to TRUE
    DataProvider = "Stanford University School of Medicine", # refers to raw data
    Maintainer = desc::desc_get_maintainer(), # refers to package maintainer
    RDataClass = "PWMatrixList", # class that is returned from hub
    DispatchClass = "FilePath", # format saved on disk; FilePath only returns file location
    # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
    RDataPath = "scMultiome/human_pwms_v2.rds"
)

metadata.tfMotifs.mouse <-
    data.frame(
        Title = "TF motifs mouse",
        Description = "Transcription factor motifs in the mouse genome from chromVARmotifs package",
        BiocVersion = "3.19",
        Genome = "",
        SourceType = "RDA", # refers to raw data
        SourceUrl = "https://github.com/GreenleafLab/chromVARmotifs/raw/master/data/human_pwms_v2.rda", # refers to raw data
        SourceVersion = "2", # no commas!
        Species = "Mus musculus", # e.g. "Homo sapiens"
        TaxonomyId = "10090", # e.g. "9606"
        Coordinate_1_based = NA, # can stay NA, defaults to TRUE
        DataProvider = "Stanford University School of Medicine", # refers to raw data
        Maintainer = desc::desc_get_maintainer(), # refers to package maintainer
        RDataClass = "PWMatrixList", # class that is returned from hub
        DispatchClass = "FilePath", # format saved on disk; FilePath only returns file location
        # Location_Prefix = "", # SKIP if data stored in the Bioconductor AWS S3
        RDataPath = "scMultiome/mouse_pwms_v2.rds"
    )


# dataset manifest metadata
# see ?listDatasets for details
manifest.tfMotifs.human <-
  data.frame(
    Call = "tfMotifs('human')",
    Author = "Greenleaf Lab",
    Title = "TF motifs human",
    Species = "Homo sapiens", # e.g. "Homo sapiens"
    Lineage = "n/a",
    CellNumber = "n/a",
    Multiome = "n/a",
    DiskSize = "256 KB",
    Version = "2"
)


manifest.tfMotifs.human <-
    data.frame(
        Call = "tfMotifs('mouse')",
        Author = "Greenleaf Lab",
        Title = "TF motifs human",
        Species = "Mus musculus", # e.g. "Homo sapiens"
        Lineage = "n/a",
        CellNumber = "n/a",
        Multiome = "n/a",
        DiskSize = "256 KB",
        Version = "2"
    )
