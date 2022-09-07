

# scMultiome

Collection of public single cell multiome data sets.

### Data

The data sets are `ArchR` projects converted to `MultiAssayExperiment` objects and saved in hdf5 files.
Experiments (parts of the MAE) are `SingleCellExperiment` objects diassembled into parts, so you can choose which experiments to load. Assays are saved as sparse arrays to save storage.

Loading the data yields SCEs, in which assays are of class `DelayedMatrix` to save memory.

_NOTE: These data sets can be quite large. See `listDatasets()` to avoid surprises._


### Installation:
```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("scMultiome")
```


### Package structure:

```
.
├── README.md                       this file
├── DESCRIPTION                     package metadata
├── NAMESPACE                       namespace information
├── developer_notes.Rmd             notes for developers
├── R/                              functions
├── man/                            help files
├── vignettes/                      vignettes
├── inst/
│   ├── extdata/                    external data, including data set metadata
│   └── scripts/                    scripts, including ones to create data set metadata
└── scMultiome.Rproj                RStudio project file           
```
    
### Contributing

To add your public data sets to the package, please contact the package mainainer.
