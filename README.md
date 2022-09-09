

# scMultiome

Collection of public single cell multiome data sets.

### Data

The data sets are `ArchR` projects converted to `MultiAssayExperiment` objects and saved in hdf5 files.
The nature of the hdf5 format allows the MAEs to be split into individual experiments and store them in one file, so you can choose freely which ones to load. 
Experiments, usually `SingleCellExperiment` objects, are disassembled into parts, which are saved in the hdf5 hierarchy. Assays are saved as sparse arrays to save storage.

Upon loading, selected SCEs are reassembled from parts and returned as experiments within an MAE object, in which assays are of class `DelayedMatrix` to save memory.

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
