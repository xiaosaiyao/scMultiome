

# scMultiome


## Motivation

Single cell data is gaining sophistication - Cells can be measured in multiple modalities including gene expression, chromatin accessibility, cell surface markers and protein expression. These orthogonal measures of the same or matched cells enable a holistic construction of the cell state. However it has been challenging to share multiomic data, especially in an integrated format that consolidates the multiple layers of measurements. The `MultiAssayExperiment` provides a framework to package the various modalities into a single dataset on a per cellbasis.

The `scMultiome` package is a collection of public single cell multiome data sets preprocessed and packaged into MultiAssayExperiment objects for downstream analysis. It also provides basic functions to save the `MultiAssayExperiment` as `.hdf5` files so that users can load only the desired modalities into memory.


## Data format

Current multiomic datasets consist of gene expression and chromatin accessibility but can be extended to include any other modalities. The datasets are either paired multiomic datasets or unpaired datasets with data integration performed by the `ArchR` [package](https://www.archrproject.com/). The `ArchR` projects were converted to `MultiAssayExperiment` objects. [MultiAssayExperiment](https://www.bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html) objects can be constructed easily from individual matrices preprocessed by users' favorite packages.

The `MultiAssayExperiment` object is saved in hdf5 files. The nature of the hdf5 format allows the MAEs to be split into individual experiments and store them in one file, so you can choose freely which ones to load. Experiments, usually `SingleCellExperiment` objects, are disassembled into parts, which are saved in the hdf5 hierarchy. Assays are saved as sparse arrays to save storage.

Upon loading, selected SCEs are reassembled from parts and returned as experiments within an MAE object, in which assays are of class `DelayedMatrix` to save memory.


_NOTE: These data sets can be quite large. See `listDatasets()` to avoid surprises._


## Installation:
```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("scMultiome")
```


## Package Contents

To list currently available data sets, use `listDatasets()` or see package help with `?scMultiome`.


## Package structure:

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
    

## Contributing

To add your public data sets to the package, review the vignette _Adding Data Sets_ and contact the package maintainer.

