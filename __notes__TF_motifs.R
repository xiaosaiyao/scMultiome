
# how to save TF binding information into a h5 file compatible with this package

devtools::load_all()

### 1. every peak collection is a GRangesList object
###    create a list of collections

# create vector of available genomes
genomes <- c("hg38", "hg19", "mm10")

# add names
genomes <- stats::setNames(genomes, genomes)

# retrieve TF motif info data for the genomes
motifInfo <- lapply(genomes, getTFMotifInfo)


### 2. save that list into a h5 file

# specify file to save to
fileName <- "inst/extdata/tfBinding_chipatlas_encode.h5"

# create file
rhdf5::h5createFile(fileName)
for (g in names(motifInfo)) {
    cat(g, "\n")
    # create group for genome
    rhdf5::h5createGroup(file = fileName, group = g)
    for (i in names(motifInfo[[g]])) {
        cat(i, ", ")
        # isolate GRanges object
        gr <- motifInfo[[g]][[i]]
        # convet to data frame and save to file (as compound type)
        rhdf5::h5write(storeGR(gr), file = fileName, name = sprintf("%s/%s", g, i))
    }
}


### 3. ready
