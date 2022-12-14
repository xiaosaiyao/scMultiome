% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scMultiome-package.R
\docType{package}
\name{scMultiome-package}
\alias{scMultiome}
\alias{scMultiome-package}
\title{scMultiome: Collection of Public Single-Cell Multiome (scATAC + scRNAseq) Datasets}
\description{
Single cell multiome data, containing chromatin accessibility (scATAC-seq) and gene expression (scRNA-seq) information analyzed with the ArchR package and presented as MultiAssayExperiment objects.
}
\details{
Single cell multiome data sets, paired and unpaired, were analyzed
with the \code{ArchR} package in order to obtain epiregulons.
\code{ArchR} projects were converted to \code{MultiAssayExperiment} objects.

The creation of all datasets is described in detail in respective help files.
Run \code{listDatasets()} to view a list of available data sets or see \code{Datasets} below.
See \code{?<DATASET_NAME>} for details on particular data sets, e.g. \code{?prostateENZ}.
}
\section{Datasets}{

\itemize{
\item \strong{colonHealthy}: Single-cell analysis of samples from healthy human colon
\item \strong{hematopoiesis}: scATAC-seq and unpaired scRNA-seq of hematopoetic cells
\item \strong{prostateENZ}: LNCaP Cells Treated with Enzalutamide
\item \strong{reprogramSeq}: Reprogram-seq of LNCaP cells
\item \strong{tfBinding_hg19}: TF Binding Info Hg19 (ChIP-Atlas and ENCODE)
\item \strong{tfBinding_hg19_cistrome}: TF Binding Info Hg19 (CistromeDB and ENCODE)
\item \strong{tfBinding_hg38}: TF Binding Info Hg38 (ChIP-Atlas and ENCODE)
\item \strong{tfBinding_hg38_cistrome}: TF Binding Info Hg38 (CistromeDB and ENCODE)
\item \strong{tfBinding_mm10}: TF Binding Info Mm10 (ChIP-Atlas and ENCODE)
\item \strong{tfBinding_mm10_cistrome}: TF Binding Info Mm10 (CistromeDB and ENCODE)
}
}

\author{
\strong{Maintainer}: Xiaosai Yao \email{yao.xiaosai@gene.com}

Authors:
\itemize{
  \item Aleksander Chlebowski \email{aleksander.chlebowski@contractors.roche.com}
  \item Aaron Lun \email{lun.aaron@gene.com}
  \item Tomasz Wlodarczyk \email{tomasz.wlodarczyk@contractors}
  \item Natalie Fox \email{natalie.fox@roche.com}
}

Other contributors:
\itemize{
  \item Shiqi Xie \email{xie.shiqi@gene.com} [contributor]
}

}
\keyword{internal}
