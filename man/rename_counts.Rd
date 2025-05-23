% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{rename_counts}
\alias{rename_counts}
\title{renaming counts to remain consistent with vesalius nomenclature}
\usage{
rename_counts(integrated, seed_cells, seed_genes, query_cells, query_genes)
}
\arguments{
\item{integrated}{seurat object containing integrated count data}

\item{seed_cells}{character vector contain seed barcodes}

\item{seed_genes}{character vector contain seed genes}

\item{query_cells}{character vector contain query barcodes}

\item{query_genes}{character vector contain query genes}
}
\value{
list of count matrices with cell and gene names added and 
each matrix is renamed to follow the vesalius nomencalture
}
\description{
renaming counts to remain consistent with vesalius nomenclature
}
