% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/embeddings.R
\name{embed_umap}
\alias{embed_umap}
\title{embed umap}
\usage{
embed_umap(counts, dimensions, features = NULL, verbose)
}
\arguments{
\item{counts}{Seurat object containing normalised counts}

\item{dimensions}{number of PCs to use for the UMAP projections}

\item{features}{custom vector of features}

\item{verbose}{logical if progress messages should be outputed}
}
\value{
normalised UMAP projection matrix
}
\description{
embed in gray scale using UMAP projections
}
\details{
Note that while you can select any number of dimensions
the number of UMAP dimensions will always be 3.
}
