% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/embeddings.R
\name{embed_lsi}
\alias{embed_lsi}
\title{embed lsi}
\usage{
embed_lsi(counts, dimensions, features = NULL, remove_lsi_1, verbose = TRUE)
}
\arguments{
\item{counts}{Seurat object containing normalised counts}

\item{dimensions}{numeric for number of latent space dimensions to use}

\item{features}{custom vector of features}

\item{remove_lsi_1}{logical if first LSI dimenions should be removed}

\item{verbose}{logical if progress messages should be outputed}
}
\value{
normalised LSI embedding matrix
}
\description{
embed in grey scale using latent semantic indexing
}
