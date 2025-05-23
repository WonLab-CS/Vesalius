% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{back_infer}
\alias{back_infer}
\title{Infer normalized counts from shared embedding values}
\usage{
back_infer(counts, embeds)
}
\arguments{
\item{counts}{list of count matrices containing normalized gene expression}

\item{embeds}{shared embedding matrix produced by Seurat.}
}
\value{
list containing inferred count values from embeddings
}
\description{
Infer normalized counts from shared embedding values
}
