% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_mod.R
\name{average_embed}
\alias{average_embed}
\title{average image stack between seed and query}
\usage{
average_embed(seed, query, dimensions)
}
\arguments{
\item{seed}{matrix - seed embedding image stack}

\item{query}{matrix - query embedding image stack}

\item{dimensions}{int vector describing which embeddings
should be selected}
}
\value{
embedding matrix containing average pixel value for both seed
query
}
\description{
average image stack between seed and query
}
\details{
Takes select embedding from seed and query and 
creates and avarage the grey scale pixel values for each spatial
location
}
