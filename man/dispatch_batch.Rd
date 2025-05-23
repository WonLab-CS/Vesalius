% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_assays.R
\name{dispatch_batch}
\alias{dispatch_batch}
\title{Dispatch cells into batches}
\usage{
dispatch_batch(cost_matrix, batch_size = 5000)
}
\arguments{
\item{cost_matrix}{matrix containing cost for each cell pair}

\item{batch_size}{int size of batch}
}
\value{
Nested list. Each element of the list will contain 
a batched cost matrix and the mapping pairs
}
\description{
Dispatch cells into batches
}
\details{
Create cell batches that will dynamically adapt to the size of
the data set with respect to batch size. Smalled data sets, cells
will be sampled to match the size of the larger data set. This
allows for multiple to multiple matching. All cells will be selected
at least once.
}
