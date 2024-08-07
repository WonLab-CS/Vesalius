% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\docType{class}
\name{vesalius_assay-class}
\alias{vesalius_assay-class}
\title{The vesalius_assay class}
\description{
The vesalius_assay class is the functional unit of vesalius. Each assay is
stored within this class and it contains all the required information to
run analysis on your assay of choice. In this object, you can find spatial
tiles, image embeddings, spatial territories, differentially expressed genes
(DEG), count matrices (raw and normalised), microscopy images (if present)
and a functional log that lets you see what had been run on this object.
}
\section{Slots}{

\describe{
\item{\code{assay}}{character assay name}

\item{\code{tiles}}{data.frame containing spatial coordinates and pixels tiles once
they have been computed}

\item{\code{embeddings}}{list containing latent space embeddings in the form of
data.frames.}

\item{\code{active}}{matrix containing active embedding data}

\item{\code{territories}}{data.frame containing spatial color segments, spatial
territories, or layers.}

\item{\code{DEG}}{list of data.frame for each differentially gene expression trial}

\item{\code{counts}}{list that containing count matrices. Raw and normalised will
be stored here and named by the normalisation method used.}

\item{\code{image}}{list containing associated microscopy images (NOT implemented)}

\item{\code{log}}{list containing analysis history of the object.}
}}

