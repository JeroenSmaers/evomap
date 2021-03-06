% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruneSample.R
\name{pruneSample}
\alias{pruneSample}
\title{Pruning data and a tree to a particular subset}
\usage{
pruneSample(data, tree, targetGroup)
}
\arguments{
\item{data}{a vector of tip values for species; should be in the same order as tiplabels in the tree}

\item{tree}{an object of class "phylo".}

\item{targetGroup}{vector listing the tip labels to which the data should be pruned}
}
\value{
pruned comparative data set
}
\description{
Prunes a comparative dataset down to a particular clade
}
