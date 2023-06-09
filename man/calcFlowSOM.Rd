% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insideFlow.R
\docType{methods}
\name{calcFlowSOM}
\alias{calcFlowSOM}
\alias{calcFlowSOM,insideFlow-method}
\title{Perform FlowSOM Clustering}
\usage{
calcFlowSOM(object, applyOnMap, trainedMap, fsXdim, fsYdim, numClust)

\S4method{calcFlowSOM}{insideFlow}(object, applyOnMap, trainedMap, fsXdim, fsYdim, numClust)
}
\arguments{
\item{object}{An insideFlow object}

\item{applyOnMap}{Boolean to apply to trained map or not}

\item{trainedMap}{Object which contains trained map}

\item{fsXdim}{Dimension 1 for grid size used for clustering}

\item{fsYdim}{Dimension 2 for grid size used for clustering}

\item{numClust}{Number of clusters to return from FlowSOM}
}
\value{
Returns object with FlowSOM clustering results
}
\description{
Perform FlowSOM Clustering
}
\examples{
\dontrun{
object <- calcFlowSOM(object, applyOnMap=FALSE, fsXdim=10, fsYdim=10, numClust=8)}
}
