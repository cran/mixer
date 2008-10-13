\name{getAdjacencyMatrix}
\docType{genericFunction}
\alias{getAdjacencyMatrix}
\alias{getAdjacencyMatrix,ANY-method}
\title{Return a Graph's Adjacency Matrix}
\description{
  This generic function returns the adjacency matrix of a graph. Non-null elements denote
  that an edge exists between the nodes corresponding to the row and column index.
}
\usage{
getAdjacencyMatrix(g, \dots)
}
\arguments{
  \item{g}{An object of a suitable class.  }
  \item{\dots}{Some methods for this generic function may take
    additional, optional arguments.  At present none do.}
}
\value{
  Returns a square-matrix where non-null elements denote the existence of edges in the graph.
}
%\author{}
%\note{}
%\examples{
%}
\keyword{graphs}

