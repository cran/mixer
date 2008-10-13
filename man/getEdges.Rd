\name{getEdges}
\docType{genericFunction}
\alias{getEdges}
\alias{getEdges,ANY-method}
\title{Return Graph Edges}
\description{
  This generic function returns the edges of a graph-like structure in a
  two-line matrix. Each column represents an edge, its first line being the
  source node and its second line being the destination node.
}
\usage{
getEdges(g, \dots)
}
\arguments{
  \item{g}{An object of a suitable class.  }
  \item{\dots}{Some methods for this generic function may take
    additional, optional arguments.  At present none do.}
}
\value{
  Returns a two line matrix with as much columns as there are edges in the graph..
}
%\author{}
%\note{}
%\examples{
%}
\keyword{graphs}

