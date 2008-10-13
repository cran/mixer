\name{spmgraph-class}
\docType{class}
\alias{spmgraph-class}

\title{Class "spmgraph"}
\description{

A class capable of reading files in "spm" format into graph-like objects

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("spmgraph", initialdata,...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"} containing the filename the data was read from.}
    \item{\code{edgesrcs}:}{Object of class \code{"numeric"} vector of indices representing edge sources.}
    \item{\code{edgedests}:}{Object of class \code{"numeric"} vector of indices representing edge destinations.}
    \item{\code{nodenames}:}{Object of class \code{"character"} vector of "strings" associating indices to symbolic node names.}
    \item{\code{nodeindexes}:}{Object of class \code{"numeric"}  }
  }
}
\section{Methods}{
\describe{
	\item{getEdges}{\code{signature(withlabels=FALSE)}: returns the edges of the graph as a two-line matrix.}
	\item{getAdjacencyMatrix}{\code{signature(withlabels=FALSE)}: returns the adjacency matrix of the graph.}
	}
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Mark Hoebeke <Mark.Hoebeke@jouy.inra.fr> }
\examples{
showClass("spmgraph")


macaque<-new("spmgraph",file.path(.path.package("mixer"),"spm","macaque_matrix.spm"))
str(macaque)
x<-getAdjacencyMatrix(macaque)
}
\keyword{classes}
