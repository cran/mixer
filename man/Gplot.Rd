\name{Gplot}
\alias{Gplot}

\title{Graph Plot}

\description{
  
  This graph-drawing function is based on low-level routines
  which have been extracted from the \code{sna} package
  (mainly \code{gplot}) and reorganized to exploit more
  conveniently its possibilies in the \code{mixer} context.

  \code{mixer} graphs have real-valued undirected edges, have 
  coloured vertices, and don't have self-loops.
  Graphical options have been chosen to try to
  be more illustrative for gene regulation networks .
  
}

\usage{
Gplot( X, cl = NULL, ...     )
}
\arguments{
  \item{X}{
    Adjacency matrix or edge list ( matrix(2,N) ).
  }
  \item{cl}{
    node classification vector.
    
    Default \code{NULL}.
  }
  \item{...}{
    Additional optional arguments are available, see \emph{details} section.
  }
}

\value{
  A matrix of node coordinates is returned, along with a plotted graph.  
}


\details{
    This function basically adapts the given precision matrix (and
  optional node classification vector) and creates a graph which is
  plotted using low-level functions extracted from from the \code{sna}
  package function (\code{gplot}). 

    Given the frame of our work, the graph is undirected,
  real-valued, with no self-loops. Edges are coloured by edge weight
  sign (red for negative partial correlation, blue for positive).
  Edge widths are proportional to absolute edge weight, the maximal
  width being set to 5. Additionally, weak edges are dotted rather than solid.
  The function returns the coordinates of the nodes in the plot,
  which can be later passed back to Gplot, effectively skipping the
  coordinate calculation phase.

    Additional arguments are :
  \describe{
    \item{\code{cols}}{
       Colour specification for each node cluster. Default \code{NULL}.
    }
    \item{\code{coord}}{
       \code{matrix} containing the 2D coordinates for the nodes (one node per row). Default \code{NULL}.
    }
    \item{\code{my.labels}}{
       Character \code{vector} : Labels for the nodes. Default \code{NULL}.
    }
    \item{\code{display.labels}}{
       Boolean : display the labels? Default \code{FALSE}.
    }
    \item{\code{degree.threshold}}{
       Scalar value : degree threshold under which nodes are considered as dustbin nodes. Default \code{.Machine$double.xmin}. %$
    }
    \item{\code{max.edges}}{
       Scalar value. Threshold on the number of edges in the graph, over which the graph will NOT BE PLOTTED AT ALL. Default \code{10000}.
    }
    \item{\code{main}, \code{sub}}{
       Graphical options corresponding to generic function \code{plot} arguments. Default \code{"Gplot"} and \code{NULL} respectively.
    }
  }
}


\seealso{
   \code{\link{mixer}}, \code{sna}, \code{gplot}, \code{\link{plot}}
}



\examples{

  # Generate a graph
  g <- graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)

  # Display adjacency matrix
  \dontrun{Gplot(g$x, g$cluster)}

  # Read a new graph
  data(blog)
  \dontrun{Gplot( as.matrix( blog$link ), blog$politicalParty)}
}



\author{ A. Smith and G. Grasseau}


\keyword{htest} 
\keyword{cluster} 