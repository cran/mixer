\name{getModel.mixer}
\alias{getModel.mixer}
\title{Get the model parameters}
\description{
  Given a number of classes, get the model parameters (\code{alphas},
  \code{Pis}, \code{Taus} ).
}
\usage{
\method{getModel}{mixer}( object, \ldots)
}
\arguments{
  \item{object}{a mixer object (output of the \code{mixer} function).}
  \item{...}{ adding arguments (depending on the object type)
   \describe{
     \item{q}{selects a q-class model. If \code{NULL}, \code{q} value is
       chosen to maximise the criterion (ICL criterion or ILvb criterion)
       (see \link{plot.mixer} and \link{mixer}).
       Default \code{NULL}.
     }
   }
  }
}


\value{
  Return a list with the following attributes:
    \item{q}{q-class model selected.}
    \item{criterion}{criterion value (ICL criterion or ILvb criterion).}
    \item{alphas}{vector of class proportion.}
    \item{Pis}{connectivity matrix of classes.}
    \item{Taus}{matrix of posterior probabilities (of the hidden
  colour knowing the graph structure).}
}
\author{G. Grasseau }
\examples{

graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2) -> g
mixer(g$x,qmin=2,qmax=6) -> xout
m <- getModel( xout )

}
\keyword{graphs}
\keyword{cluster}
