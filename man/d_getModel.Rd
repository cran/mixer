\name{getModel}
\alias{getModel}
\title{Get the parameters of a model}
\description{Generic method to get the parameters of a model. Type ?getModel.mixer for  more details}
\usage{getModel( object, \ldots)}
\arguments{
  \item{object}{an object representing a model}
  \item{...}{ adding arguments (depending on the object type)}
}
\value{
  Return parameters of the model.
}

\examples{

graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2) -> g
mixer(g$x,qmin=2,qmax=6) -> xout
m <- getModel( xout )

m <- getModel( xout, q=5 )


}

\keyword{graphs}
\keyword{cluster}
