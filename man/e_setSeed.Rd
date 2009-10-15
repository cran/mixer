\name{setSeed}
\alias{setSeed}
\title{Set internal seed}
\description{
This utility sets the internal random seed used by mixer. 
}
\usage{
setSeed( seed = 1 )
}
\arguments{
  \item{seed}{sets the seed of the internal random generator in the mixer
    C/C++ libraries (integer value).}
}
\details{
  
  Sets the seed of the random generator (\code{srand}) in the C standard
  library. This random generator is used inside the mixer initialization
  stage. This function is useful to generate exactly the same initial
  conditions before two \code{mixer} runs. 
}



\author{G. Grasseau }
\examples{

graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
setSeed(777)
mixer(g$x,qmin=2,qmax=6)->xout
\dontrun{plot(xout)}

# Produce strictly the same result
setSeed(777)
mixer(g$x,qmin=2,qmax=6)->xout
\dontrun{plot(xout)}
}
\keyword{graphs}
\keyword{cluster}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
