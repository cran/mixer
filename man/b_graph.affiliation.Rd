\name{graph.affiliation}
\alias{graph.affiliation}
\title{Simulation of an Affiliation Graph}

\description{Simulate an affiliation network with a given number of
  clusters, specific class proportions and within/between connection
  probabilities.
}
\usage{
graph.affiliation( n=100, alphaVect=c(1/2,1/2), lambda=0.7,
                   epsilon=0.05, directed=FALSE ) 
}

\arguments{
  \item{n}{number of nodes of the simulated network.}
  \item{alphaVect}{vector of cluster proportions.}
  \item{lambda}{within-cluster probability of edge.}
  \item{epsilon}{between-clusters probability of edge.}
  \item{directed}{\code{TRUE} for directed graphs.}
 }
\details{\code{graph.affiliation} simulates a simple Erdös-Rényi
  Mixture of Graph  model, using the same within-cluster edge probability
  for all clusters and a unique between-cluster edge probability.
}
\value{\code{graph.affiliation} 
returns a list of 2 objects:
\item{x}{an adjacency matrix of size n by n,}
\item{cluster}{a vector of integers indicating the cluster to which each node is allocated.}
}

% \references{}

\author{Christophe Ambroise }
\examples{
graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
str(g)
}
\keyword{graphs}
