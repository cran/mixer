\encoding{latin1}
\name{graph.affiliation}
\alias{graph.affiliation}
\title{Simulation of an Affiliation Graph}

\description{Simulation an affiliation Graph with a given number of
cluster, specific proportions and within between probabilies of connection
}
\usage{
graph.affiliation(n=100,alphaVect=c(1/2,1/2),lambda=0.7,epsilon=0.05) 
}

\arguments{
  \item{n}{Number of nodes of the simulated Graph}
  \item{alphaVect}{Vector of proportions of the clusters}
  \item{lambda}{Within-cluster probability of edge}
  \item{epsilon}{Between-clusters probability of edge}
 }
\details{\code{graph.affiliation} simulate a simple Erdos-Renyi
  Mixture of Graph  model, using the same within-cluster probability
  for all cluster and a unique between-cluster edge probability.
}
\value{\code{graph.affiliation} 
returns a list of 2 objects
\item{x}{A adjacency matrix of size n by n}
\item{cluster}{A vector of integers indicating the cluster to which each node is allocated}
}

\references{}

\author{Christophe Ambroise }
\examples{
graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
str(g)
}
\keyword{graphs}
