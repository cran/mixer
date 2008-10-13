\name{mixer}
\alias{mixer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture of Erdös-Rényi Random Graphs}
\description{
Estimate the parameters and hidden class variable of a Mixture of
Erdös Rényi Random Graphs. The estimation is performed for a 
Bernoulli Mixture 
}
\usage{
mixer(x,qmin=2,qmax=NULL,method="variational",nbiter=10,improve=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Adjacency matrix or a matrix of egdes (each column gives
  the nodes defining an edge) or spm filename (a .spm
  file describes the network as a sparse matrix}
  \item{qmin}{Number of classes estimated (if used alone, it is the
  minimum number of classes which will be tested)}
  \item{qmax}{Maximimum number of classes (use with qmin)}
  \item{method}{strategy used for the  estimation: "classification", 
     "variational", or "bayesian"}
  \item{nbiter}{Number of EM iterations for the E step  (Default: 10)}
  \item{improve}{Flag to choose between improved or basic strategies
  for chosing the model}
}
\details{
 \code{mixer}  is a R wrapper for the c++ software package mixnet
 developped by Vincent Miele (2006). 

Erdös-Rényi Mixture Model for Graph (MixNet), which 
has been proposed by Daudin et. al (2008) with an associated
EM  estimation algorithm, and is not to be confused with Exponential Random Graph Models for 
Network Data (ERGM) which consider distributions ensuing from the
exponential  family to model the edge distribution. The MixNet model
allows  to capture the structure of a network and in particular to detect communities.

There exists a strong connection between Mixnet and block clustering..  
Block clustering searches for homogeneous blocks in a data matrix by
simultaneous  clustering of rows and columns. 

The proposed estimation strategies deals with undirected graphs. They
are of two type: 
\itemize{
\item{variational}{ which refers to the paper of Daudin et. al
(2008). It is the default method. }
\item{classification}{ which implements the method described in Zanghi
et. al (2008). This method is faster than the variational approach and
is able to deal with bigger network but can produce biased estimates.}
 \item{bayesian}{ which implements the method described in Latouche
et. al (2008). It improves over variational and classification
strategies when dealing with small networks (less than 50
nodes). It should not be used with networks having more than 100 nodes. } 

}


}
\value{
 \code{mixer}
returns an object of class mixer, which is basically a list of item,
each item containing the result of the estimation for a given number
of class q. It has  the following  fields
  \item{ICL }{Integrated Classification Likelihood, which is the
  criterion used for model selection}
  \item{alphas}{The vector of proportion, whose length is the number
  of component}
  \item{Pis}{The connectivity matrix}
  \item{Taus}{The matrix of posterior probabilities (of the hidden
  color knowing the graph structure)}
}
\references{
Jean-Jacques Daudin, Franck Picard and Stephane Robin 
June (2008) , 
\emph{A mixture model for random graphs.}
Statistics and Computing, 18, 2, 151--171.

Hugo Zanghi, Christophe Ambroise and Vincent Miele (2008),
\emph{Fast online graph clustering via Erdös-Rényi mixture.}
Pattern Recognition, 41, 3592-3599.


Hugo Zanghi, Franck Picard, Vincent Miele, and Christophe Ambroise (2008),
\emph{Strategies for Online Inference of Network Mixture},
 \url{http://genome.jouy.inra.fr/ssb/preprint/AMPZ07-SSBpreprint.pdf}

Pierre Latouche,  Etienne Birmele, and Christophe Ambroise (2008),
\emph{Bayesian methods for graph clustering},
 \url{http://genome.jouy.inra.fr/ssb/preprint/SSB-RR-17.bayesianMixNet.pdf}

Vincent Miele. MixNet C++ package,
\url{http://stat.genopole.cnrs.fr/sg/software/mixer/}.

}


\author{Christophe Ambroise }
\examples{

graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
mixer(g$x,qmin=2,qmax=6)->xout
\dontrun{plot(xout)}

graph.affiliation(n=50,c(1/3,1/3,1/3),0.8,0.2)->g
mixer(g$x,qmin=2,qmax=5, method="bayesian")->xout
\dontrun{plot(xout)}

\dontrun{plot(xout)}
data(blog)
mixer(x=blog$links,qmin=2,qmax=12)->xout
\dontrun{plot(xout)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}
\keyword{cluster}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
