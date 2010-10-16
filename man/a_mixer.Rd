\name{mixer}
\alias{mixer}
\title{MIXtures of Erdos-Renyi random graphs}
\description{
Estimate the parameters, the hidden class variables, as well as the
number of classes \code{q} of a MIXture of Erdos Renyi random graphs.
The estimation is performed for binary graphs
(edges are assumed to be drawn from Bernoulli distributions). 
}
\usage{
mixer( x, qmin=2, qmax=NULL, method="variational", directed = NULL,
nbiter=10, fpnbiter=5, improve=FALSE, verbose=TRUE)
}
\arguments{
  \item{x}{an adjacency matrix or a matrix of edges (each column gives
  the two node indexes defining an edge) or a spm file name (a \code{.spm}
  file describes the network as a sparse matrix).}
  \item{qmin}{minimum number of classes.}
  \item{qmax}{maximum number of classes (if \code{NULL},
    only \code{q=qmin} is considered).}
  \item{method}{strategy used for the  estimation: 
     "variational", "classification", or "bayesian"}
  \item{directed}{ \code{TRUE}/\code{FALSE} for directed/undirected
    graph. 
    Default is \code{NULL}, i.e. according to the input array \code{x},
    \code{mixer} identifies whether the graph is directed or undirected. }
  \item{nbiter}{maximum number of EM iterations (default: 10).}
  \item{fpnbiter}{maximum number of internal iterations for the E step
    (default: 5).}
  \item{improve}{selects between improved or basic strategies
    (default: \code{FALSE}).}
  \item{verbose}{display warning messages (default: \code{TRUE}).}
}
\details{
  \code{mixer} implements the Erdos-Renyi mixture model for graphs
  (called MixNet) which has been proposed by Daudin et. al (2008)
  with an associated EM  estimation algorithm.
  The MixNet model is well suited to capture the structure of
  a network and in particular to detect communities.

  MixNet must not to be confused with Exponential Random Graph Models for 
  Network Data (ERGM) which considers distributions ensuing from the
  exponential family to model the edge distribution.
  
  There exists a strong connection between Mixnet and block clustering.
  Block clustering searches for homogeneous blocks in a data matrix by
  simultaneous  clustering of rows and columns. 
  
  The \code{mixer} package implements three different estimation strategies
  which were developed to deal with directed and undirected graphs: 
  
   \describe{
    \item{variational}{refers to the paper of Daudin et. al
      (2008). It is the default method.}
    \item{classification}{implements the method described in Zanghi
      et. al (2008). This method is faster than the variational approach and
      is able to deal with bigger networks but can produce biased
    estimates.}
    \item{bayesian}{implements the method described in Latouche
      et. al (2008).}
    
  }

  The implementation of the two first methods consists of an R wrapper of
  the c++ software package \emph{mixnet} developed by Vincent Miele
  (2006).
  
  The mixer routine uses the estimation strategy described in
  \code{method} and computes a model selection criterion for each value
  of \code{q} (the number of classes) between \code{qmin} and
  \code{qmax}. The ICL criterion is used for the \code{variational} and
  \code{classification} methods. It corresponds to an asymptotic
  approximation of the Integrated Classification Likelihood. The other
  criterion, so called  ILvb (Integrated Likelihood variational
  Bayes), is used for the \code{bayesian} method. It is based on a variational
  (non-asymptotic) approximation of the Integrated observed Likelihood.
 
  \code{mixer} is an user-friendly package with a reduced number of functions.
  For R-developers in statistical networks a more complete set, called
  \code{mixer-dev}, is provided (see below).
  

}
\value{
 \code{mixer} returns an object of class mixer. Below the main attributes of this
 class:
 \item{nnodes}{number of connected nodes.}
 \item{map}{mapping from connected nodes to the whole set of nodes.}
 \item{edges}{edge list.}
 \item{qmin, qmax}{number of classes.}
 \item{output}{output list of \code{qmax-qmin+1} items. Each item
   contains the result of the estimation for a given number
   of class q. Details of output field:}
 \item{output[[i]]$criterion}{ICL criterion or ILvb criterion used for
   model selection (see details section for more).}
  \item{output[[i]]$alphas}{vector of proportion, whose length is the number
  of component.}
  \item{output[[i]]$Pis}{class connectivity matrix.}
  \item{output[[i]]$Taus}{matrix of posterior probabilities (of the hidden
  color knowing the graph structure).}
}
\references{
Jean-Jacques Daudin, Franck Picard and Stephane Robin 
June (2008), \emph{A mixture model for random graphs.}
Statistics and Computing, 18, 2, 151--171.

Hugo Zanghi, Christophe Ambroise and Vincent Miele (2008),
\emph{Fast online graph clustering via Erdos-Renyi mixture.}
Pattern Recognition, 41, 3592-3599.


Hugo Zanghi, Franck Picard, Vincent Miele, and Christophe Ambroise (2008),
\emph{Strategies for Online Inference of Network Mixture},\cr
\url{http://arxiv.org/abs/0910.2034v1}

Pierre Latouche,  Etienne Birmele, and Christophe Ambroise (2008),
\emph{Bayesian methods for graph clustering},\cr
\url{http://genome.jouy.inra.fr/ssb/preprint/SSB-RR-17.bayesianMixNet.pdf}

Vincent Miele, MixNet C++ package,\cr
\url{http://stat.genopole.cnrs.fr/sg/software/mixer/}.

\code{mixer-dev} tool: see \url{http://ssbgroup.fr/mixnet/mixer.html}
}


\author{ Christophe Ambroise, Gilles Grasseau, Mark Hoebeke, Pierre
  Latouche, Vincent Miele, Franck Picard }
\examples{

graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
mixer(g$x,qmin=2,qmax=6)->xout
\dontrun{plot(xout)}

graph.affiliation(n=50,c(1/3,1/3,1/3),0.8,0.2)->g
mixer(g$x,qmin=2,qmax=5, method="bayesian")->xout
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
