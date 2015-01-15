\name{plot.mixer}
\alias{plot.mixer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of mixer object }
\description{
\code{plot.mixer} can display five kinds of figure: model selection criterion curve, the adjacency matrix map, the degree distribution histogram, the
connectivity matrix graph and the adjacency matrix graph. By default the four
first plots are displayed.}

\usage{
\method{plot}{mixer}(x, q=NULL, frame=1:4, classes=NULL,
classes.col=NULL, quantile.val=0.1, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a mixer object (output of the  mixer function).}
  \item{q}{the q-class model to display. By default, the \code{q}
    is set to the value which maximizes the criterion (see frame 1).} 
  \item{frame}{a vector of frame numbers to display (5 kinds of plots,
    see details section for more).}
  \item{classes}{an external classification used for frame 4 (pie chart):
    vector as \code{factor} of node elements (the number of external
    class levels corresponds to the number of levels). 
  }
  \item{classes.col}{a vector of user colors used to identify the classes.}
  \item{quantile.val}{ filters the connectivity matrix values
    (Pis) in frame 4. Display the upper part (specified by
    \code{quantile.val}) of the distribution.
  }
  
  \item{...}{further graphical arguments.}
}
\details{

  Frame values:
  \describe{
    \item{\bold{1}}{
      criterion (ICL or ILvb) versus the number of classes
      (see  \code{\link{mixer}}).
    }
    \item{\bold{2}}{
      adjacency matrix reorganized according to the estimated partition for
  a given number of classes \code{q}.
    }
    \item{\bold{3}}{
      degree distribution (histogram) and theoretical degree
      distribution (blue curve) computed from the q-class model
      parameters (\code{alphas}, \code{Pis}). 
    }
    \item{\bold{4}}{
      matrix connectivity between classes (Pis) given a number of
      classes \code{q}. The thickest edges
      identify the highest values of the connectivity probabilities
      and the largest nodes point out the most populated classes.

      Providing external classes (see \bold{classes} argument) each
      node displays a pie chart pointing out the classification relevance.
      
    }
    \item{\bold{5}}{
      graph display of the adjacency matrix. 
    }
  }
}
\author{C. Ambroise, G. Grasseau}

\seealso{ \code{\link{mixer}},  \code{\link{getModel}}}
\examples{
#
#  Simple example : display the 4 frames for the best class number estimation
#
g <- graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)
xout <- mixer(g$x,qmin=2,qmax=6)
\dontrun{plot(xout)}

#
#  Display the same for 4 classes with no filtering
#
\dontrun{ plot(xout, q=4, quantile.val=0) }

#
#  Display a pie chart for 4 classes
#
data(blog)
xout <- mixer(x=blog$links,qmin=2,qmax=12)
#  Unconnected nodes have been removed by mixer.
#  xout$map contains the mapping from connected nodes to the whole set 
ext.classes <-  blog$politicalParty
\dontrun{ plot( xout, frame=4, classes=ext.classes )}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphs}

