\name{blog}
\docType{data}
\alias{blog}
\title{French Political Blogosphere network}
\description{
French Political Blogosphere network dataset consists of a single day snapshot
of over 1,100 political blogs automatically extracted the 14 october 2006 and manually classified
by the ''Observatoire Presidentielle'' project. This project is the result of a collaboration
between RTGI SAS and Exalead and aims at analyzing the French presidential campaign on the web.

}
\usage{data(blog)}
\format{
 196 political blogs described by their political origin and connections. \cr
 There are 2 items in the list blog
 \describe{
  \item{links}{ Adjacency matrix describing if there exist an
  hyperlink between two blog}
  \item{politicalParty}{ political group of the blogs}
}
}
\details{
In this data set, nodes represent hostnames (a hostname contains a set of pages) and edges
represent hyperlinks between different hostnames. If several links exist between two different
hostnames, we collapse them into a single one. Note that intra domain links can be considered if
hostnames are not identical. Finally, in this experimentation we consider that edges are not
oriented which is not realistic but which does not affect the interpretation of the groups. This network
presents an interesting communities organization due to the existence of several political parties and
commentators. We assume that authors of these blogs tend to link, by
political affinities, 
blogs
with similar political positions. 
Six known communities compose this network : \bold{G}auche (''french
democrat''), \bold{D}ivers \bold{C}entre (Moderate party),
\bold{D}roite (french republican),  \bold{E}cologiste (green),
\bold{L}iberal (supporters of economic-liberalism)
and finally \bold{A}nalysts. 
}
\references{
 \url{http://www.observatoire-presidentielle.fr/}
}
\examples{
data(blog)
mixer(x=blog$links,qmin=2,qmax=12)->xout
\dontrun{plot(xout)}
}
\keyword{datasets}
