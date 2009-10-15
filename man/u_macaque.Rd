\name{macaque}
\docType{data}
\alias{macaque}
\title{Connection of macaque brain cortical regions}
\description{
The dataset of consists in 47 brain cortical regions connected
by 505 inter-regional pathways in the Macaque Cortex. 
}
\usage{
data(macaque)
}
\format{A data frame desribing the adjacency matrix of the connection
of the 47 brain cortical regions of the Macaque Cortex
}

\details{
As brain function is based on
inter-regional connections, studying the way cortical regions interact may offer new perspectives
in the comprehension of information flows within the brain. It appears that particular brain
regions may play different roles: some regions can be at the ''center'' of a particular part of
the network, meaning that a lot of information will pass through them, whereas other parts
of the network may be more ''peripherica''. Consequently, identifying central zones would be
important, as their lesion may compromise the integrity of the whole network.
From a topological view, finding those ''hubs'' as focused much attention, with a definition
based on degree only. However, there exists many ways for a node to be a hub, and degree is
one criteria. As there is no definition of what a hub is, there are many different hubs (provincial
and central). This is why [1] developed a multi-criteria strategy to find nodes that can be called
''hubs''. From a methodological point of view, this approach seems to be limited as the resuting
hubs will be criteria-dependent.
The gain of Mixer is that the model can be used to find those hubs. Indeed, using the
underlying missing data framework, MixNet will find nodes that connect heavily to other nodes
in the network, and that share this connectivity pattern (a class of hubs for instance).
}

\references{
[1] Sporns O., Honey C., and Kotter R. Identification and classification of hubs in brain networks.
PLOS one, 10:1-14, 2007.
}
\examples{
data(macaque)
mixer(macaque,qmin=8)->xout
\dontrun{plot(xout)}
}
\keyword{datasets}
