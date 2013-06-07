pkgname <- "mixer"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mixer')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("a_mixer")
### * a_mixer

flush(stderr()); flush(stdout())

### Name: mixer
### Title: MIXtures of Erdos-Renyi random graphs
### Aliases: mixer
### Keywords: graphs cluster

### ** Examples


graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
mixer(g$x,qmin=2,qmax=6)->xout
## Not run: plot(xout)

graph.affiliation(n=50,c(1/3,1/3,1/3),0.8,0.2)->g
mixer(g$x,qmin=2,qmax=5, method="bayesian")->xout
## Not run: plot(xout)

data(blog)
mixer(x=blog$links,qmin=2,qmax=12)->xout
## Not run: plot(xout)



cleanEx()
nameEx("b_graph.affiliation")
### * b_graph.affiliation

flush(stderr()); flush(stdout())

### Name: graph.affiliation
### Title: Simulation of an Affiliation Graph
### Aliases: graph.affiliation
### Keywords: graphs

### ** Examples

graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
str(g)



cleanEx()
nameEx("c_plot.mixer")
### * c_plot.mixer

flush(stderr()); flush(stdout())

### Name: plot.mixer
### Title: Plot of mixer object
### Aliases: plot.mixer
### Keywords: graphs

### ** Examples

#
#  Simple example : display the 4 frames for the best class number estimation
#
g <- graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)
xout <- mixer(g$x,qmin=2,qmax=6)
## Not run: plot(xout)

#
#  Display the same for 4 classes with no filtering
#
## Not run:  plot(xout, q=4, quantile.val=0) 

#
#  Display a pie chart for 4 classes
#
data(blog)
xout <- mixer(x=blog$links,qmin=2,qmax=12)
#  Unconnected nodes have been removed by mixer.
#  xout$map contains the mapping from connected nodes to the whole set 
ext.classes <-  blog$politicalParty
## Not run:  plot( xout, frame=4, classes=ext.classes )



cleanEx()
nameEx("d_getModel.mixer")
### * d_getModel.mixer

flush(stderr()); flush(stdout())

### Name: getModel.mixer
### Title: Get the model parameters
### Aliases: getModel.mixer
### Keywords: graphs cluster

### ** Examples


graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2) -> g
mixer(g$x,qmin=2,qmax=6) -> xout
m <- getModel( xout )




cleanEx()
nameEx("e_setSeed")
### * e_setSeed

flush(stderr()); flush(stdout())

### Name: setSeed
### Title: Set internal seed
### Aliases: setSeed
### Keywords: graphs cluster

### ** Examples


graph.affiliation(n=100,c(1/3,1/3,1/3),0.8,0.2)->g
setSeed(777)
mixer(g$x,qmin=2,qmax=6)->xout
## Not run: plot(xout)

# Produce strictly the same result
setSeed(777)
mixer(g$x,qmin=2,qmax=6)->xout
## Not run: plot(xout)



cleanEx()
nameEx("t_blog")
### * t_blog

flush(stderr()); flush(stdout())

### Name: blog
### Title: French Political Blogosphere network
### Aliases: blog
### Keywords: datasets

### ** Examples

data(blog)
mixer(x=blog$links,qmin=2,qmax=12)->xout
## Not run: plot(xout)



cleanEx()
nameEx("u_macaque")
### * u_macaque

flush(stderr()); flush(stdout())

### Name: macaque
### Title: Connection of macaque brain cortical regions
### Aliases: macaque
### Keywords: datasets

### ** Examples

data(macaque)
mixer(macaque,qmin=8)->xout
## Not run: plot(xout)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
