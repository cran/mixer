mixer<-function(x,qmin=2,qmax=NULL,method="variational",nbiter=10,improve=FALSE)
  {

  ## How the graph is coded ?
  if (is.character(x)){
    ## spm file
       g <- new("spmgraph",x)
       m <- getEdges(g)
     }
  else if (dim(x)[1]==dim(x)[2]){
    ## Adjacency matrix: convert the adjacency matrix to the required edge list format  
        AdjMat2Edges(x)->m
      }
  else {
    ## Edge list
       m<-x
    }

  nbrNodes<-max(m)  
  
  ## How to find the proper order.....
  readingOrder<-unique(as.numeric(m));
  

  ## prepare the arguments
  undirected<-TRUE
  loop<-FALSE
  kmeansnbclass<- 0   # Accelerate the initialization (used to start the HAC (sould be between NbrNodes and qmax)
  kmeansnbiter<-30     #  
  emeps<-1e-10          # tolerance for em (compare the likelihood)
  fpeps<-1e-4             # tolerance for the fixed point internal loop (on the taus...)
  nokmeans<-TRUE     # No acceleration via kmeans
  fpnbiter<-5               # fixed point nbiter
  
  silent<-TRUE         # no verbose
  initnbv<-0           # size of the initial network for the online version
  improvenbiter<-3     # number of iteration for the improvment phase


  ## ensure the options compatibility
  if  (method=="classification") {
    classif<-TRUE; stochastique<-FALSE; online<-TRUE}
  else    {
    stochastique<-FALSE; classif<-FALSE; online<-FALSE}


  if (undirected==TRUE){
    symetrize<-TRUE}
  else {
    symetrize<-FALSE} 
    
  ## Ensure number of classes coherence
  if (is.null(qmax)){
    qmax<-qmin+0
  }
  
  ## compute the size of the returned array from the .C call
  
  nbrClasses <- qmax - qmin + 1
  span       <- qmin:qmax
  nbrICL     <- nbrClasses;         elts <- c(nbrICL)
  nbrAlphas  <- sum(span);          elts <- c(elts, nbrAlphas)
  nbrPis     <- sum(span*span);     elts <- c(elts, nbrPis)
  nbrTaus    <- nbrNodes*nbrAlphas; elts <- c(elts, nbrTaus)
  nbrValues  <- sum(elts)

  ##Chose the method for the parameter estimation
  if (method=="bayesian"){
    bout <- VariationalBayes(m, qmin, qmax, nbiter, fpnbiter, emeps, fpeps)
  }
  else if (method=="classification"){
    xout <- .C("main_ermgo",
          as.integer(loop),
          as.integer(silent),
          as.integer(initnbv),
          as.integer(improvenbiter),
          as.integer(nbiter),
          as.integer(improve),
          as.integer(classif),
	  as.integer(stochastique),
          as.integer(qmax),
          as.integer(qmin),
          nbrEdges = as.integer(length(m)/2),# size of the given array
          size = as.integer(nbrValues),  # size of the returned array
          lNodes = as.integer(m),       # given array
          res = double(nbrValues))  # returned array
     y <- vector("list", length(span))
   } else {
         xout <- .C("main_ermg",
          as.integer(symetrize),
          as.integer(loop),
          as.integer(undirected),
          as.integer(silent),
          as.integer(improvenbiter),
          as.integer(kmeansnbclass),
          as.integer(kmeansnbiter),
          as.integer(nbiter),
          as.double(emeps),
          as.integer(fpnbiter),
          as.double(fpeps),
          as.integer(improve),
          as.integer(classif),
          as.integer(nokmeans),
          as.integer(qmax),
          as.integer(qmin),
          nbrEdges = as.integer(length(m)/2),# size of the given array
          size = as.integer(nbrValues),  # size of the returned array
          lNodes = as.integer(m),       # given array
          res = double(nbrValues))  # returned array
         
         y <- vector("list", length(span))
       }

  if (method != "bayesian") {
    j <- 1
    cur <- 1
    for (i in span){
      ## format : y[[j]]$name <- dataFormat(x$res[cur:end]);
      ##          cur <- (offset equal to the size of dataFormat)
      y[[j]]$criterion    <- xout$res[cur]; cur <- cur+1
      y[[j]]$alphas <- xout$res[cur:(cur-1+i)]; cur <- cur+i
      y[[j]]$Pis    <- matrix(xout$res[cur:(cur-1+(i*i))], i,i);
                       cur <- cur+(i*i)
      y[[j]]$Taus   <- matrix(xout$res[cur:(cur-1+(i*nbrNodes))], i,
                              nbrNodes,byrow=TRUE); cur <- cur+(i*nbrNodes)
      y[[j]]$Taus[,readingOrder] <- y[[j]]$Taus
      j <- j+1
    }
    result<-list(method=method,edges=m,qmin=qmin,qmax=qmax,output=y)
  } else {
    result<-list(method=method,edges=m,qmin=qmin,qmax=qmax,output=bout)
  }
  class(result)<-"mixer"
  return(result)
}
############################################################
# Plot the icl criterion
############################################################

ploticl<-function(x,q,...)
  {
    if (x$method == "bayesian" ){
      title = "Bayesian criterion versus class number"
      y.lab = "Bayesian criterion"
    } else {
      title = "Integrated Classification Likelihood"
      y.lab = "ICL"
    }
    Q<-unlist(lapply(x$output,ICL<-function(x) length(x$alphas)))
    ICL<-unlist(lapply(x$output,ICL<-function(x) x$criterion))
    plot(Q,ICL,xlab="Number of classes",ylab=y.lab,main=title)
    lines(Q,ICL)
    abline(v=q,col="red",lty=2)
  }

############################################################
# Plot the reorganized adjacency matrix
############################################################

plotam<-function(edges,cluster)
  {
    neworder<-order(cluster)
    max(edges)->n
    m<-t(matrix(order(neworder)[as.numeric(edges)],2))
    plot(1, 1, xlim = c(0, n + 1), ylim = c(n + 1, 0), type = "n", axes= FALSE,xlab="classes",ylab="classes",main="Reorganized Adjacency matrix")
    rect(m[,2]-0.5,m[,1]-0.5,m[,2]+0.5,m[,1]+0.5,col=1)
    rect(m[,1]-0.5,m[,2]-0.5,m[,1]+0.5,m[,2]+0.5,col=1)
    table(cluster)->limits # find the class limits
    cumsum(limits)[1:(length(limits)-1)]+0.5->limits
    abline(v=c(0.5,limits,n+0.5),h=c(0.5,limits,n+0.5),col="red")
  }

############################################################
# Plot the Pis matrix and alphas vector using spectral decomposition
############################################################
plotparam<-function(Pis,alphas,q=NULL){
length(alphas)->q
if (q==1) {D<-list(vector=data.frame(1,1)); a<-b<-1} else {
if (q==2) {a<-b<-1} else {a<-2; b<-3}
D<-colSums(Pis)
L<-diag(rep(1,q)) -  diag(D^(-1/2)) %*% Pis %*% diag(D^(-1/2))
eigen(L)->D
}
plot(D$vector[,a],D$vector[,b],cex=1/min(alphas^(1/2))*alphas^(1/2)*3,axes=FALSE,xlab="",ylab="",main="Specral view of the connection matrix",pch=19,col="red")
points(D$vector[,a],D$vector[,b],cex=1/min(alphas^(1/2))*alphas^(1/2)*3)

text(D$vector[,a],D$vector[,b],label=1:q)  
#gplot((Pis>median(Pis))*Pis,vertex.cex=1/min(alphas^(1/2))*alphas^(1/2)*3,edge.lwd=(Pis>median(Pis))*Pis*1/min(median(Pis)),label=1:length(alphas),label.pos=6)
}


############################################################
# Plot the reorganized adjacency matrix
############################################################


mixture<-function(x,alphas,lambdaq){fx<-0; for (q in 1:length(alphas)) fx<-fx+alphas[q]*dpois(x,lambda=lambdaq[q])}

plotmixture<-function(degrees,Pis,alphas,n){
  colSums(Pis*alphas)*(n-1)->lambdaq
  min(degrees):max(degrees)->x
  mixture(x,alphas,lambdaq)->y
  histo<-hist(degrees,plot=FALSE)
  plot(histo,ylim=c(0,max(histo$density,y)),freq=FALSE,col=7,main="Degree distribution",)
  lines(x,y,lwd=2,col="blue")
  points(x,y)
  }
  

############################################################
# Plot the estimated degree distribution
############################################################
is.mixer<-function(x){if (class(x)=="mixer") TRUE else FALSE}

plot.mixer<-function(x,q=NULL,...){  
x->mixer.res
par(mfrow=c(1,2))
n<-dim(mixer.res$x)[1]
if (!is.mixer(mixer.res)) stop("Not a mixer object")
 
if (is.null(q)) {# find the best number of classes according ICL
                  ICL<-unlist(lapply(mixer.res$output,
                                     ICL<-function(x) x$criterion))
                  which.max(ICL)->i
                  q<-length(mixer.res$output[[i]]$alphas)
                }

apply(mixer.res$output[[q-mixer.res$qmin+1]]$Taus,2,which.max)->cluster
Pis<-mixer.res$output[[q-mixer.res$qmin+1]]$Pis
alphas<-mixer.res$output[[q-mixer.res$qmin+1]]$alphas

ploticl(mixer.res,q)
plotam(mixer.res$edges,cluster)
x11()
Gplot(mixer.res$edges, cl=cluster, main="Graph")
}


 
############################################################
# Simulation of an affiliation graph
############################################################

class.ind<-function (cl)
{ 
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
}


graph.affiliation<-function(n=100,alphaVect=c(1/2,1/2),lambda=0.7,epsilon=0.05) {
      # INPUT  n: number of vertex
      #           alphaVect : vecteur of class proportion
      #           lambda: proba of edge given  same classe
      #           epsilon: proba of edge given two different classes
      # OUTPUT x: adjacency matrix
      #              cluster: class vector
      #           
     
      x<-matrix(0,n,n);
      Q<-length(alphaVect);
      rmultinom(1, size=n, prob = alphaVect)->nq;
      Z<-class.ind(rep(1:Q,nq));
      Z<-Z[sample(1:n,n),];
      for (i in 1:n)
        for (j in i:n)
            {
            # if i and j in same class
            if (which.max(Z[i,])  == which.max(Z[j,])) p<-lambda else  p<-epsilon
            if ((rbinom(1,1,p))&(i != j)) {x[i,j]<-1; x[j,i]<-1}
            }
       return(list(x=x,cluster=apply(Z,1,which.max)) )   
  }


##############################################################
#  Spectral Clustering using normalized Laplacian
##############################################################
spectralkmeans<-function(x,q=2){
  #INPUT:
  #    x is an adjacency matrix
  #OUTPUT:
  #    An object of class "kmeans" which is a list with components:
  n<-dim(x)[1]
  D<-colSums(x)
  L<-diag(rep(1,n)) -  diag(D^(-1/2))%*% x %*% diag(D^(-1/2))
  eigen(L)->D
  kmeans(as.matrix(D$vectors[,max(1,(n-q)): (n-1)]),q)->res         
}

##############################################################
#  Compute the rand index between two partition
##############################################################
randError<-function(x, y) {
  # function to calculate the adjusted rand statistic
  # x and y are vectors containing the two partitions to be compared
  # first, get crosstabs
  ctab <- table(x,y);

  # now calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2

  # use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1,ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1,nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)
  # now put them together
  adj.rand <- (cellsum - (rowsum*colsum/totsum))/(.5*(rowsum +colsum)-(rowsum*colsum/totsum))
  return (adj.rand);
}

##############################################################
#  transform of an adjacency matrix  into an array of edges  
##############################################################

AdjMat2Edges<-function(x)
  { if (dim(x)[1]==dim(x)[2]){
        # 1) Adjacency matrix: convert the adjacency matrix to the required edge list format  
        aloof<-which(colSums(x)==0)
        if (length(aloof)>0) {
              x<-x[-aloof,-aloof]
              warning("Some nodes are not connected to the network",call. = FALSE)
            }
        nbrNodes<-dim(x)[1]
        m<-t(which((x==1) & (upper.tri(x)),arr.ind=TRUE))
      }
    return(m)
  }




