## __________________________________________________________
##
## Gplot
##
## Public
##
## author : A. Smith and G. Grasseau
##
## INPUT : 
##	Matrix   : graph precision matrix
##
##	cl        : node classification vector 
##	cols      : colours to use (per class if cl is specified) 
##	coord	  : matrix of node coordinates for gplot
##	my.labels  : vector of node labels (to override dimnames(PrecMat))
##	display.labels : boolean : display labels?
##
##	degree.threshold : threshold under which degrees are considered null ;
##		null-degree nodes are considered as dutbin nodes
##	plot.dustbin.nodes : boolean : plot dustbin nodes?
##
##	max.edges : max number of edges below which the graph is plotted
##
##	main	  : graphical parameter main
##	sub       : graphical parameter sub
##
## OUTPUT :
##	gplot matrix of node coordinates
##	gplot graphic
##
## Wrapper for gplot in our case :
## Our graphs have real-valued undirected edges, have 
## coloured vertices, and don't have self-loops.
## __________________________________________________________
##
Gplot <- function ( X, 
		    cl=NULL,
		    ... ) {
 
  #    Defaults for hidden optional arguments 
  #    --------------------------------------
  
  #  Colours for each class
  cols			<- sub.param("cols"		, NULL	, ...)
  
  #  Coordinates for nodes
  coord			<- sub.param("coord"		, NULL	, ...)
  
  #  Labels for nodes                                       
  my.labels		<- sub.param("my.labels"	, NULL	, ...)
  
  #  Display labels ?
  display.labels	<- sub.param("display.labels"	, !is.null(my.labels),	
                                     ...)
  # Degree threshold under which edges are considered as
  degree.threshold	<- sub.param("degree.threshold"	, .Machine$double.xmin,	
                                     ...)
  # Plot dustbin nodes (isolated nodes)
  plot.dustbin.nodes	<- sub.param("plot.dustbin.nodes", FALSE, ...)
  
  # Max. number of plottable edges
  max.edges		<- sub.param("max.edges"	, 10000	, ...)
  
  # Graph title
  main			<- sub.param("main"		, "Gplot", ...)
  
  # Graph subtitle
  sub			<- sub.param("sub"		, NULL	, ...)	


  if (sum(abs(X))==0) {
    cat("Gplot warning : no edges to plot \n")
    return(NULL)
  }

  # Edge matrix case
  if ( dim(X)[1] == 2) {
    if ( dim(X)[2] == 2 ) {
      cat("Gplot warning : Matrix(2,2), treated as an adjacency matrix \n")
    } else {
        NNodes <- max(X)
        Y <- matrix(0, NNodes, NNodes )
        Y[ cbind( X[1,], X[2,] ) ] <- 1
        Y[ cbind( X[2,], X[1,] ) ] <- 1
        X <- Y
    }
  }
  # Dimension check
  p <- dim(X)[1]
  if (p != dim(X)[2]) {
    cat("Gplot warning : X matrix must be square \n")
    return(coord)
  }
  old.p <- p
  
  if (! isSymmetric(X)) {
    cat("Gplot warning : X matrix must be symmetric \n")
    cat("\t automatically symmetrizing using weak rule. \n")
    X <- Symmetrize(X)
  }


  # Identify dustbin nodes
  b.dust <- rep(FALSE, p)
  if (!is.null(cl)) {
    b.dust[cl=="D"] <- TRUE
  } else {
    b.dust  <- ( CalcDegrees(abs(X)) <= degree.threshold )
  }

  # If plot.dustbin.nodes == FALSE, we'll just drop those nodes!
  # (gplot handles dustbin nodes poorly)
  if (!plot.dustbin.nodes) {
    X <- X[!b.dust,!b.dust]
    p <- dim(X)[1]
    if (!is.null(coord)) {
      coord <- coord [!b.dust,]
    }
    if (!is.null(cl)) {
      cl <- cl [!b.dust]
    }
    if (!is.null(my.labels)) {
      my.labels <- my.labels [!b.dust]
    }
  } 


  # gplot threshold and line weight handling
  e.col <- X
  e.col[X<0] <- "red"  # red : negative partial correlation (inhibition)
  e.col[X>0] <- "blue" # blue : positive partial correlation (activation)
  e.col[X==0]<- "black"	

  lwd <- abs(X)
  lwd <- lwd * 5 / max( lwd )	        # max line weight is 5
	
  e.lty <- lwd
  e.lty[abs(lwd)>1]  <- "solid"
  e.lty[abs(lwd)<=1] <- "dotted"
  e.lty[abs(lwd)==0] <- "blank"


  # Get node labels
  if (is.null(my.labels)) {
   my.labels <- dimnames(X)[[1]]
  }


  # Get node colours
  if (is.null(cl)) {
    if (is.null(cols)) {
      vertex.col <- rep(1, p)
    } else {
      vertex.col <- cols
    }
  } else {
    cl <- as.factor(cl)
    if (is.null(cols)) {
      cols <- rainbow(length(levels(cl)))
    }
    vertex.col <- cols[as.numeric(cl)]
  }
  if (plot.dustbin.nodes) {
    vertex.col[b.dust] <- NA
  }

  # Strictly undirected graph, matrix considered symmetric
  X <- abs(X)
  X[lower.tri(X)] <- 0

  # Only plot if not too many edges
  num.edges <- sum(X>.Machine$double.xmin)
  if (num.edges>max.edges) {

    cat( "Gplot warning : too many edges to plot (",
         num.edges,">",max.edges,") \n")
    g  <- coord
    g2 <- g

  } else {
    gr <- Gplot.graphics( X,
                          thresh		= 0,	    # no edge blocking
                          margin		= 0.1,
                          main		= main,
                          sub		= sub
                        )
    gr <- Gplot.network( gr,
                         drawloops       = FALSE,
                         displayisolates = plot.dustbin.nodes,
                         coord           = coord
                       )
    gr <- Gplot.edge( gr,
                      col       = e.col,
                      lty       = e.lty,
                      lwd       = lwd,
                      usearrows = FALSE
                    )
    gr <- Gplot.vertex( gr, col= vertex.col )

    if ( display.labels )
      gr <- Gplot.vertex.label( gr,
                                label=my.labels,
                                cex=0.7,
                                useboxes=FALSE
                               )
    if (plot.dustbin.nodes) {
      g2 <- gr$network$coord
    } else {
      g2 <- matrix(0, old.p, 2)
      g2[!b.dust,] <- gr$network$coord
    }
  }

  invisible(g2)

}



## __________________________________________________________
##
## CalcDegrees
##
## Internal
##
## author : A. SMITH
##
## INPUT : 
##	K : precision matrix
##	MARGIN : change this to get in- and out- degrees
##		note : on a symmetric matrix, 
##		makes no difference
## OUTPUT :
##	vector of node degrees
## __________________________________________________________
##
CalcDegrees <- function ( K, MARGIN=1 ) {

	K <- as.matrix(K)
	diag(K) <- 0
	return( apply( abs(K), MARGIN, sum ) )
	
}

Symmetrize <- function (X) {

  n <- dim(X)[1]
  for (j in 1:n) X[j:n, j] <- X[j, j:n]

  return(X)
}

sub.param <- function (param, default=NULL, ... ) {

  if (missing(param)) { return(default) }

  l <- list(...)
  res <- l[[  which( names(l) %in% param)[1] ]] 

  if (is.null(res)) { res <- default }

  return ( res )

}
