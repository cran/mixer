##########################################################################
#  
#  Routines to draw graphs.
#
#  These routines have been extracted from the SNA package (mainly gplot)
#  and reorganized to exploit more conveniently its possibilies.
#
#  Author : G.Grasseau and A.Smith (Statistics and Genome laboratory)
#
#  Entry points:
#  ------------
#
#    Gplot.graphics:     set default graphic parameters.
#                        Main steps:
#                        + store initial matrix as boolean one,
#                        + build the structure (list) which handle all
#                        graphics parameters.
#
#    Gplot.network:      specify the network representation (up to now only
#                        undirect graphs are taken into account).
#                        Main steps:
#                        + compute vertex position,
#                        + graph box limits and scaling factor,
#                        + define which vertices to display.
#
#    Gplot.vertex:       draw vertices.
#
#    Gplot.vertex.label: draw vertex labels. 
#
#    Gplot.edge:         draw edges
#                        Main steps:
#                        + identify edges to be drawn.
#                        + remove loops
#                        + draw edges/arrows.
#
#  Utilities:
#  ---------
#   + draw.edge:        (used by Gplot.edge) draw edges/arrows
#   + isolate.vertices: (used by Gplot.network) tag isolate vertices (TRUE)
#
#  To improve:
#  ----------
#   + sides of some boxed labels are not displayed
#   + 2 variables (displayisolates and use.isolates) should be introduced
#     to avoid taking into account of isolate vertices in the vertex coordinate
#     computation (4 cases have to be studied). 
#
#  To implement:
#  ------------
#   + drawing the loops
#
##########################################################################

Gplot.graphics<-function( mat, thresh=0, xlim=NULL, ylim=NULL, scale=0.01,
                          margin=0.2, main="", sub=""
                        ) {
# -------- Arguments -----------------------------------------------------
#
# mat        (matrix): initial matrix
# thresh     (scalar): matrix values are set to zero if they are < threshold
#                      (Pb: should be "abs( mat ) < 0")
# xlim, ylim (vector): x minimun and x maximum (same for y)
# scale  (scalar)    : scaling factor for the whole graphics.
# margin (scalar)    : margin value to add to all the graphics box sides.
# main  (char): title  
# sub   (char): subtitle
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters (default)
#
# ------------------------------------------------------------------------
  
   n<-dim(mat)[1]

   # Replace NAs with 0s
   mat[is.na(mat)]<-0
   
   # Save a copy of mat
   mat.raw<-mat
   
   # Binary matrix
   mat<-matrix(as.numeric(mat>thresh),n,n)

   l.graphics <- list( xlim=xlim, ylim=ylim, scale=scale, margin=margin,
                       main=main, sub=sub ) 

   l.network <- list( mode=NULL, drawloops=FALSE, vertex.pos.mode=NULL,
                      coord= NULL, displayisolates=TRUE, use=NULL )
   
   l.label   <- list( label=c(1:dim(mat)[1]), cex=1, col=1, pos=0,
                      useboxes=TRUE, box.margin=0.5, box.col=1,
                      box.bg="white", box.lty=NULL, box.lwd=par("lwd") )

   l.vertex <- list( cex=1, sides=8, col=2, border=1, lty=1 ,
                     label = l.label, baserad=0 )

   l.edge  <- list( col=1, lty=1, lwd=0,
                    usearrows=TRUE, arrow.cex=1, loop.cex=1 )

   graph   <- list( mat=mat, thresh=thresh, mat.raw = mat.raw,
                    graphics=l.graphics, network=l.network,
                    vertex=l.vertex, edge=l.edge )

   return( graph )
}


Gplot.network <-function( graph,
                          mode=NULL, drawloops=FALSE, vertex.pos.mode=NULL,
                          coord= NULL,
                          displayisolates=TRUE, ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph (list): structure handling all graphics parameters
# mode  (char): kind of network to draw. 
# drawloops       (bool): draw network loops (not implemented).
# vertex.pos.mode (char): vertex positionning mode (not used).
# coord         (vector): vertex position.
#                         If NULL these coordinates are computed
# displayisolates (bool): display isolate vertices.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------

  # Arguments mode and vertex.pos.mode are not implemented
  if( ! missing(drawloops) )
    graph$network$drawloops = drawloops
  if( ! missing(coord) )
    graph$network$coord = coord
  if( ! missing(displayisolates) )
    graph$network$displayisolates = displayisolates

  if(is.null(graph$network$coord)){      
    
    # Provide default settings
    n           <- dim(graph$mat)[1]
    niter       <- 500
    max.delta   <- n
    area        <- n^2
    cool.exp    <- 3
    repulse.rad <- area*n
    
    # Set initial positions randomly on the circle
    tempa<-sample((0:(n-1))/n)
    x<-n/(2*pi)*sin(2*pi*tempa)
    y<-n/(2*pi)*cos(2*pi*tempa)

    layout<-.C( "vertex_coord_C",
                as.integer(graph$mat),
                as.double(n), as.integer(niter),
                as.double(max.delta), as.double(area),
                as.double(cool.exp), as.double(repulse.rad),
                x=as.double(x), y=as.double(y),
                PACKAGE="mixer"
               )

    graph$network$coord <- cbind(layout$x,layout$y)
  }
  
  # Remove isolated vertex (if displayisolates FALSE)
  use <- displayisolates | ( ! isolate.vertices( graph$mat ) ) 

  graph$network$use = use
  
  x <- graph$network$coord[,1]
  y <- graph$network$coord[,2]

                         
  # Set limits for plotting region
  xlim = graph$graphics$xlim                
  ylim = graph$graphics$ylim
  margin = graph$graphics$margin
  if(is.null(xlim))
    xlim<-c(min(x[use])-margin,max(x[use])+margin)  # Save x, y limits
  if(is.null(ylim))
    ylim<-c(min(y[use])-margin,max(y[use])+margin)
  xrng<-diff(xlim)          
  yrng<-diff(ylim)
  xctr<-(xlim[2]+xlim[1])/2                 # Get center of plotting region
  yctr<-(ylim[2]+ylim[1])/2
  
  # Force scale to be symmetric
  if(xrng<yrng)
    xlim<-c(xctr-yrng/2,xctr+yrng/2)
  else
    ylim<-c(yctr-xrng/2,yctr+xrng/2)

  graph$graphics$xlim = xlim                
  graph$graphics$ylim = ylim               

  # Extract "base radius"
  graph$vertex$baserad <- min(diff(xlim),diff(ylim))* graph$graphics$scale
  
  # Configure the graphic box
  plot( 0,0,
        xlim=graph$graphics$xlim,
        ylim=graph$graphics$ylim, type="n", xlab="",ylab="", asp=1, axes=FALSE,
        main= graph$graphics$main, sub=graph$graphics$sub,
        ...
       )
  return( graph )
}


Gplot.vertex <-function( graph,
                         cex, sides, col, border, lty, ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# cex    (scalar/vector): scaling factors.
# sides  (scalar/vector): number of node sides.
# col    (scalar/vector): vertex colors.
# border (scalar/vector): color of node (vertex) borders.
# lty    (scalar/vector): type of nodes borders.
#                         (Pb: would be better with the line width
#                         border parameter "lwd".
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters 
#
# ------------------------------------------------------------------------

  if( ! missing(cex) )
   graph$vertex$cex    = cex
  if( ! missing(sides) )
   graph$vertex$sides  = sides
  if( ! missing(col) )
   graph$vertex$col    = col
  if( ! missing(border) )
   graph$vertex$border = border
  if( ! missing(lty) )
   graph$vertex$lty    = lty

   n <- dim(graph$mat)[1]
  
   # Build vectors describing vertex
   v.cex    <- rep( graph$vertex$cex                       , length=n)
   v.radius <- graph$vertex$baserad * v.cex
   v.sides  <- rep( graph$vertex$sides                     , length=n)
   v.col    <- rep( graph$vertex$col                       , length=n)
   v.border <- rep( graph$vertex$border                    , length=n)
   v.lty    <- rep( graph$vertex$lty                       , length=n)

   # remove unused
   use = graph$network$use
   v.radius <-  v.radius[use]
   v.sides  <-  v.sides[use]
   v.col    <-  v.col[use]
   v.border <-  v.border[use]
   v.lty    <-  v.lty[use]

   x <- graph$network$coord[use,1]
   y <- graph$network$coord[use,2]   
   n <- length(x)

 # Compute the coordinates
  coord<-vector()

  for(i in 1:n){
    ang <- (1:v.sides[i])/v.sides[i]*2*pi
    dx <- v.radius[i]*cos(ang)
    dy <- v.radius[i]*sin(ang)
    XY = rbind( cbind( x[i]+dx, y[i]+dy ), c(NA,NA) )
    coord<-rbind(coord, XY)
  }
  # Plot the vertices
  polygon(coord, col=v.col, border=v.border, lty=v.lty, ...)

  return( graph )
}


Gplot.vertex.label <-function( graph,
                               label, cex, col, pos,
                               useboxes, box.margin, box.col, box.bg,
                               box.lty, box.lwd,
                               ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# label         (vector): label titles (if NULL a default is provided)
# cex    (scalar/vector): scaling factors.
# col    (scalar/vector): label colors.
# pos           (scalar): label positionning mode
#                         + 0 labels are placed away from the graph
#                         + 1 labels are placed below the vertices      
#                         + 2 labels are placed on the vertex left.      
#                         + 3 labels are placed above the vertices      
#                         + 4 labels are placed on the vertex right.
# useboxes        (scalar): frame (box) the labels
# box.margin      (scalar): margin between the label titles and their boxes
#                           (in character size unit).
# box.col  (scalar/vector): box colors.
# box.bg   (scalar/vector): box backgroung color.
# box.lty  (scalar/vector): boxe line type.
# box.lwd  (scalar/vector): boxe line width.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------
  if( ! missing( label ) )
    if( ! is.null( label ) )
      graph$vertex$label$label = label
  if( ! missing( cex ) )
   graph$vertex$label$cex   = cex
  if( ! missing( col ) )
   graph$vertex$label$col   = col
  if( ! missing( pos ) )
   graph$vertex$label$pos   = pos
  if( ! missing( useboxes ) )
   graph$vertex$label$useboxes   = useboxes
  if( ! missing( box.margin ) )
   graph$vertex$label$box.margin = box.margin
  if( ! missing( box.col ) )
   graph$vertex$label$box.col    = box.col
  if( ! missing( box.bg ) )
   graph$vertex$label$box.bg     = box.bg
  if( ! missing( box.lty ) )
   graph$vertex$label$box.lty    = box.lty
  if( ! missing( box.lwd ) )
   graph$vertex$label$box.lwd    = box.lwd

  # Plot vertex labels
  use  <- graph$network$use
  x <- graph$network$coord[use,1]
  y <- graph$network$coord[use,2]
   
  if((!all(graph$vertex$label$label==""))&(!all(use==FALSE))){

    # Label display mode
    if ( graph$vertex$label$pos == 0 ){
      
      # Labels are placed away from the graph
      xoff <- x - mean(x)
      yoff <- y - mean(y)
      roff <- sqrt(xoff^2+yoff^2)
      xhat <- xoff/roff
      yhat <- yoff/roff

    } else if (graph$vertex$label$pos<5) {
      
      # below (0,-1), left (-1,0), top (0,1) , right (1,0)
      xhat <- switch( graph$vertex$label$pos,  0,-1, 0, 1)
      yhat <- switch( graph$vertex$label$pos, -1, 0, 1, 0)

    } else {
      xhat <- 0
      yhat <- 0
    }
     
    # Get character size
    l.cex = graph$vertex$label$cex
    char.len <- par()$cxy * l.cex

    # Get label width and height
    label <- graph$vertex$label$label[use]
    lw <- strwidth( label,cex=l.cex) / 2
    lh <- strheight(label,cex=l.cex) / 2
    b.size   = 1 + graph$vertex$label$box.margin
    
    v.radius <- graph$vertex$baserad * rep( graph$vertex$cex,
                                            dim(graph$mat)[1] )

    # Draw boxes
    if( graph$vertex$label$useboxes ){
      rect(x - lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
           y - lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
           x + lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
           y + lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
           col    = graph$vertex$label$box.bg,
           border = graph$vertex$label$box.border,
           lty    = graph$vertex$label$box.lty,
           lwd    = graph$vertex$label$box.lwd)
    }

    # Draw labels
    text(x + xhat * ( lw*(b.size+0.2) + v.radius ),
         y + yhat * ( lh*(b.size+0.2) + v.radius ),
         label, cex=l.cex, col=graph$vertex$label.col, offset=0, ...)
    
  }
  return ( graph )
}

Gplot.edge <-function( graph,
                       col=1, lty=1, lwd=0,
                       usearrows=TRUE, arrow.cex=1, loop.cex=1,
                       ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph              (list): structure handling all graphics parameters
# col  (scalar/vector/matrix): edge colors.
# lty  (scalar/vector/matrix): edge line types.
# lwd  (scalar/vector/matrix): edge line widths.
# usearrows          (bool): draw arrows.
# arrow.cex (scalar/vector): arrow head scaling factors.
# loop.cex  (scalar/vector): loop scaling factors.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------
  
   # Remark : arrow.cex and loop.cex could be fused

   if ( ! missing(col) )
     graph$edge$col = col
   if ( ! missing(lty) )
     graph$edge$lty = lty
   if ( ! missing(lwd) )
     graph$edge$lwd = lwd
   if ( ! missing(usearrows) )
     graph$edge$usearrows = usearrows
   if ( ! missing(arrow.cex) )
     graph$edge$arrow.cex = arrow.cex
   if ( ! missing(loop.cex) )
     graph$edge$loop.cex  = loop.cex

   n <- dim( graph$mat )[1]

   # Build vectors describing edges
   # Each edge is a polygon
   px0<-vector()   
   py0<-vector()
   px1<-vector()
   py1<-vector()
   
   e.lwd<-vector()  # Create edge attribute vectors
   e.type<-vector()
   e.col<-vector()
   e.hoff<-vector() # Offset radii for heads
   e.toff<-vector() # Offset radii for tails
   e.diag<-vector() # Indicator for self-ties

   # Coerce edge.col/edge.lty to array form
   if(!is.array(graph$edge$col))   
     col <- array(graph$edge$col, dim=dim(graph$mat))
   else
     col <- graph$edge$col
   if(!is.array(graph$edge$lty))
     lty<-array(graph$edge$lty, dim=dim(graph$mat))
   else
     lty = graph$edge$lty
   if(!is.array( graph$edge$lwd)){
     if( graph$edge$lwd>0 )
       lwd<-array( graph$edge$lwd * graph$mat.raw, dim=dim(graph$mat))
     else
       lwd<-array(1, dim=dim(graph$mat))
   }
   
   v.radius <- graph$vertex$baserad * rep( graph$vertex$cex, dim(graph$mat)[1] )

   # Select edges between vertices
   # -----------------------------
   x <- graph$network$coord[,1]
   y <- graph$network$coord[,2]

   for(i in (1:n)[graph$network$use]) {   
     for(j in (1:n)[graph$network$use]) {
       if( graph$mat[i,j] ){        # Edge exists
         px0 <- c(px0,as.real(x[i]))  # Store endpoint coordinates
         py0 <- c(py0,as.real(y[i]))
         px1 <- c(px1,as.real(x[j]))
         py1 <- c(py1,as.real(y[j]))
         e.toff <-c ( e.toff, v.radius[i] ) # Store endpoint offsets
         e.hoff <-c ( e.hoff, v.radius[j] )
         e.col  <-c ( e.col , col[i,j])     # Store other edge attributes
         e.type <-c ( e.type, lty[i,j])
         e.lwd  <-c ( e.lwd , lwd[i,j])
         e.diag <-c ( e.diag, i==j)         # Set to true if diagonal 
       }
     }
   }

   # Remove loops
   # ------------
   if(length(px0)>0){
     px0 <- px0[!e.diag] 
     py0 <- py0[!e.diag]
     px1 <- px1[!e.diag]
     py1 <- py1[!e.diag]
     e.lwd  <- e.lwd[!e.diag]
     e.type <- e.type[!e.diag]
     e.col  <- e.col[!e.diag]
     e.hoff <- e.hoff[!e.diag]
     e.toff <- e.toff[!e.diag]
   }

   # Draw edges
   if(length(px0)>0)
     draw.edges(
            as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1),
            width=e.lwd*graph$vertex$baserad/10, col=e.col, lty=e.type,
            o.head=e.hoff, o.tail=e.toff,
            arrow=usearrows, a.len=2*graph$vertex$baserad*arrow.cex, a.angle=20,
            ...
          )

  return( graph )
}
  
isolate.vertices<-function( mat. ) {
# -------- Arguments -----------------------------------------------------
#
# mat. (matrix): binary (0/1) square matrix
#
# -------- Return value --------------------------------------------------
#
# isolate (vector): isolated (if TRUE) vertex vector.
#
# ------------------------------------------------------------------------
  
  n <- dim(mat.)[1];
  mat <- mat.
  if ( n > 1 ){
    # Set to zero NA and diagonal terms
    for(i in 1:n){
      mat[i,i] = 0
      mat[i,1:n] ==  as.numeric( ! is.na(mat[i,1:n]) )
    }
    isolate <- vector()
    for(i in 1:n) {
      isolate = c( isolate, all(( mat[i,] == 0 )) & all(( mat[,i] == 0 )) )
    }
  }
  return( isolate )
}

draw.edges<-function( x0, y0, x1, y1,
                      width=0.01, col=1, lty=1,
                      o.head=0, o.tail=0,
                      arrow=TRUE, a.len=0.1, a.angle=2,
                      ... )
{
# -------- Arguments -----------------------------------------------------
#
# x0, y0       (vector): start coordinates of edges to draw.
# x1, y1       (vector): end coordinates of edges to draw.
# width (scalar/vector): edge line widths.
# col   (scalar/vector): edge colors.
# lty   (scalar/vector): edge line types.
# o.head (scalar/vector): offset (vertex size shift) at the start points.
# o.tail (scalar/vector): offset (vertex size shift) at the end points.
# arrow           (bool): arrows are drawn.
# a.len   (scalar/vector): arrow head lengths.
# a.angle (scalar/vector): arrow head angle (in degree).
#  
# -------- Return value --------------------------------------------------
#
# No value
#
# ------------------------------------------------------------------------

  if(length(x0)==0)   #Leave if there's nothing to do
    return;

  n<-length(x0)

  # Transform scalars into vectors
  width <- rep(width,length=n)
  col   <- rep(col,length=n)
  lty   <- rep(lty,length=n)
  
  # Offsets
  o.head  <- rep(o.head,length=n) 
  o.tail  <- rep(o.tail,length=n)

  # Arrow parameters
  a.angle <- rep(a.angle,length=n)/360*2*pi
  a.len   <- rep(a.len,length=n)

  # Debug point
  # cat("xy  :",x0, y0, x1,y1, "\n")
  # cat("width :", width, "\n")
  # cat("col :", col, "\n")
  # cat("lty :", lty, "\n")
  # cat("o.tail :",o.tail, "\n")
  # cat("o.head :", o.head , "\n")
  # cat("a.angle :", a.angle, "\n")
  # cat("a.len   :", a.len, "\n")

  # Computes edges/arrows coordinates
  coord<-vector()
  XY <- vector()
  for(i in 1:n) {  

    slen<-sqrt((x0[i]-x1[i])^2+(y0[i]-y1[i])^2)  #Find the total length

    if(arrow){
      
      #  With Arrows
      a.sin = sin( a.angle[i] )
      a.cos = cos( a.angle[i] )
      XY<-rbind(                    
        c( - width[i]/2      , o.tail[i]),
        c( - width[i]/2      , slen - 0.5*a.len[i] - o.head[i]),
        c( - a.len[i] * a.sin, slen - a.len[i]*a.cos - o.head[i]),
        c(   0               , slen - o.head[i]),
        c(   a.len[i] * a.sin, slen - a.len[i]*a.cos - o.head[i]),
        c(   width[i]/2      , slen - 0.5*a.len[i] - o.head[i]),
        c(   width[i]/2      , o.tail[i] ),
        c(   NA              , NA)
      )
    }else{
        
      #  Without Arrows
      XY<-rbind(                    
                 c( - width[i]/2, o.tail[i]       ),
                 c( - width[i]/2, slen - o.head[i]),
                 c(   width[i]/2, slen - o.head[i]),
                 c(   width[i]/2, o.tail[i]       ),
                 c(   NA,      NA)
                )
    }
    # Rotate
    theta <- atan2(y1[i]-y0[i],x1[i]-x0[i])-pi/2     
    rmat  <- rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
    XY    <- XY %*% rmat
    # Translate
    XY[,1] <- XY[,1]+x0[i]            
    XY[,2] <- XY[,2]+y0[i]

    coord<-rbind( coord, XY)
  }
  
  # GG - ???
  #coord<-coord[-NROW(coord),]

  #Draw polygons
  polygon(coord,col=col,border=col,lty=lty,...)
}
