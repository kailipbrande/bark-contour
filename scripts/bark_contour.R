## R code to calculate bark thickness distributions from a digitized contour

## D.W. Schwilk 2013. Revisions for github repo 2020 in reposnse to recent
## interest in the code.

## (See Adams and Jackson 1995 for description of the method)


library(ggplot2)
#library(grid)
library(extrafont)
#font_import(pattern="Arial")
loadfonts()

## Constants
DATA_FOLDER <- "../data" # assume working directory is location of script in
# ./scripts
RESULTS_FOLDER <- "../results"

CURVE_METHOD <- "smoothSpline"
##PX_PER_MM <- 11.7
SSPAR = 0.7  # spline parameter for smoothing.

### read master summary file. This global data frame is used by several
### function to look up tree diamter, etc.
bark.sum <- read.csv(file.path(DATA_FOLDER, "BarkData.csv"))



## ggplot theme
bestfit <- geom_smooth(method="lm",se = F, color = "black", size=1.5)
isoline <- geom_abline(slope = 1, color="black", size=1.5, linetype="dashed") 
textsize <- 12
smsize <- textsize-2
pt2mm <- 0.35146
smsize.mm <- smsize*pt2mm
fontfamily = "Arial"
themeopts <-   theme(axis.title.y = element_text(family=fontfamily,
                       size = textsize, angle = 90, vjust=0.3),
               axis.title.x = element_text(family=fontfamily, size = textsize, vjust=-0.3),
               axis.ticks = element_line(colour = "black"),
               panel.background = element_rect(size = 1.6, fill = NA),
               panel.border = element_rect(size = 1.6, fill=NA),
               axis.text.x  = element_text(family=fontfamily, size=smsize, color="black"),
               axis.text.y  = element_text(family=fontfamily, size=smsize, color = "black"),
               strip.text.x = element_text(family=fontfamily, size = smsize, face="italic"),
               strip.text.y = element_text(family=fontfamily, size = smsize, face="italic"),
               legend.title = element_text(family=fontfamily, size=textsize),
               legend.text = element_text(family=fontfamily, size=smsize, face="italic"),
               legend.key = element_rect(fill=NA),
               panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(), #element_line(colour = "grey90", size = 0.2),
               strip.background = element_rect(fill = "grey80", colour = "grey50")      
               #panel.grid.major = element_line(colour = NA)
                )


############################################################
## Helper functions
###########################################################

## Euclidian distance
dist <- function(x1,y1,x2,y2){
  return(sqrt( (x1-x2)^2 + (y1-y2)^2))
}

## Intersection of 2 circles
cintersect <- function(x1,y1,x2,y2,r){
 e <- x2-x1
 f <- y2-y1
 p <- sqrt(e^2 + f^2)
 k <- (p^2) / (2*p)
 x3 <- x1 + ((e*k) /p) + ((f/p) * sqrt(r^2 - k^2))
 y3 <- y1 + ((f*k) /p) - ((e/p) * sqrt(r^2 - k^2))
 x4 <- x1 + ((e*k) /p) - ((f/p) * sqrt(r^2 - k^2))
 y4 <- y1 + ((f*k) /p) + ((e/p) * sqrt(r^2 - k^2))
 ## only return the upper (higher y) intersection
 ##print(c(x1,y1,x2,y2,x3,y3,x4,y4))
 if (is.nan(y3)) {
   x <- NA
   y <- NA
 } else if (y3 > y1) {
   x<-x3
   y<-y3
 } else {
    x<-x4
    y<-y4
  }
   return(c(x,y))
}

# x, y: the x and y coordinates of the hull points
# n: the number of points in the curve.
bez <- function(x, y, t)
	{
	outx <- 0
	outy <- 0
	n <- length(x)-1
	for (i in 0:n)
		{
		outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
		outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
		}
 
	return (list(x=outx, y=outy))
	}
bezierCurve <- function(x, y, n=10)
	{
	outx <- NULL
	outy <- NULL
 
	i <- 1
	for (t in seq(0, 1, length.out=n))
		{
		b <- bez(x, y, t)
		outx[i] <- b$x
		outy[i] <- b$y
 
		i <- i+1
		}
 
	return (list(x=outx, y=outy))
	}
  
# Example usage
##  <- c(4,6,4,5,6,7)
## y <- 1:6
## plot(x, y, "o", pch=20)
## points(bezierCurve(x,y,20), type="l", col="red")


### Center of tree bole find circle intersections fo each point on tree bole
### and average as best estimate of center location. x, y and r are vectors.
### One value for each bore depth
bole.center <- function(x,y,r){
  cx <- c()
  cy <- c()
  for(p1 in 1:(length(x)-1)) {
    for(p2 in (p1+1):length(x)) {
      origin <- cintersect(x[p1],y[p1],x[p2],y[p2],r)
      cx <- append(cx,origin[1])
      cy <-append(cy,origin[2])
    }
  }
  return(c(mean(cx),mean(cy)))
}


cart2polar <- function(df) {
   return(data.frame( r = sqrt(df$x^2 + df$y^2), theta = atan2(df$y,df$x))  )
 }

polar2cart <- function(df) {
  return(data.frame(x = df$r * cos(df$theta), y = df$r * sin(df$theta)))
}

local.maxima <- function(x,y,n=5) {
  out.x <- c()
  out.y <- c()
  s <- seq(min(x),max(x),length.out=n+1)
  for(i in seq(1,n) ) {
    j <- which.max(y[x>s[i] & x<s[i+1]])
    out.x <- append(out.x,  x[x>s[i] & x<s[i+1]][j]  )
    out.y <- append(out.y,  y[x>s[i] & x<s[i+1]][j])
  }
  return(list(x = out.x, y = out.y))
}


#################################################################################
###  FUNCTION: get.curves Returns data frame with inner (cambial curve
###             approximation, and outer contour), both truncated to estimable
###             arc
 # Note: we do any unit conversions prior. For now we are assuming all coords and
  # lengths are in same units.  
  # Arguments:
  # - bore.x, bore.y and blengths are three vectors that contain the x,y and
  #   depth information. These obviously must all be in the correct order as the
  #   values really should be lines in a data frame
  # - diam: the diameter of the tree at contour height
  # - contour.x/y: long vector describing the actual outer contour from the
  # - carpenters tool

get.curves <- function(bore.x, bore.y, blengths, diam, contour.x, contour.y,
                       method="bezierPolar") {
  # get tree radius (including bark)
  r <- diam/2.0
  #print(r)
  # trim contour to arc needed based on bores
  rowsneeded = contour.x < max(bore.x) & contour.x > min(bore.x)
  contour.x <- contour.x[rowsneeded]
  contour.y <- contour.y[rowsneeded]
    
  # solve for center of bole -- first rough pass
  #print(bore.x)
  #print(bore.y)
  #plot(bore.x,bore.y)
  origin <- bole.center(bore.x,bore.y,r)
  center.x <- origin[1]
  center.y <- origin[2]
  ## print( center.x)
  ## print( center.y)
  
  ## convert all coordinates to new origin (center of bole)
  bore.x <- bore.x-center.x
  bore.y <- bore.y- center.y
  contour.x <- contour.x - center.x
  contour.y <- contour.y - center.y

     # convert all pts, contour and bores pts to polar coords with origin center
  # of bole
  bore.R <- sqrt(bore.x^2 + bore.y^2)
  bore.theta <- atan2(bore.y,bore.x)
  bore.inner.R <- bore.R - blengths
  contour.R <- sqrt(contour.x^2 + contour.y^2)
  contour.theta <-atan2(contour.y,contour.x)

  #print(bore.R)
  #print(bore.theta)

  ### now get more accurate center:
  cmax <- local.maxima(contour.theta, contour.R , 6)
   # new origin (in cart coords)

  tx <- cmax$y * cos(cmax$x)
  ty <- cmax$y * sin(cmax$x)
  ## print(tx)
  ## print(ty)
  ## print(bore.x)
  ## print(bore.y)
  new.origin <- bole.center(cmax$y * cos(cmax$x), cmax$y * sin(cmax$x), r )
  
  center.x <- new.origin[1]
  center.y <- new.origin[2]
  #print(center.x)
  #print(center.y)

  # transform to new origin:
  bore.x <- bore.x-center.x
  bore.y <- bore.y- center.y
  contour.x <- contour.x - center.x
  contour.y <- contour.y - center.y

    # convert all pts, contour and bores pts to polar coords with origin center
  # of bole
  bore.R <- sqrt(bore.x^2 + bore.y^2)
  bore.theta <- atan2(bore.y,bore.x)
  bore.inner.R <- bore.R - blengths
  contour.R <- sqrt(contour.x^2 + contour.y^2)
  contour.theta <-atan2(contour.y,contour.x)

  ## for graphing:
  bore.inner.x <- bore.inner.R *cos(bore.theta)
  bore.inner.y <- bore.inner.R * sin(bore.theta)


  ## smoothing methods each produces a smooth innner bole based on bore pts and
  ## outputs result in cartesian and polar coords methods smooth on polar
  ## coordinates ath x and y -- seems to work best to enfroce overal
  ## circularity
  if(method == "bezier") {
    t <- bezierCurve(bore.theta,bore.inner.R, length(contour.theta))
    inner.theta <- t$x
    inner.r <- t$y
    # and cartesian version
    inner.x <- inner.r * cos(inner.theta)
    inner.y <- inner.r * sin(inner.theta)  
  } else if (method=="circle") {
    ### method 1: circular bole
    ## ## average inner radius (for  inner as circle rather than loess)
    ## then get inner contour of bores as cartesian coordinates (lx,ly)
    radius.ave <- mean(bore.R-blengths)
    inner.theta = contour.theta
    inner.r = radius.ave
    # and cartesian version
    inner.x <- radius.ave * cos(contour.theta)
    inner.y <- radius.ave * sin(contour.theta)
} else if (method=="spline") {
    fspline <- splinefun(x=bore.theta,y=bore.inner.R)
    inner.theta <- contour.theta
    inner.r <- fspline(inner.theta)
    inner.x <- inner.r * cos(inner.theta)
    inner.y <- inner.r * sin(inner.theta)      
} else if (method=="smoothSpline") {
    fspline <- smooth.spline(x=bore.theta,y=bore.inner.R,spar = SSPAR, all.knots=TRUE)
    inner.theta <- contour.theta
    inner.r <- predict(fspline, inner.theta)$y
    inner.x <- inner.r * cos(inner.theta)
    inner.y <- inner.r * sin(inner.theta)    
  } else {
    print("Unknown method")

  }
  
  
  ndata <- data.frame(outer.theta=contour.theta, outer.r=contour.R,
                      inner.r=inner.r, outer.x = contour.x,
                      outer.y=contour.y, inner.x=inner.x,
                      inner.y=inner.y, bt=contour.R-inner.r)

  bdata <- data.frame(bore.inner.x=bore.inner.x,bore.inner.y=bore.inner.y)
  
  return( list(ndata, bdata, data.frame(x=bore.x, y = bore.y)) )
}

  
## read.contour <- function(fname) {
##     t <- read.table(fname, header=FALSE, col.names = c("X","Y"))
##     t$Y <- 10000 - t$Y
##     t <- t[order(t$X),]
##     t <- ddply(t, .(X), summarize, Y = mean(Y))
##     return(t)
## }

read.contour <- function(fname) {
    t <- read.table(fname, header=FALSE, col.names = c("X","Y"))
    t$Y <- 10000 - t$Y
    t <- t[order(t$X),]
    res <- aggregate(t$Y, list(X=t$X), mean)
    names(res) <- c("X","Y")
    return(res)
}

read.pts <- function(fname) {
    t <- read.table(fname, header=FALSE, col.names = c("X","Y","depth"))
    t$Y <- 10000 - t$Y
    t <- t[order(t$X),]
    return(t)
}


corrected.thickness <- function(df){
  theta.out <- seq(min(df$outer.theta), max(df$outer.theta), length.out=500)
  thick.fun <- approxfun(df$outer.theta, df$bt, rule=2)
  return(thick.fun(theta.out))
}


## Create contour for tree ID "id" looking diameter and depths in the bark
## summary data table and finding the contour points in the coord-data
## directory. Returns a list as porduced by get.curves()
makeContour <- function(id, method=CURVE_METHOD) { 
  diam <- 10 * bark.sum$D60[bark.sum$TreeID==id]  # convert to mm from cm
  res <- bark.sum$Pix.mm[bark.sum$TreeID==id]

  pfile = file.path(DATA_FOLDER, "coord-data", paste(id, "pts.txt", sep=""))
  cfile = file.path(DATA_FOLDER, "coord-data", paste(id, "xy.txt", sep=""))

  if (file.exists(pfile) & file.exists(cfile)) {
    print(id)
    points <- read.pts(pfile)
    contour <- read.contour(cfile)
    ## correct units
    cx <- contour$X / res
    cy <- contour$Y / res
    px <- points$X / res
    py <- points$Y / res

    # Now hand this over the to get.curves function
    result <- get.curves(px,py, points$depth, diam,cx,cy, method=method)
  } else {
    warning(paste("No file found for tree ID:", id))
    result <- NULL
  }
  return(result)
}


# return a vector of bark thicknesses for the get.curve object
getThicknesses <- function(x) {
    thick.df <- x[[1]]
    bore.df <- x[[2]]
    ma <- function(x,n=3){stats::filter(x,rep(1/n,n), sides=2)}
    thick.df$outer.y2 <- ma(thick.df$outer.y)
    thicknesses <- corrected.thickness(thick.df)
  return(thicknesses)
}


# make a pretty figure of a bark contour. Argument 'x' is an object as returned
# by the get.curves function.
makeContourPlot <- function(x) {
  ma <- function(x,n=3){stats::filter(x,rep(1/n,n), sides=2)}
  thick.df <- x[[1]]
  bore.df <- x[[2]]
  outpoints <- x[[3]]
  thick.df$outer.y2 <- ma(thick.df$outer.y)

  boresegs <- data.frame(xo = outpoints$x,
                         yo = outpoints$y,
                         xi = bore.df$bore.inner.x,
                         yi = bore.df$bore.inner.y)
  
  bark.polygon <- data.frame(x = c(thick.df$outer.x, thick.df$inner.x),
                             y = c(thick.df$outer.y, thick.df$inner.y))
  
  result <- ggplot(thick.df, aes(x=outer.x, y=outer.y2)) +
    geom_segment(aes(x=xo,y=yo,xend=xi,yend=yi), data = boresegs,
                 size = 1,
                 alpha = 0.9,
                 arrow = arrow(length = unit(0.4,"cm"))) +
    geom_line(aes(inner.x, inner.y), size=1, linetype="dashed") + # smoothed inner bole
    geom_line(size=1)  +  # outer bark contour line
    #      geom_point(aes(bore.inner.x, bore.inner.y),data=bore.df,  size=3) +
    scale_x_continuous("X dimension (mm)") +
    scale_y_continuous("Y dimension (mm)") +
    coord_fixed(ratio=1) + themeopts

  return(result)
}


###############################################################################
#### main script  
###############################################################################

# create empty results data frame
results.df <- data.frame(TreeID = NULL, bark.thicknesses = NULL)

# go through each tree id in bark.sum (read in at top of script)
for( id in bark.sum$TreeID) {
  contour <- makeContour(id)
  if(is.null(contour)) {
    print(paste("File missing for id", id))
  } else {
    thicknesses <- getThicknesses(contour)
    t.df <- data.frame(TreeID = id, bark.thicknesses = thicknesses)
    results.df <- rbind(results.df,t.df)
  }
}
## now write results.df to a csv file if you want to save the info
write.csv(results.df, file = file.path(RESULTS_FOLDER, "barkthick.csv"),
          row.names = FALSE)

## Now we can summarize that thickness data if we wish:
library(dplyr)
results.sum <- summarize(group_by(results.df, TreeID),
                         bark.mean = mean(bark.thicknesses),
                         bark.sd = sd(bark.thicknesses),
                         q10 = quantile(bark.thicknesses,0.10, na.rm=TRUE),
                         q90 = quantile(bark.thicknesses,0.90, na.rm=TRUE))

  
# and make example figure
makeContourPlot(makeContour("QEMO515"))
