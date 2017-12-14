#function to calculate pi by comparing the area of the unit
#circle to that of the unit sphere (should be pi/4)
# inputs:
#   n.pts    = number of points to test
#
mcpi <- function(n.pts) {
  #generate n.pts (x,y) points in the unit square
  m = matrix(runif(2*n.pts),n.pts,2)
  #determine if they are in the unit circle
  in.ucir <- function(x) { as.integer((x[1]^2 + x[2]^2) <= 1) }
  in.cir = apply(m, 1, in.ucir )
  #return the proportion of points in the unit circle * 4
  return (4*sum(in.cir)/n.pts)
}
