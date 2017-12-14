vecadd.loop <- function(x,y) {
  l = length(x)
  z = rep(NA,l)
  for (i in 1:l) {
    z[i] = x[i] + y[i]
  }
  return(z)
}

vecadd.array <- function(x,y) {
  return(x+y)
}

x = rnorm(1e6); y = rnorm(1e6)

cat("looped vector add:\n")
z <- vecadd.loop(x,y)
system.time(z <- vecadd.loop(x,y))

cat("array vector add:\n")
z <- vecadd.array(x,y)
