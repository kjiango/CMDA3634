#Script to run Metropolis-Hastings sampling for various problem sizes and numbers of cores

#Load the parallel package
library(parallel)

#problem size (number of points)
n.pts <- 1e7
#cluster size (number of cores)
ncores <- 16

tm.tot.st = Sys.time()

tm.tot.st = Sys.time()

#function to determine if a point is in the unit circle
in.ucir <- function(x) { as.integer((x[1]^2 + x[2]^2) <= 1) }

n.pts.cl <- ceiling(n.pts / ncores)  # Strong scaling
m = matrix(runif(2*n.pts),n.pts,2)
# n.pts.cl <- n.pts                    # Weak scaling
# m = matrix(runif(2*ncores*n.pts.cl),ncores*n.pts.cl,2)
cat('running parallel version with n.pts = ',n.pts,' and ncores = ',ncores,' (',n.pts.cl,' pts per core) \n')

#start up and initialize the cluster
cl <- makeForkCluster(ncores)
#cl <- makeCluster(ncores)

tm.comp.st = Sys.time()

#Compute in parallel (parallel row apply)
cir = parRapply(cl, m, in.ucir )

#return the proportion of points in the unit circle * 4

pi.approx = 4*mean(cir)

tm.comp = as.numeric(Sys.time() - tm.comp.st, units="secs")

#stop the cluster
stopCluster(cl)

cat('     pi estimate = ',pi.approx,'\n')

pi.err = abs(pi - pi.approx)/pi
cat('  relative error = ',pi.err,'\n')

tm.tot = as.numeric(Sys.time() - tm.tot.st, units="secs")
cat('  computational time = ',tm.comp,'\n')
cat('          total time = ',tm.tot,'\n')
