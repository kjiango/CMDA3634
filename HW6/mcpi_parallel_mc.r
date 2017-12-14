#Script to run Metropolis-Hastings sampling for various problem sizes and numbers of cores

#Load the parallel package
library(parallel)

#problem size (number of points)
n.pts <- 1e7
#cluster size (number of cores)
ncores <- 16

tm.tot.st = Sys.time()

n.pts.cl <- ceiling(n.pts / ncores)  # Strong scaling
cat('running parallel version with n.pts = ',n.pts,' and ncores = ',ncores,' (',n.pts.cl,' pts per core) \n')

tm.comp.st = Sys.time()

#Load the function to calculate pi
source('mcpi.r')

#call mcpi in parallel
n.pts.cl = rep(n.pts.cl,ncores) #array of inputs to mcpi
pi.approx = mean( unlist( mclapply(n.pts.cl, mcpi, mc.cores=ncores) ) )

tm.comp = as.numeric(Sys.time() - tm.comp.st, units="secs")

cat('     pi estimate = ',pi.approx,'\n')

pi.err = abs(pi - pi.approx)/pi
cat('  relative error = ',pi.err,'\n')

tm.tot = as.numeric(Sys.time() - tm.tot.st, units="secs")
cat('  computational time = ',tm.comp,'\n')
cat('          total time = ',tm.tot,'\n')

