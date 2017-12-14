# File name: mcpi_pbdr.r
# Description: Calculate pi using Monte Carlo
# Run: mpiexec -np 8 Rscript --vanilla mcpi_pbdr.r

# Load pbdMPI and initialize the communicator
library(pbdMPI, quiet = TRUE)
init()

#Load mcpi() function (on each rank)
source("mcpi.r")

# Number of total points to use
n.pts <- 1e7
# Number of points per rank
n.ranks <- comm.size()
#n.pts <- n.pts/n.ranks
n.pts <- n.pts*n.ranks
n.pts.rank <- ceiling(n.pts / n.ranks)

comm.print(paste('running pbdMPI version with n.pts = ',n.pts,' and nranks = ',n.ranks,' (',n.pts.rank,' pts per rank) \n'))

tm.comp.st = Sys.time()

# Compute pi() on each rank
pi.rank <- mcpi(n.pts.rank)

# Use reduce() to total across ranks
pi.approx <- reduce(pi.rank, op="sum") / n.ranks

tm.comp = as.numeric(Sys.time() - tm.comp.st, units="secs")
pi.err <- abs(pi - pi.approx)/pi

# Output from RANK 0 since reduce() will dump only to 0 by default
comm.print(paste('approximate value of pi is: ',pi.approx,'\n'))
comm.print(paste('            relative error: ',pi.err,'\n'))

comm.print(paste('  computational time = ',tm.comp,'\n'))

# Wrap up
finalize()

