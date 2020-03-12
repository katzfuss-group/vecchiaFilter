setwd("~/HVLF")
source('aux-functions.r')
resultsDir = "simulations-linear"



## define the temporal evolution function
evolAdvDiff = function(state, adv=0, diff=0){
  
  # we assume that there is the same number of grid points
  # along each dimension
  N = dim(as.matrix(state))[1]
  Ny = Nx = sqrt(N)
  
  dx = dy = 1/Nx
  d = diff/(dx**2)
  
  c1 = 1 + 2*(d + d) - adv*(1/dx + 1/dy)
  c2 = - d + adv*(1/dy)
  c3 = - d
  
  diags = list(rep(c2, N-Nx), rep(c2, N-1), rep(c1, N), rep(c3, N-1), rep(c3,N-Nx) )
  E = Matrix::bandSparse(N, k=c(-Nx, -1, 0, 1, Nx), diag=diags)
  
  if(class(state)=='matrix' || class(state)=='dgCMatrix') return( E %*% state )
  else as.numeric(E %*% as.matrix(state))
}



######### set parameters #########
set.seed(1996)
spatial.dim = 2
n = 300**2
m = 50

frac.obs = 0.05
Tmax = 100
diffusion = 0.0000001
advection = 0.001

evolFun = function(X) evolAdvDiff(X, adv = advection, diff = diffusion)
max.iter = 1


## covariance parameters
sig2 = 0.0; range = 10; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)

## likelihood settings
me.var = 0.25;
data.model = "gauss"

lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha = 2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 
inds = intersect(which(locs[,1] < 0.6 & locs[,1] > 0.4 ), which(locs[,2] < 0.6 & locs[,2] > 0.4))
x0 = matrix(rep(0, dim(locs)[1]), ncol = 1)
x0[inds] = 1
#x0 = RandomFields::RFsimulate(model = RMmatern(nu = smooth, scale = range),
#                              x = locs[,1], y = locs[,2])[[1]]



RRMSPE = list(); LogSc = list()


for (iter in 1:max.iter) {  

    XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, b = 1, smooth = smooth, range = range, locs = locs)

}

########## plot results ########## 
m = M = 0
for( t in 1:Tmax ){
  zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  #zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  m = min(m, zrange[1])
  M = max(M, zrange[2])
}
for(t in (Tmax-1):Tmax){
  defpar = par(mfrow=c(1, 2), oma=c(0, 0, 2, 0))
  nna.obs = which(!is.na(XY$y[[t]]))
  fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim=c(m, M), nx=sqrt(n), ny=sqrt(n), main="obs" )
  fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim=c(m, M), nx=sqrt(n), ny=sqrt(n), main="truth" )
  #fields::quilt.plot( locs, predsVL[[t]]$state, zlim=zrange, nx=sqrt(n), ny=sqrt(n), main="prediction" )
  par(defpar)
}