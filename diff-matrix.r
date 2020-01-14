function = diffAdvMat(Nx, Ny, a=1, diff=0.5, advection=0.0, scheme='forward'){
  
  x = seq(from=0, to=1, length.out=Nx)
  y = seq(from=0, to=1, length.out=Ny)
  
  dx = x[2] - x[1]
  dy = y[2] - y[1]
  
  # Mesh Fourier numbers in each direction
  diff_x = a*float(diff)/(dx**2)
  diff_y = a*float(diff)/(dy**2)
  
  Ix = seq(Nx+2)
  Iy = seq(Ny+2)
  
  # Data structures for the linear system
  N = (Nx+1)*(Ny+1)
  E = matrix(zeros(N, N)
  
  
  m = function(i, j){ j*(Nx+1)+i)}
  
  for(j in Iy){
    for(i in Ix){
      p = m(i,j)
      if(j>0){
        stop("function not implemented yet") 
      }
    }
  }
  
}

function = Adv1d(Nx, adv=0.0){
  
  E = matrix(rep(0, (Nx+2)**2), ncol=(Nx+2))
  diag(E)=-1
  
}