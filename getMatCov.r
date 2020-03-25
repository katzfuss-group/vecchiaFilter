getMatCov = function(V, covariances, factor=FALSE){

  if (factor){
    if (methods::is(covariances, "sparseMatrix")) {
      L = methods::as(covariances, "CsparseMatrix")[V$ord, V$ord]
      M = getMatCovFromFactorCpp(L, V$U.prep$revNNarray)
      M[ M==0 ] = NA
      return(M)
    } else {
      getMatCovFromFactor(V, covariances) 
    }
  } else if (class(covariances)=='function') {
    getMatCovFromFunction(V, covariances)
  } else if (class(covariances)=='matrix') {
    getMatCovFromMat(V, covariances)
  } else {
    stop("Wrong covariance format passed")
  }
}


getMatCovFromFunction = function(V, covfun){
  revNNarray = V$U.prep$revNNarray
  rows = c()
  cols = c()
  for(i in 1:nrow(revNNarray)){
    r = revNNarray[i,];
    newrows = rep(i, sum(!is.na(r)))
    newcols = r[!is.na(r)]
    rows = c(rows, newrows)
    cols = c(cols, newcols)
  }
  n = nrow(V$locsord)
  
  Sig.sel = matrix(rep(NA, length(revNNarray)), nrow=ncol(revNNarray))
  inds_to_fill=which(!is.na(t(revNNarray)))

  d = matrix(sqrt(rowSums((V$locsord[rows,] - V$locsord[cols,])**2)))
  vals = as.numeric(covfun(d))

  Sig.sel[inds_to_fill] = vals
  Sig.sel = t(Sig.sel)
  return(Sig.sel)
}


getMatCovFromMat = function(V, Sigma){
  revNNarray = V$U.prep$revNNarray
  rows = c()
  cols = c()
  for(i in 1:nrow(revNNarray)){
    r = revNNarray[i,];
    newrows = rep(i, sum(!is.na(r)))
    newcols = r[!is.na(r)]
    rows = c(rows, newrows)
    cols = c(cols, newcols)
  }
  n = dim(Sigma)[1]
  inds = cbind(rows, cols)
  inds = as.vector(sapply(seq(nrow(inds)), function(r) inds[r,2]-1+n*(inds[r,1]-1)+1))
  
  Sig.sel = matrix(rep(NA, length(revNNarray)), nrow=ncol(revNNarray))
  inds_to_fill = which(!is.na(t(revNNarray)))
  Sigma.ord = Sigma[V$ord, V$ord]
  
  Sig.sel[inds_to_fill] = Sigma.ord[inds]
  Sig.sel = t(Sig.sel)
  return(Sig.sel)  
}



getMatCovFromFactor = function(V, L.org){
  revNNarray = V$U.prep$revNNarray
  
  Sig.sel = Matrix::Matrix(rep(NA, nrow(revNNarray)*ncol(revNNarray)), ncol = ncol(revNNarray))
  
  if( class(L.org)=="dgCMatrix" ){
    L = methods::as(L.org, "RsparseMatrix")[V$ord, V$ord]
  } else {
    L = L.org[V$ord, V$ord]
  }
  
  for(i in 1:nrow(revNNarray)){
    r = revNNarray[i,]
    r.inds = r[which(!is.na(r))]
    this.row = Matrix::Matrix(L[i,], nrow=1)
    nrow = length(r.inds)
    if(nrow==1){
      submatrix = Matrix::Matrix(L[r.inds,], nrow=length(r.inds))  
    } else {
      submatrix = Matrix::Matrix(L[r.inds,])  
    }
    Sig.sel[i, which(!is.na(r))] = this.row %*% Matrix::t(submatrix)
  }
  return(Sig.sel)
}
