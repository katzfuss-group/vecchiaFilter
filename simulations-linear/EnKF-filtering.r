source('covFunts.r')
source('getMatCov.r')
library(Matrix)



saveEnKFResults = function(ensemble, Tmat, saveUQ){

    results = list(state = apply(ensemble, 1, mean))
    
    if(saveUQ == "W"){
        results[["W"]] = solve( Matrix(Tmat * tcrossprod( ensemble, ensemble ), sparse=TRUE), sparse=TRUE ) 
    } else if(saveUQ == "L") {
        L.tt = Cholesky( Matrix(Tmat * tcrossprod( ensemble, ensemble ), sparse=TRUE) )
        results[["L"]] = L.tt
    }
    return( results )
}



########## filtering ##########
filterEnKF = function(Nens, XY, saveUQ = ""){
   
    preds = list()
    n = length(XY$y[[1]])

    Q = GPvecchia::MaternFun(fields::rdist(locs), c(sig2, range, smooth))
    Q.chol = t(chol(Q))

    ensemble = Q.chol %*% matrix(rnorm(Nens*n), ncol=Nens)
    
    Tmat = KanterCovFun(locs, Nens)
    
    for (t in 1:Tmax) {

        cat(sprintf("EnK filtering, time t=%d\n", t))
        w = Q.chol %*% matrix(rnorm(Nens*n), ncol=Nens)
        fensemble = apply(ensemble, 2, evolFun) + w
        fcovmat = Matrix(Tmat * tcrossprod( fensemble, fensemble ), sparse=TRUE)
        
        obs.inds = which(!is.na(XY$y[[t]]))
        y = XY$y[[t]][obs.inds]
        n.obs = length(y)
        
        H = Diagonal(n)[obs.inds,]
        R = Diagonal(n.obs)*me.var
        K = fcovmat %*% t(H) %*% solve( H %*% fcovmat %*% t(H) + R )

        v = Matrix(me.var * rnorm(n.obs * Nens), ncol=Nens)
        y.tilde = apply(v, 2, '+', y)
        ensemble = fensemble + K %*% ( y.tilde - H %*% fensemble)

        preds[[t]] = saveEnKFResults(ensemble, Tmat, saveUQ)
      
    }
    return( preds )   
}
