setwd("~/vecchiaFilter")
library(Matrix)
source('aux-functions.r')
source('covFunts.r')
source('getMatCov.r')
Rcpp::sourceCpp('src/getMatCovFromFactor.cpp')
source('scores.r')


saveResults = function(preds.aux, L.tt, saveUQ){

    results = list(state = matrix(preds.aux$mean, ncol = 1))
    if(saveUQ == "W"){
        results[["W"]] = preds.aux$W
    } else if(saveUQ == "L") {
        results[["L"]] = L.tt
    }
    return( results )
}



########## filtering ##########
filterEnKF = function(Nens, XY){
   
    preds = list()
    n = length(XY$y[[1]])
    ensemble = w = v = Matrix::Matrix(rep(0, n*Nens), ncol=Nens)
    for(i in 1:Nens){
        ensemble[,i] = RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range),
                                              x = locs[,1], y = locs[,2], spConform = FALSE) 
    }
    
    Tmat = KanterCovFun(locs, Nens)
    
    for (t in 1:Tmax) {


        cat(sprintf("EnK filtering, time t=%d\n", t))
        
        for(i in 1:Nens){
          w[,i] = RandomFields::RFsimulate(model = RandomFields::RMmatern(nu = smooth, scale = range, var = sig2),
                                         x = locs[,1], y = locs[,2], spConform = FALSE) 
        }
        
        E = evolFun(Matrix::Diagonal(n))
        fensemble = E %*% ensemble + w
        
        fcovmat = Matrix::Matrix(Tmat * tcrossprod( fensemble, fensemble ), sparse=TRUE)
        
        obs.inds = which(!is.na(XY$y[[t]]))
        y = XY$y[[t]][obs.inds]
        n.obs = length(y)
        
        H = Matrix::Diagonal(n)[obs.inds,]
        R = Matrix::Diagonal(n.obs)*me.var
        K = fcovmat %*% t(H) %*% solve( H %*% fcovmat %*% t(H) + R )

        v = Matrix(me.var * rnorm(n.obs * Nens), ncol=Nens)
        y.tilde = apply(v, 2, '+', y)
        ensemble = fensemble + K %*% ( y.tilde - H %*% fensemble)
        
        preds[[ t ]] = list(state = apply(ensemble, 1, mean),
                            W = Matrix::solve( Matrix::Matrix(Tmat * tcrossprod( ensemble, ensemble ), sparse=TRUE) ), sparse=TRUE )
      
    }
    return( preds )   
}
