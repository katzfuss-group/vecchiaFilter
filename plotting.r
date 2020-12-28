plotResults2d = function(XY, preds1, preds2, data.model, resultsDir, name1, name2){

    Tmax = lenght(XY$x)
    m = M = 0
    for (t in 1:Tmax) {
        zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE)))))
        m = min(m, zrange[1])
        M = max(M, zrange[2])
    }
    for (t in 1:Tmax) {
        
        if(t<10){
            number = paste("0", t, sep="")  
        } else {
                number = t
        }
        pdf(paste(resultsDir, "/", data.model, "/", number, ".pdf",sep=""), width=8, height=8)
        defpar = par(mfrow = c(2, 2), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
        nna.obs = which(!is.na(XY$y[[t]]))
        fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], nx = sqrt(n), ny = sqrt(n), main = "obs" )
        fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "truth" )
        fields::quilt.plot( locs, preds1[[t]]$state, zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = paste("prediction", name1) )
        fields::quilt.plot( locs, preds2[[t]]$state, zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = paste("prediction", name2) )
        par(defpar)
        dev.off()
    }

}

