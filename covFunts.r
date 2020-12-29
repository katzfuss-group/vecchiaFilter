determine_radius = function(locs, N){
    
    h = abs(locs[1] - locs[2])
    s = floor(sqrt(N))
    
    if( s %% 2 == 0) {
        sf = s-1
    } else {
        sf = s    
    }
    
    if( N==sf**2 ){
        return( h*1.01*(sf-1)/2*sqrt(2) )
    }
    
    base = (sf-1)/2
    
    intervals = c(sf**2)
    while( tail(intervals, 1)<(sf+2)**2 ){
        if(length(intervals)==1 || (sf+2)**2 - tail(intervals, 1)==4){
            intervals = c(intervals, tail(intervals, 1) + 4)
        } else {
            intervals = c(intervals, tail(intervals, 1) + 8)
        }
    }
    
    ind = findInterval(N, intervals)
    middle = (intervals[ind] + intervals[ind+1])/2
    
    if( N<=middle ){
        app_ind = ind-1
    } else {
        app_ind = ind
    }
    
    if( app_ind == 0 ){
        return( h*base*sqrt(2) + h*0.01 )
    } else {
        return( h*sqrt( (base+1)**2 + (app_ind-1)**2 ) + h*0.01 )
    }
}

KanterCovFun = function(locs, condSetSize){

    radius = determine_radius(locs, condSetSize)
    D = Matrix::Matrix(fields::rdist(locs)/radius)
    piD2 = 2*pi*D
    R = (1-D) * sin(piD2)/piD2 + (1/pi) * (1-cos(piD2))/piD2
    R[D>1]=0
    R[D==0]=1
    
    return( Matrix::Matrix(R, sparse=TRUE) )
}


#rm(list=ls())

#n = 20**2
#N = 22
#nrow = 130


#grid.oneside = seq(0,1,length = round(sqrt(n)))
#locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 

#M = KanterCovFun(locs, N)
#nnz = max(apply(M, 1, function(r) length(which(r>0))))
#image(matrix(ceiling(M[nrow,]), ncol=sqrt(n)), main=paste("nnz =",nnz))
