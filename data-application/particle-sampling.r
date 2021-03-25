### Pilot data
ranges = c(6.92, 6.898, 6.898, 5.873, 6.899)
as = c(0.1249, 0.1246, 0.1249, 0.1249, 0.1253)
sigs = c(27.12, 21.54, 21.54, 39.81, 21.54)
params = matrix(c(ranges, as, sigs), byrow=FALSE, ncol=3)
colnames(params) = c("range", "a", "sig2")


prior    = list(a     = list(mean = log(2),    sd = 0),
                c     = list(mean = 0.8,       sd = 0),
                sig2  = list(mean = log(0.36), sd = 0),
                range = list(mean = log(0.15), sd = 0.3),
                nu    = list(mean = log(0.5),  sd = 0))


prop     = list(a     = list(sd = 0),
                c     = list(sd = 0),
                sig2  = list(sd = 0),
                range = list(sd = 0.2),
                nu    = list(sd = 0))



log.dist.eval = function( name, particles, distributions, l ){

    value = particles[[name]]
    distribution = distributions[[name]]

    if( distribution[["sd"]]==0 ){
        return(1e8)
    }
    if( name %in% c("a", "sig2", "range", "nu") ) {
        value = log(value)
    }
    if( length(distribution$mean) > 1 ){
        m = distribution$mean[l]
    } else {
        m = distribution$mean[1]
    }

    if( length(distribution$sd) > 1 ){
        s = distribution$sd[l]
    } else {
        s = distribution$sd[1]
    }
    
    dnorm( value, m, s, log=TRUE )
    
}



update.sampling = function( p, prop ){
    
    prop$a$mean = log(p[,"a"])
    prop$c$mean = p[,"c"]
    prop$sig2$mean = log(p[,"sig2"])
    prop$range$mean = log(p[,"range"])
    prop$nu$mean = log(p[,"nu"])
    return( prop )
}



sample.particles = function( Np, distr ){
    
    a      = exp( rnorm( Np, distr[["a"]][["mean"]], distr[["a"]][["sd"]] ) )
    c      = rnorm( Np, distr[["c"]][["mean"]], distr[["c"]][["sd"]] )
    sig2   = exp( rnorm( Np, distr[["sig2"]][["mean"]], distr[["sig2"]][["sd"]] ) )
    range  = exp( rnorm( Np, distr[["range"]][["mean"]], distr[["range"]][["sd"]] ) )
    nu     = exp( rnorm( Np, distr[["nu"]][["mean"]], distr[["nu"]][["sd"]] ) )
    new.parts = matrix( c(a, c, sig2, range, nu), ncol=5, byrow=FALSE )    
    #np = Np**(1/2)
    #c = seq(0.5, 1.5, length.out = np)
    #sig2 = seq(10, 35, length.out = np)
    #range = seq(0.1, 3, length.out = np)
    #new.parts = expand.grid(c, sig2, range)
    #new.parts = expand.grid(sig2, range)
    #new.parts = as.matrix(cbind(2, c, sig2, range, 1.5))   
    
    colnames( new.parts ) = c("a", "c", "sig2", "range", "nu")
    return( new.parts )
    
}
#cat("Test: 10 sample particles \n")
#print(sample.particles( 10, prior ))
