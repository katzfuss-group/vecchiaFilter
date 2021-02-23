### Pilot data
ranges = c(6.92, 6.898, 6.898, 5.873, 6.899)
as = c(0.1249, 0.1246, 0.1249, 0.1249, 0.1253)
sigs = c(27.12, 21.54, 21.54, 39.81, 21.54)
params = matrix(c(ranges, as, sigs), byrow=FALSE, ncol=3)
colnames(params) = c("range", "a", "sig2")


prior    = list(a = list(mean = log(2), sd = 0.1),
                c = list(mean = 0.8, sd = 0.11),
                sig2 = list(mean = log(5.0), sd = 0.1),
                range = list(mean = log(7), sd = 0.1),
                smooth = list(mean = log(1.5), sd = 0))


prop     = list(a = list(sd = 0.1),
                c = list( sd = 0.1 ),
                sig2 = list( sd = 0.1 ),
                range = list( sd = 0.1 ),
                smooth = list(sd = 0))



log.dist.eval = function( name, particles, distributions, l ){

    value = particles[[name]]
    distribution = distributions[[name]]

    if( distribution[["sd"]]==0 ){
        return(1e8)
    }
    if( name %in% c("a", "sig2", "range", "smooth") ) {
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
    prop$smooth$mean = log(p[,"smooth"])

    return( prop )
}


sample.particles = function( Np, distr ){

    
    a      = exp( rnorm( Np, distr[["a"]][["mean"]], distr[["a"]][["sd"]] ) )
    c      = rnorm( Np, distr[["c"]][["mean"]], distr[["c"]][["sd"]] )
    sig2   = exp( rnorm( Np, distr[["sig2"]][["mean"]], distr[["sig2"]][["sd"]] ) )
    range  = exp( rnorm( Np, distr[["range"]][["mean"]], distr[["range"]][["sd"]] ) )
    smooth     = exp( rnorm( Np, distr[["smooth"]][["mean"]], distr[["smooth"]][["sd"]] ) )
    
    new.parts = matrix( c(a, c, sig2, range, smooth), ncol=5, byrow=FALSE )
    colnames( new.parts ) = c("a", "c", "sig2", "range", "smooth")
    
    return( new.parts )
    
}


#cat("Test: 10 sample particles \n")
#print(sample.particles( 10, prior ))
