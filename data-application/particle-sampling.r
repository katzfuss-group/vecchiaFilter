### Pilot data
ranges = c(6.92, 6.898, 6.898, 5.873, 6.899)
as = c(0.1249, 0.1246, 0.1249, 0.1249, 0.1253)
sigs = c(27.12, 21.54, 21.54, 39.81, 21.54)
params = matrix(c(ranges, as, sigs), byrow=FALSE, ncol=3)
colnames(params) = c("range", "a", "sig2")


prior    = list(a = list(mean = log(2), sd = 0),
                c = list(mean = 1.0, sd = 0),
                sig2 = list(mean = log(1.0), sd = 0),
                range = list(mean = log(0.15), sd = 0.01),
                nu = list(mean = log(0.5), sd = 0))


prop     = list(a = list(sd = 0),
                c = list( sd = 0 ),
                sig2 = list( sd = 0 ),
                range = list( sd = 0.01 ),
                nu = list(sd = 0))



log.dist.eval = function( name, particles, distributions ){

    value = particles[[name]]
    distribution = distributions[[name]]

    if( distribution[["sd"]]==0 ){
        1e8
    } else if( name %in% c("a", "sig2", "range", "nu") ) {
        sum(dnorm( log(value), distribution$mean, distribution$sd, log=TRUE ))
    } else if( name=="c" ) {
        sum(dnorm( value, distribution$c$mean, distribution$c$sd, log=TRUE ))
    }
    
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
    range  = c(0.148, 0.15)#exp( rnorm( Np, distr[["range"]][["mean"]], distr[["range"]][["sd"]] ) )
    nu     = exp( rnorm( Np, distr[["nu"]][["mean"]], distr[["nu"]][["sd"]] ) )
    
    new.parts = matrix( c(a, c, sig2, range, nu), ncol=5, byrow=FALSE )
    colnames( new.parts ) = c("a", "c", "sig2", "range", "nu")

    return( new.parts )
    
}


#cat("Test: 10 sample particles \n")
#print(sample.particles( 10, prior ))
