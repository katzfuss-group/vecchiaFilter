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
    prop$c$mean = log(p[,"c"])
    prop$sig2$mean = log(p[,"sig2"])
    prop$range$mean = log(p[,"range"])
    prop$nu$mean = log(p[,"nu"])
    return( prop )
}



sample.particles = function( Np, distr ){
    
    a      = exp(rnorm(Np, distr[["a"]][["mean"]], distr[["a"]][["sd"]]))
    c      = exp(rnorm(Np, distr[["c"]][["mean"]], distr[["c"]][["sd"]]))
    sig2   = exp(rnorm(Np, distr[["sig2"]][["mean"]], distr[["sig2"]][["sd"]]))
    range  = exp(rnorm(Np, distr[["range"]][["mean"]], distr[["range"]][["sd"]]))
    nu     = exp(rnorm(Np, distr[["nu"]][["mean"]], distr[["nu"]][["sd"]]))
    new.parts = matrix(c(a, c, sig2, range, nu), ncol=5, byrow = FALSE)    
    
    colnames( new.parts ) = c("a", "c", "sig2", "range", "nu")
    return( new.parts )
    
}
