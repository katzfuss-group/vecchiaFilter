### Pilot data
ranges = c(6.92, 6.898, 6.898, 5.873, 6.899)
as = c(0.1249, 0.1246, 0.1249, 0.1249, 0.1253)
sigs = c(27.12, 21.54, 21.54, 39.81, 21.54)
params = matrix(c(ranges, as, sigs), byrow=FALSE, ncol=3)
colnames(params) = c("range", "a", "sig2")


dist.params = list(log.a = list(mean = mean(log(params[,"a"])), ## we assue that log(a) is normal
                                sd = diff(range(log(params[,"a"])))),

                   ## c is normal
                   c = list(mean = 0.8, sd = 0.3),

                   ## we assume that inverse of sigma is gamma
                   inv.sig2 = list(shape = mean(1/params[,"sig2"])**2/var(1/params[,"sig2"]),
                                   rate = mean(1/params[,"sig2"])/var(1/params[,"sig2"])),

                   ## we assume that log(range) is normal
                   log.range = list(mean = mean(log(params[,"range"])),
                                    sd = diff(range(log(params[,"range"])))),

                   ## we assume that nu is constant at 1.5
                   log.nu = list(mean = log(1.5), sd = 0))
                   


getNewParticles = function( particles, d.params ) {

    Np = nrow(particles)
    
    a = exp( rnorm(Np, mean=log(particles[,"a"]), sd=d.params[["log.a"]][["sd"]] ) )

    c = rnorm( Np, mean=particles[,"c"], sd=d.params[["c"]][["sd"]] )

    # we use the particles as means and keep original variances
    or.mean = d.params[["inv.sig2"]][["shape"]] / d.params[["inv.sig2"]][["rate"]]
    or.var  = 1/d.params[["inv.sig2"]][["rate"]] * or.mean
    shape = (1/particles**2) / or.var
    rate = (1/particles) / or.var
    sig2 = 1 / rgamma( Np, shape = shape, rate = rate )
    
    range = exp( rnorm(Np, mean=log(particles[,"range"]), sd=d.params[["log.range"]][["sd"]]) )

    nu = exp( rnorm( particles[,"nu"], d.params[["log.nu"]][["mean"]], d.params[["log.nu"]][["sd"]] ) )

    new.parts = matrix( c(a, c, sig2, range, nu), ncol=5, byrow=FALSE )
    colnames( new.parts ) = c("a", "c", "sig2", "range", "nu")
    
    return( new.parts )
    
}



particles = matrix( rep(c(0.124, 0.8, 25, 6.5, 1.5), 10), ncol=5, byrow=TRUE )
colnames(particles) = c("a", "c", "sig2", "range", "nu")

new.particles = getNewParticles( particles, dist.params )
print(new.particles)
