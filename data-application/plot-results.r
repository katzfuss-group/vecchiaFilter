plotFields = function(XY, preds, locs, x0 = NULL) {

    Tmax = length(XY$x)
    N = nrow(locs)

    nx = length(unique(locs[1,]))
    ny = length(unique(locs[2,]))
    
    if (round(sqrt(N)) != sqrt(N)) {
        warning("The size of the grid is not a square of an integer. Plots might look weird")
    }

    #inds = sapply(1:TMAX, function(t) which(sapply(predsMRA$preds[[t]], function(s) length(s))>1)[1])


    zlim_x = range(sapply(1:Tmax, function(t) range(XY$x[[t]], na.rm = TRUE)), XY$x0)
    zlim_y = range(sapply(1:Tmax, function(t) range(XY$y[[t]], na.rm = TRUE)))
    #zlim_p = range(sapply(1:Tmax, function(t) range(preds$preds[[t]][[inds[t]]]$state)))
    zlim_p = range(sapply(1:Tmax, function(t) range(preds[[t]]$state)))

    if (!is.null(x0)) {
        pdf(sprintf("data-application/tests/init-field.pdf"))
        fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$x0), zlim = zlim_x, nx = sqrt(N), ny = sqrt(N), main = "init")
        dev.off()
    }

    for (t in 1:Tmax) {
        pdf(sprintf("data-application/tests/test-field-%d.pdf", t), width=18, height=6)
        oldpar = par(mfrow = c(1, 3))
        fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$x[[t]]), zlim = zlim_x, nx = sqrt(N), ny = sqrt(N), main = sprintf("truth %d", t))
        fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(XY$y[[t]]), zlim = zlim_y, nx = sqrt(N), ny = sqrt(N), main = sprintf("data %d", t))
        fields::quilt.plot(locs[, 1], locs[, 2], as.numeric(preds[[t]]$state), zlim = zlim_p, nx = sqrt(N), ny = sqrt(N), main = "reconstruction")
        par(oldpar)
        dev.off()
    }
    par(oldpar)
}



plotParamPaths = function(particles, indices, data.model, truth) {

    Tmax = length(particles)
    time = 1:Tmax
    
    ## parameters
    pdf(sprintf("data-application/tests/%s-particles-over-time.pdf", data.model), width = 6, height = 4)
    oldpar = par(mfrow = c(2, 3))

    for (t in time) {
        inds = indices[[t]]
        particles[[t]] = particles[[t]][inds, ]
    }

    
    for (par.name in colnames(particles[[1]])) {    

        uq = sapply(particles, function(t) as.numeric(quantile(t[inds, par.name], 0.95)))
        means = sapply(particles, function(t) as.numeric(mean(t[inds, par.name])))
        lq = sapply(particles, function(t) as.numeric(quantile(t[inds, par.name], 0.05)))
        
        if (all(uq == lq)) {
            cat(sprintf("%s was not sampled\n", par.name))
            next()
        }
    
        ylim = range(c(uq, lq))
        ylim[2] = 1.1 * ylim[2]
        plot(time, uq, type = "l", lty = 2, ylim = ylim, ylab = par.name, main = par.name)
        lines(time, means, type = "l")
        lines(time, lq, type = "l", lty = 2)
        abline(h = truth[[par.name]], col = "red")
    }
    par(oldpar)
    dev.off()
}




plot1dLogLikelihood = function(all_particles, logliks, truth) {

    particles = all_particles[[length(all_particles)]]
    all_params = colnames(particles)
    
    params = c()
    for (par_name in all_params) {
        unique_particles = unique(particles[, par_name])
        if (length(unique_particles)>1) {
            params = c(params, par_name)
        }
    }

    oldpar = par(mfrow = c(1, length(params)))

    for (param in params) {
        this_param = particles[, param]
        ord = order(this_param)
        loglik = logliks[[length(logliks)]][ord]
        plot(this_param[ord], loglik, type = "l", xlab = param, ylab = "loglikelihood",
             main = sprintf("likelihood at time t=%d", length(logliks)))
        abline(v = truth[[param]], col = "red")
    }
    par(oldpar)
}


plotMarginalDistribution = function(all_particles, all_indices, truth) {

    inds = all_indices[[length(all_indices)]]
    particles = all_particles[[length(all_particles)]]
    all_params = colnames(particles)

    params = c()
    for (par_name in all_params) {
        unique_particles = unique(particles[, par_name])
        if (length(unique_particles)>1) {
            params = c(params, par_name)
        }
    }

    oldpar = par(mfrow = c(1, length(params)))
    for (param in params) {
        hist(particles[inds, param], xlab = param, main = sprintf("Histogram of %s", param))
        abline(v = truth[[param]], col = "red", lwd = 2, lty = 2)
    }
    par(oldpar)
    
}
