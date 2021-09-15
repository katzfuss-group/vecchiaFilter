#### Covparm and shape estimation ####
update_a = function(a_init, covparms,  vecchia.approx, vecchia.approx.IW, XB, obs){
  cat("\tEstimating a\n")
  a_prev = a_init

  for(i in 1:10){
    t_start = Sys.time()
    posterior = GPvecchia::calculate_posterior_VL(obs,
                                       vecchia.approx,
                                       likelihood_model="gamma",
                                       covparms=covparms,
                                       likparms = list("alpha"=a_prev),
                                       prior_mean = XB)
    mu = posterior$mean + XB
    llh = function(a) -sum( -a * exp(-mu) * obs + (a-1) * log(obs) + a * log(a) - a * mu - log(gamma(a))) # concave in a and XB
    param_est = optim(a_prev, llh, method = "Brent", lower = .01, upper = 1e2)
    a = param_est$par
    if(abs(a-a_prev) < 1e-5) {
        cat(sprintf("\tconvergence criteria met (fitting shape parameter)\n"))
        cat(sprintf("\ta = %f\n", a))
        break
    }
    a_prev = a
  }
  return(a)
}




fit_me_var = function(me_var_init, vecchia_approx, XB, obs, covparms, max_iter = 20) {
    cat("\tEstimating m.e. var.:\n")
 
    me_var = 0
    me_var_prev = me_var_init
    n = nrow(XB)

    iter = 0
    while (iter < max_iter) {
        posterior = GPvecchia::calculate_posterior_VL(obs,
                                                      vecchia_approx,
                                                      likelihood_model = "gauss",
                                                      covparms = covparms,
                                                      likparms = list(sigma = sqrt(me_var_prev)),
                                                      prior_mean = XB)
        mu = posterior$mean
        y = matrix(obs)
        me_var = mean((obs - mu)**2)
        if (abs(me_var - me_var_prev) < 1e-5)
            break
        me_var_prev = me_var
        iter = iter + 1
    }

    if (iter < maxiter) {
        cat(sprintf("\t\tconvergence reached after %d iterations\n", iter))
    } else {
        cat(sprintf("\t\tno convergence after %d iterations\n", iter))
    }
    cat(sprintf("\t\testimated parameter: %f\n", me_var))
    return (me_var)
}




fit_covparms = function(a, covparms_init, vecchia.approx, vecchia.approx.IW, XB, obs){
    cat("\tEstimating covariance parameters:\n")
    
    vl_likelihood = function(x0){
        covparms = c(exp(x0), SMOOTH) ## sigma range smoothness
        cat(sprintf("\t\tEvaluating covparms = (%.4f %.4f %.4f)\n", covparms[1], covparms[2],covparms[3]))
        default_lh_params = list("alpha" = a, "sigma" = sqrt(.1))
        ## Perform inference on latent mean with Vecchia Laplace approximation

        vll = GPvecchia::vecchia_laplace_likelihood(obs,
                                                    vecchia.approx,
                                                    likelihood_model="gamma",
                                                    covparms = covparms,
                                                    return_all = FALSE,
                                                    likparms = default_lh_params,
                                                    prior_mean = as.numeric(XB))#,
                                                    #vecchia.approx.IW=vecchia.approx.IW,
                                                    #y_init = NA )
        cat(sprintf("\t\tLikelihood for covparms = (%.4f %.4f %.4f): %.4f\n",
                    covparms[1], covparms[2],covparms[3], vll))
        return(-vll)
    }

    x0 = log(covparms_init)
    vl = vl_likelihood(x0)

    #cat(sprintf("\t\tInitial likelihood %.6f at values (%f, %f, %f)\n", vl, exp(x0[1]), exp(x0[2]), exp(x0[3])))
    cat(sprintf("\t\tInitial likelihood %.6f at values (%f, %f)\n", vl, exp(x0[1]), exp(x0[2])))#, exp(x0[3])))
    if( vl==Inf || vl==-Inf )
        stop()
    
    res = optim(x0,
                vl_likelihood,
                method = "Nelder-Mead",
                control = list("trace" = 0, "maxit" = 500, "reltol" = 1e-5, parscale = x0))

    cat(sprintf("\tconvergence result: %s\n", res$convergence))
    cat(sprintf("\testimated parameter vector: %f, %f, %f\n", exp(res$par)[1], exp(res$par)[2], exp(res$par)[3]))

    return(exp(res$par))
}



## this is a newer and more general version of the function fit_covparms allowing for different
## likelihood families
fit_covparams = function(covparms_init, vecchia_approx, XB, obs, lik_params){

    cat("\tEstimating covariance parameters:\n")
    
    vl_likelihood = function(x0){
        covparms = c(exp(x0), SMOOTH) ## sigma range smoothness
        cat(sprintf("\t\tEvaluating covparms = (%.4f %.4f %.4f)\n", covparms[1], covparms[2], covparms[3]))

        ## Perform inference on latent mean with Vecchia Laplace approximation
        vll = GPvecchia::vecchia_laplace_likelihood(obs,
                                                    vecchia_approx,
                                                    likelihood_model = lik_params[["data.model"]],
                                                    covparms = covparms,
                                                    return_all = FALSE,
                                                    likparms = lik_params,
                                                    prior_mean = as.numeric(XB))
        cat(sprintf("\t\tLikelihood for covparms = (%.4f %.4f %.4f): %.4f\n",
                    covparms[1], covparms[2], covparms[3], vll))
        return(-vll)
    }

    x0 = log(covparms_init)
    vl = vl_likelihood(x0)

    #if( vl==Inf || vl==-Inf )
    #    stop()

    res = optim(x0,
                vl_likelihood,
                method = "Nelder-Mead",
                control = list("trace" = 0, "maxit" = 500, "reltol" = 1e-5, parscale = x0))

    cat(sprintf("\tconvergence result: %s\n", res$convergence))
    cat(sprintf("\testimated parameter vector: %f, %f, %f\n", exp(res$par)[1], exp(res$par)[2], exp(res$par)[3]))

    return(exp(res$par))
}
