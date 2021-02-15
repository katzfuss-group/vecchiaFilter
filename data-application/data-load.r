suppressPackageStartupMessages(library(tidyverse))

estimate = function( data, locs, m, theta.ini, reltol, conditioning ){


    ## specify vecchia approximation
    vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning=conditioning)


    ## initial covariance parameter values

     if(missing(theta.ini) || is.na(theta.ini)){         
      stop("Initial cov. parameter values must be specified")
     }
  

    ## specify vecchia loglikelihood
    n.par=length(theta.ini)

    negloglik.vecchia=function(logparms){
        if(exp(logparms[3])>10 && is.character(covmodel) && covmodel=='matern'){
            stop("The default optimization routine to find parameters did not converge. Try writing your own optimization.")
        }
        GPvecchia::vecchia_laplace_likelihood( data,
                                              vecchia.approx,
                                              'gamma',
                                              exp(logparms[-n.par]),
                                              covmodel = 'matern',
                                              likparms = list(shape=exp(logparms[n.par])) )
    }

    ## find MLE of theta (given beta.hat)
    output.level = 9
    opt.result=stats::optim(par=log(theta.ini),fn=negloglik.vecchia,
                            control=list(
                                trace=output.level,maxit=300, parscale=log(theta.ini),
                                reltol=reltol
                            )) # trace=1 outputs iteration counts
    theta.hat=exp(opt.result$par)
    names(theta.hat) = c("variance", "range", "smoothness", "nugget")

    return( theta.hat )
    
}





day = 2
start = lubridate::with_tz(lubridate::now(), "America/Chicago")

data.dir = '~/vecchiaFilter/data-application/data'
results.file = paste("~/vecchiaFilter/data-application/param-estimation-results", sep="")

full.data = readr::read_csv(sprintf("%s/TPW.csv", data.dir), n_max = 10000)
locations = full.data[,c("x", "y")]

seed = 1998
set.seed(seed)

m = 50
sample.ns = c(5000, 10000)
tolerances = c(1e-4, 1e-5)

theta.ini = c(1e6, 0.15, 1.5, 200)




for(reltol in tolerances){
    for(n.sample in sample.ns){

        cat(sprintf("Running the optimization for m=%d, n.sample=%d and rel.tol=%f\n", m, n.sample, reltol))
        
        data = full.data %>%
            sample_n(n.sample) %>%
            select(c(x, y, !!as.character(day))) %>%
            rename(values=!!as.character(day)) %>%
            filter(!is.na(values)) %>% mutate(values = values - min(values))
        locations = as.matrix(select(data, x, y))

        if (!exists("theta.ini")){
            theta.ini = as.numeric(rep(NA, 4))
        }
        t0 = proc.time()
            
        vecchia.est = estimate(data$values,
                               locations,
                               m=m,
                               theta.ini = theta.ini,
                                        #cond.yz = 'y',
                               conditioning = 'mra',
                               reltol = reltol)
        elapsed = as.numeric(proc.time()[3] - t0[3])
        cat(sprintf("Estimation took %.2f seconds\n", elapsed))
        end = lubridate::with_tz(lubridate::now(), Sys.timezone())
            
        this.result = tibble(
            set = day,
            nsample = n.sample,
            m = m,
            elapsed = elapsed,
            start = start,
            stop = end,
            var_ini = theta.ini[1],
            range_ini = theta.ini[2],
            smooth_ini = theta.ini[3],
            nugget_ini = theta.ini[4],
            var_est = vecchia.est$theta.hat[1],
            range_est = vecchia.est$theta.hat[2],
            smooth_est = vecchia.est$theta.hat[3],
            nugget_est = vecchia.est$theta.hat[4],
            ic0 = TRUE,
                                        #method = "MRA",
            method = 'MRA',
            seed = seed,
            reltol = reltol)
        
        
        results = read_csv(results.file, col_types=cols())
        results = rbind(results, this.result)
        write_csv(results, results.file)
        theta.ini = vecchia.est$theta.hat
    }
    rm("theta.ini")
}


