#!/usr/bin/env Rscript
library(GPvecchia)
source("~/vecchiaFilter/simulations-spatial/run_scenario.R")

##########################################################################
##################### Compare models:  MSE and Log score  ########################
##########################################################################

## likelihood settings
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
    data.model = 'gauss'
} else if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
    stop("One of the models has to be passed as argument")
} else {
    data.model = args[1]
}

model_no = which(data.model==c("gauss", "logistic", "poisson", "gamma"))
filename = paste(data.model, "results", sep=".")

########  Setup Simulations




niter = 110
d_vals = c(1)  # domain, [0,1]
s_vals = c(2500)
seed_vals = 1:niter
smoothness_vals = .5
range_vals = .15
nbrs = c(1, 3, 5, 10, 20, 40, 60, 80, 100)
dimen = 2 
models_tested = model_no
run_laplace = TRUE



scenario_table = c()

#2d Scenarios
#nbrs = c(1, 2, 3, 4, 5, 10, 20, 50, 100, 250, 500)
#dimen = 2


# Compare log likelihood of true y for posterior density of approx p(y | z)
for (seed_r in seed_vals) {
  for (domn in d_vals) {
    for (samp_size in s_vals) {
      for (neighbors in nbrs ) {
        for (smth in smoothness_vals) {
          for (rnge in range_vals) {
              for (mod_type in models_tested) {
                  run_laplace = neighbors==1
                  scenario_table = rbind(scenario_table, c(seed_r, domn, dimen, samp_size, neighbors, smth, mod_type, rnge, TRUE, run_laplace))
            }
          }
        }
      }
    }
  }
}

aggregated_data = c()



header = c("Mod", "Domain", "Dimen", "Sample", "C_Smoothness", "C_Range","Neighbors","Seed_off",
           "MSE_Laplace", "MSE_VL", "MSE_VL_z",  "MSE_LowRank",
           "LS_Laplace", "LS_VL",  "LS_VL_z", "LS_LowRank",
           "Time_Laplace", "Time_VL", "Time_VL_z",  "Time_LowRank",
           "Iter_Laplace", "Iter_VL", "Iter_VL_z",  "Iter_LowRank")

scen_params = c("seed_r", "domn", "dimen", "samp_size", "neighbors", "smth", "mod_type", "rnge", "show_output", "run_algo")



t_start = Sys.time()



scenario_runner = create_scenario_tester(header, filename)# filename= "delete_me.csv"
for (i in 1:length(scenario_table[,1])) {
    params = setNames( as.list(scenario_table[i,]), c(scen_params))
    aggregated_data = rbind(aggregated_data, do.call(scenario_runner, params))
}


t_end = Sys.time()

colnames(aggregated_data) = header
write.csv(aggregated_data, file = filename, row.names = F)
