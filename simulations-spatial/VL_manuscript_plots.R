require(GPvecchia)
require(ggplot2)
require(directlabels)
require(reshape)
library(tidyverse)
require(scales) # transparent plots
setwd("~/HVLF")
rm(list = ls())

#######################################################################################
#####  Contains all plots that do not involve HMC/MCMC, LGCP, or satellite data  ######
#######################################################################################

map_mod = function(x) c("Gaussian", "Logistic", "Poisson","Gamma")[x]

create_mod_factor = function(agg_df){
  agg_df[,1] = map_mod(agg_df[,1])
  agg_df$Mod = factor(agg_df$Mod, levels = c("Gaussian", "Logistic", "Poisson","Gamma") )
  return(agg_df)
}




####### MSE Plot #######
create_MSE_plot = function(data_df, dim, colScale, shapeScale, exclude_rf = FALSE){

  # aggregate data
  agg_mse_2d = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_VL_z, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness, data = data_df, mean )
  colnames(agg_mse_2d)[5:7] = c("MSE_Hierachical Vecchia", "MSE_VL-RF", "MSE_Low Rank")
  if (exclude_rf) agg_mse_2d = agg_mse_2d[,-6]
  agg_mse = agg_mse_2d

  if ( ncol(agg_mse) > 6 )  agg_mse[,7] = sqrt(agg_mse[,7]/agg_mse[,4])
  agg_mse[,5] = sqrt(agg_mse[,5]/agg_mse[,4])
  agg_mse[,6] = sqrt(agg_mse[,6]/agg_mse[,4])
  agg_mse[,4] = agg_mse[,4]/agg_mse[,4]

  agg_mse = create_mod_factor(agg_mse)
  MSE_melted = melt(agg_mse, id = c("Mod", "Neighbors", "C_Smoothness"))
  colnames(MSE_melted)[4:5] <- c("Algorithm", "RRMSE")
  colnames(MSE_melted)[2] <- c("m")
  MSE_melted$Algorithm <- substr(MSE_melted$Algorithm,5, 99)
  
  p1 = ggplot(MSE_melted, aes(x = m, y = RRMSE, color = Algorithm, shape = Algorithm, linetype = Algorithm)) +
        geom_line() + 
        geom_point() +
        theme_bw() + 
        theme(legend.position = "top") + 
        scale_shape(solid = FALSE) +
        facet_grid(cols = vars(Mod)) + #,labeller = label_bquote(nu == .(C_Smoothness))) + 
        colScale# + 
  #      shapeScale

  # if(exclude_gauss){
  #   MSE_melted = MSE_melted[which(MSE_melted$Mod!="Gaussian"),]
  #   p1 = ggplot(MSE_melted, aes(x = m, y = RRMSE, color = Algorithm, shape =Algorithm, linetype = Algorithm))+
  #     geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
  #     facet_grid(Mod, scales = "free_y",labeller = label_bquote(nu == .(C_Smoothness))) +
  #     colScale+shapeScale
  #   #facet_wrap(.~C_Smoothness, scales = "free_y",labeller = label_bquote(nu == .(C_Smoothness))) +colScale+shapeScale
  # }

  p1
}




create_LS_plot = function(data_df, dim, colScale, shapeScale, exclude_rf=FALSE){
  
  # aggregate and negate logscores (for consistence appearance wrt MSE)
  agg_ls = aggregate( cbind(LS_Laplace, LS_VL, LS_VL_z, LS_LowRank) ~ Mod+Neighbors+C_Smoothness, data = data_df, mean )
  agg_ls[,4:7] = -agg_ls[,4:7]
  
  colnames(agg_ls)[5:7] = c("LS_Hierachical Vecchia", "LS_VL-RF", "LS_Low Rank")
  agg_ls = agg_ls[,-6]
  
  if (ncol(agg_ls) > 6) agg_ls[,7] = (agg_ls[,7] - agg_ls[,4])
  agg_ls[,6] = (agg_ls[,6] - agg_ls[,4])
  agg_ls[,5] = (agg_ls[,5] - agg_ls[,4])
  agg_ls[,4] = agg_ls[,4] - agg_ls[,4]
  
  agg_ls = create_mod_factor(agg_ls)
  LS_melted = melt(agg_ls, id = c("Mod", "Neighbors", "C_Smoothness"))
  colnames(LS_melted)[4:5] <- c("Algorithm", "dLS")
  colnames(LS_melted)[2] <- c("m")
  LS_melted$Algorithm <- substr(LS_melted$Algorithm,4,99)
  
  p1 = ggplot(LS_melted, aes(x = m, y = dLS, color = Algorithm, shape = Algorithm, linetype = Algorithm)) +
        geom_line() + 
        geom_point() + 
        theme_bw() + 
        theme(legend.position = "top") +
        scale_shape(solid = FALSE) +
        facet_grid(C_Smoothness~Mod, labeller = label_bquote(nu == .(C_Smoothness))) + 
        colScale + 
        shapeScale
  
  # if (FALSE) {
  #   colnames(LS_melted)[2] <- c("m")
  #   LS_melted = LS_melted[which(LS_melted$Mod != "Gaussian"),]
  #   LS_melted = LS_melted[which(LS_melted$Algorithm != "VL_I"),]
  #   ggplot(LS_melted, aes(x = m, y = LS, color = Algorithm, shape = Algorithm, linetype = Algorithm)) +
  #     geom_line() + geom_point() + theme_bw() + theme(legend.position = "bottom") + scale_shape(solid = FALSE) +
  #     facet_wrap(.~C_Smoothness, scales = "free_y", labeller = label_bquote(nu == .(C_Smoothness))) + colScale + shapeScale
  #   #ggsave("LGCP_LS_pois_v4.pdf", device= "pdf",width = 7, height = 4)
  #   
  #   pois_subset = LS_melted[which(LS_melted$Mod == "Poisson" & LS_melted$C_Smoothness == .5) ,]
  #   #pois_subset[which(pois_subset$Algorithm=="VL_ZY"),][,4] = "LOCA"
  #   pois_subset = pois_subset[which(pois_subset$Algorithm != "VL_I"),]
  #   
  #   
  #   ptit = expression("2D Poisson Samples, n=2500,"~nu == .5)
  #   ggplot(pois_subset, aes(x = Neighbors, y = LS, color = Algorithm, shape = Algorithm, linetype = Algorithm)) +
  #     geom_line() + geom_point() + theme_bw() + theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5)) +
  #     scale_shape(solid = FALSE) + colScale
  #   #ggsave("2D_poisson_alone.pdf", device= "pdf",width = 4, height = 4)
  # }
  p1
  
}



data_df = read_csv("simulations-spatial/results") 

MSE_Laplace_filled = data_df %>% select(MSE_Laplace) %>% as_vector("numeric")
LS_Laplace_filled = data_df %>% select(LS_Laplace) %>% as_vector("numeric")
cursor_mse = NULL
cursor_ls = NULL
for (i in 1:length(MSE_Laplace_filled)) {
  if (MSE_Laplace_filled[i] > 0) {
    cursor_mse = MSE_Laplace_filled[i]
  } else {
    MSE_Laplace_filled[i] = cursor_mse
  }
  if (LS_Laplace_filled[i] != 0) {
    cursor_ls = LS_Laplace_filled[i]
  } else {
    LS_Laplace_filled[i] = cursor_ls
  }
}

data_df = data_df %>%
          mutate( MSE_Laplace = MSE_Laplace_filled) %>% 
          mutate( LS_Laplace = LS_Laplace_filled) %>%
          filter(Seed_off < 11) %>%
          select(-Domain, -Dimen, -Sample, -C_Range, -Time_Laplace,
                 -Time_VL, -Time_VL_z, -Time_LowRank, -Iter_Laplace,
                 -Iter_VL, -Iter_LowRank, -Iter_VL_z)

colScale_2D <- scale_colour_manual(name = "Algorithm",values = c(1,2,4,3))
shapeScale_2D <- scale_shape_manual(name = "Algorithm",values = c(NA,0,4,2))
create_MSE_plot(data_df, dim = 2, colScale = colScale_2D, shapeScale = shapeScale_2D, exclude_rf = TRUE)
create_LS_plot(data_df, dim = 2, colScale = colScale_2D, shapeScale = shapeScale_2D, exclude_rf = TRUE)
