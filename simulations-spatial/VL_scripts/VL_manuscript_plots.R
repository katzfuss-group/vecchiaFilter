require(GPvecchia)
require(ggplot2)
require(directlabels)
require(reshape)
require(scales) # transparent plots
setwd("~/HVLF")
rm(list=ls())

#######################################################################################
#####  Contains all plots that do not involve HMC/MCMC, LGCP, or satellite data  ######
#######################################################################################

map_mod = function(x) c("Gaussian", "Logistic", "Poisson","Gamma")[x]

create_mod_factor = function(agg_df){
  agg_df[,1] = map_mod(agg_df[,1])
  agg_df$Mod<- factor(agg_df$Mod, levels =c("Gaussian", "Logistic", "Poisson","Gamma") )
  return(agg_df)
}




####### MSE Plot #######
create_MSE_plot=function(data_df, dim, colScale, shapeScale, exclude_gauss=FALSE, exclude_rf = FALSE){

  # aggregate data
  agg_mse_2d = aggregate( cbind(MSE_Laplace, MSE_VL, MSE_VL_z, MSE_LowRank) ~ Mod+Neighbors+C_Smoothness,    data = data_df, mean )
  colnames(agg_mse_2d)[5:6] = c("MSE_VL-IW", "MSE_VL-RF")
  if(exclude_rf) agg_mse_2d = agg_mse_2d[,-6]
  agg_mse = agg_mse_2d

  if(ncol(agg_mse)>6)  agg_mse[,7] = sqrt(agg_mse[,7]/agg_mse[,4])
  agg_mse[,5] = sqrt(agg_mse[,5]/agg_mse[,4])
  agg_mse[,6] = sqrt(agg_mse[,6]/agg_mse[,4])
  agg_mse[,4] = agg_mse[,4]/agg_mse[,4]

  agg_mse = create_mod_factor(agg_mse)
  MSE_melted = melt(agg_mse, id = c("Mod", "Neighbors", "C_Smoothness"))
  colnames(MSE_melted)[4:5] <- c("Algorithm", "RRMSE")
  colnames(MSE_melted)[2] <- c("m")
  MSE_melted$Algorithm <- substr(MSE_melted$Algorithm,5, 99)
  
  p1 = ggplot(MSE_melted, aes(x = m, y = RRMSE, color = Algorithm, shape =Algorithm, linetype = Algorithm)) +
        geom_line() + 
        geom_point() +
        theme_bw() + 
        theme(legend.position = "top") + 
        scale_shape(solid=FALSE) +
        #facet_grid(Mod, scales = "free_y"),labeller = label_bquote(nu == .(C_Smoothness)))# + 
        colScale + 
        shapeScale

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


data_df = read.csv("VL_scripts/saved_data/simulations-marcin")
colScale_2D <- scale_colour_manual(name = "Algorithm",values = c(1,2,4,3))
shapeScale_2D <- scale_shape_manual(name = "Algorithm",values = c(NA,0,4,2))
create_MSE_plot(data_df, dim=2, colScale=colScale_2D, shapeScale=shapeScale_2D, exclude_rf = TRUE, exclude_gauss = FALSE )
