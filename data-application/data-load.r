suppressPackageStartupMessages(library(tidyverse))

start = 2
end = start + 7

data.dir = '~/vecchiaFilter/data-application/data'
results.file = paste("~/vecchiaFilter/data-application/param-estimation-results", sep="")

full.data = readr::read_csv(sprintf("%s/TPW.csv", data.dir)) %>%
    filter( x < -21.76952 ) #%>%
    #filter( x < quantile(x, 0.1), y < quantile(y, 0.1) )
locations = full.data[,c("x", "y")]
nx = length(unique(locations[,1] %>% as_vector()))
ny = length(unique(locations[,2] %>% as_vector()))


seed = 1998
set.seed(seed)
 

#zlim = c(min(full.data %>% summarize( across('2':'30', function(t) min(log(t), na.rm=TRUE))) %>% as_vector()),
#         max(full.data %>% summarize( across('2':'30', function(t) max(log(t), na.rm=TRUE))) %>% as_vector()))


oldpar = par(mfrow=c(2, 4), oma=c(1, 1, 0, 0) + 1, mar=c(0, 0, 1, 1) + 1)
for( i in c(start:end) ){
    data = full.data %>% select(as.character(i)) %>% pull()
    data = -log(data+1e-5)+10
    if( all(is.na(data)) ) next()
    
    shape = mean(data, na.rm=TRUE)**2/var(data, na.rm=TRUE)
    scale = var(data, na.rm=TRUE)/mean(data, na.rm=TRUE)

    mu = mean(data, na.rm=TRUE)
    sigma = sd(data, na.rm=TRUE)
    
    xs = seq(0, max(data[is.finite(data)], na.rm=TRUE), length.out=100)
    dg = dgamma(xs, shape=shape, scale=scale)
    dn = dnorm(xs, mu, sigma)
    
    hist(data, main=sprintf("histogram of data at t=%d", i), freq=FALSE, ylim=c(0, max(dg, dn)))
    lines(xs, dg, type="l", col="red")
    lines(xs, dn, type="l", col="blue")
    
    #fields::quilt.plot( locations, log(data), zlim = zlim, nx = nx, ny = ny, main = sprintf("observations on Dec %d",i), axes=FALSE )
}
    






