rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results\\results\\p0=.2varre=.2')

aaa <- NULL
for( r in 5:25 ){
	load(paste(getwd(),'/output',r,'.RData',sep=''))
	beta.sas1 <- betas.sas
	beta.uni1 <- betas.uni
	beta.eff1 <- betas.eff
	beta.pois1 <- betas.pois
	
	load(paste(getwd(),'/test/output',r,'.RData',sep=''))
	aaa <- c(aaa,all(beta.sas1 == betas.sas))
	aaa <- c(aaa,all(beta.uni1 == betas.uni))
	aaa <- c(aaa,all(beta.eff1 == betas.eff))
	aaa <- c(aaa,all(beta.pois1 == betas.pois))
	
}
