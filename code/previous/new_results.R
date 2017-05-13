rm(list=ls())

bias.eff = matrix(0,20,5)
v.eff = matrix(0,20,5)
bias.eff.pois = matrix(0,20,5)
v.eff.pois = matrix(0,20,5)
bias.rr = matrix(0,20,5)
v.rr = matrix(0,20,5)
bias.pois = matrix(0,20,5)
v.pois = matrix(0,20,5)
bias = matrix(0,20,5)
v = matrix(0,20,5)

index <- 1:20
for( r in index ){
   load(paste(getwd(),'/results/clustersize5/larger mean/output',500+r,'.RData',sep=''))
   # efficient
   bias.eff[r,] <- colMeans(betas.eff) - b
   v.eff[r,]    <- apply(betas.eff,2,var)
   # eff / poisson
   bias.eff.pois[r,] <- colMeans(betas.pois.eff) - b
   v.eff.pois[r,]    <- apply(betas.pois.eff,2,var)
   # independent
   bias[r,] <- colMeans(betas.uni) - b
   v[r,]    <- apply(betas.uni,2,var)
   # gen mod
   bias.rr[r,] <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
   v.rr[r,]    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
   # poisson
   bias.pois[r,] <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
   v.pois[r,]    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
}



######################
bias.eff = matrix(0,41,5)
v.eff = matrix(0,41,5)
bias.eff.pois = matrix(0,41,5)
v.eff.pois = matrix(0,41,5)
bias.rr = matrix(0,41,5)
v.rr = matrix(0,41,5)
bias.pois = matrix(0,41,5)
v.pois = matrix(0,41,5)
bias = matrix(0,41,5)
v = matrix(0,41,5)




index <- 1:41
for( r in index ){
   load(paste(getwd(),'/results/vary_single_rr/output',r,'.RData',sep=''))
   eff = betas.pois.eff
   pois = betas.pois
   uni = betas.uni
   load(paste(getwd(),'/results/vary_single_rr/1/output',r,'.RData',sep=''))
   betas.pois.eff = rbind(eff , betas.pois.eff)
   betas.pois = rbind(pois, betas.pois)
   betas.uni = rbind(uni , betas.uni)
   # efficient
   bias.eff[r,] <- colMeans(betas.eff) - b
   v.eff[r,]    <- apply(betas.eff,2,var)
   # eff / poisson
   bias.eff.pois[r,] <- colMeans(betas.pois.eff) - b
   v.eff.pois[r,]    <- apply(betas.pois.eff,2,var)
   # independent
   bias[r,] <- colMeans(betas.uni) - b
   v[r,]    <- apply(betas.uni,2,var)
   # gen mod
   #bias.rr[r,] <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
   #v.rr[r,]    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
   # poisson
   bias.pois[r,] <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
   v.pois[r,]    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
}

v.eff.pois/v.pois

plot(v.eff.pois[,1]/v.pois[,1],type='l',ylim=c(.8,1),ylab='Relative efficiency',axes=F , xlab='Risk ratio',lty='dashed')
lines(v.eff.pois[,1]/v[,1])
axis(2)
axis(1,at=seq(1,41,length.out=11),labels=seq(1,1.1,.01))
box()
legend('topright',legend=c('Independent','Modified Poisson'),lty=c(1,2))

rbind((v.eff.pois/v.pois)[36,2:5],
(v.eff.pois/v)[36,2:5],
(v.pois/v)[36,2:5])




######################
bias.eff = matrix(0,41,5)
v.eff = matrix(0,41,5)
bias.eff.pois = matrix(0,41,5)
v.eff.pois = matrix(0,41,5)
bias.rr = matrix(0,41,5)
v.rr = matrix(0,41,5)
bias.pois = matrix(0,41,5)
v.pois = matrix(0,41,5)
bias = matrix(0,41,5)
v = matrix(0,41,5)




index <- 1:41
for( r in index ){
   load(paste(getwd(),'/results/vary_single_rr/household/output',r,'.RData',sep=''))
   eff = betas.pois.eff
   pois = betas.pois
   uni = betas.uni
   # efficient
   bias.eff[r,] <- colMeans(betas.eff) - b
   v.eff[r,]    <- apply(betas.eff,2,var)
   # eff / poisson
   bias.eff.pois[r,] <- colMeans(betas.pois.eff) - b
   v.eff.pois[r,]    <- apply(betas.pois.eff,2,var)
   # independent
   bias[r,] <- colMeans(betas.uni) - b
   v[r,]    <- apply(betas.uni,2,var)
   # gen mod
   #bias.rr[r,] <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
   #v.rr[r,]    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
   # poisson
   bias.pois[r,] <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
   v.pois[r,]    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
}

v.eff.pois/v.pois

plot(v.eff.pois[,1]/v.pois[,1],type='l',ylim=c(.8,1),ylab='Relative efficiency',axes=F , xlab='Risk ratio',lty='dashed')
lines(v.eff.pois[,1]/v[,1])
axis(2)
axis(1,at=seq(1,41,length.out=11),labels=seq(1,1.1,.01))
box()
legend('topright',legend=c('Independent','Modified Poisson'),lty=c(1,2))

rbind((v.eff.pois/v.pois)[36,2:5],
      (v.eff.pois/v)[36,2:5],
      (v.pois/v)[36,2:5])
