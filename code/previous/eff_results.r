# eff results
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results')

bias.out <- NULL
mse.out <- NULL
cov.out <- NULL
fail.out <- NULL

for (n in c(50,100,200)){
	load(paste(getwd(),'/output',n,'.RData',sep=''))
	# efficient
	bias <- colMeans(betas.eff) - b
	v    <- apply(betas.eff,2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	cov.out <- rbind(cov.out,coverage)
	
	# independent
	bias <- colMeans(betas.uni) - b
	v    <- apply(betas.uni,2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	cov.out <- rbind(cov.out,coverage.u)
	
	# gen mod
	fail.out <- c(fail.out , sum(betas.sas==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
	v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	#cov.out <- rbind(cov.out,coverage)
	
	# poisson
	fail.out <- c(fail.out , sum(betas.pois==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
	v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
}

c1 <- rep(c('efficient','independent'),4)

header <- c('estimator','$\\beta_1$','$\\beta_2$','$\\beta_3$','$\\beta_4$','$\\beta_5$')
caption <- 'Bias'
f <- paste(getwd(),'/bias.tex',sep='')
label <- NULL
rowname <- c('50','','100','','200','','500','')
my.tex(file=f,table=cbind(c1,bias.out),sideways=F,center=T,align=NULL,rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='sample size')

header <- c('estimator','$\\beta_1$','$\\beta_2$','$\\beta_3$','$\\beta_4$','$\\beta_5$')
caption <- 'MSE $10^{-3}$'
f <- paste(getwd(),'/mse.tex',sep='')
label <- NULL
my.tex(file=f,table=cbind(c1,round(10^3*mse.out,2)),sideways=F,center=T,align=NULL,rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='sample size')

header <- c('estimator','$\\beta_1$','$\\beta_2$','$\\beta_3$','$\\beta_4$','$\\beta_5$')
caption <- 'Coverage'
f <- paste(getwd(),'/coverage.tex',sep='')
label <- NULL
my.tex(file=f,table=cbind(c1,round(100*cov.out,2)),sideways=F,center=T,align=NULL,rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='sample size')



###############################
#######3 yelland ###############

# eff results
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results\\yelland')
source('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\my_tex.r')

fail.gee <- fail.pois <- NULL
fail.row <- c(7,11,5,9)
for (ii in 0:8){
	for( r in fail.row ){
		load(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
	}
}
fail.gee <- matrix(fail.gee,nrow=9,byrow=T)
fail.pois <- matrix(fail.pois,nrow=9,byrow=T)


header <- c('Covariate RR','Interventional','','Observational','')
caption <- 'Convergence rate'
f <- paste(getwd(),'/convergence.tex',sep='')
label <- NULL
rowname <- c('','1','','','1.25','','','2','','')
my.tex(file=f,table=rbind(c('','20 Clusters' , '50 Clusters','20 Clusters' , '50 Clusters'),cbind(rep(c(1,1.25,2),3),1-fail.gee)),sideways=F,center=T,align='|l|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(2,5,8,11),caption=caption,label=label,rownames=rowname,row.header='Treatment RR')





bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
bias.pois.eff <- mse.pois.eff <- cov.pois.eff <- NULL
rows <- c(8,7,6,5)
for (ii in 0:8){
	for ( r in rows ){
		load(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))
		# efficient
		bias <- colMeans(betas.eff) - b
		v    <- apply(betas.eff,2,var)
		mse  <- bias^2 + v 
		bias.eff <- c(bias.eff,bias[1])
		mse.eff <- c(mse.eff,mse[1])
		cov.eff <- c(cov.eff,coverage[1])
		
		# pois efficient
		bias <- colMeans(betas.pois.eff) - b
		v    <- apply(betas.pois.eff,2,var)
		mse  <- bias^2 + v 
		bias.pois.eff <- c(bias.pois.eff,bias[1])
		mse.pois.eff <- c(mse.pois.eff,mse[1])
		cov.pois.eff <- c(cov.pois.eff,coverage[1])
		
		# independent
		bias <- colMeans(betas.uni) - b
		v    <- apply(betas.uni,2,var)
		mse  <- bias^2 + v 
		bias.uni <- rbind(bias.uni,bias[1])
		mse.uni <- rbind(mse.uni,mse[1])
		cov.uni <- rbind(cov.uni,coverage.u[1])
	
		# gen mod
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.gee <- rbind(bias.gee,bias[1])
		mse.gee <- rbind(mse.gee,mse[1])
	
		# poisson
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.pois <- rbind(bias.pois,bias[1])
		mse.pois <- rbind(mse.pois,mse[1])
	}
}
for (ii in 0:8){
	for ( r in rows ){
		load(paste(getwd(),'/output',r+8*ii+4,'.RData',sep=''))
		# efficient
		bias <- colMeans(betas.eff) - b
		v    <- apply(betas.eff,2,var)
		mse  <- bias^2 + v 
		bias.eff <- c(bias.eff,bias[1])
		mse.eff <- c(mse.eff,mse[1])
		cov.eff <- c(cov.eff,coverage[1])
		
		# pois efficient
		bias <- colMeans(betas.pois.eff) - b
		v    <- apply(betas.pois.eff,2,var)
		mse  <- bias^2 + v 
		bias.pois.eff <- c(bias.pois.eff,bias[1])
		mse.pois.eff <- c(mse.pois.eff,mse[1])
		cov.pois.eff <- c(cov.pois.eff,coverage[1])
		
		# independent
		bias <- colMeans(betas.uni) - b
		v    <- apply(betas.uni,2,var)
		mse  <- bias^2 + v 
		bias.uni <- rbind(bias.uni,bias[1])
		mse.uni <- rbind(mse.uni,mse[1])
		cov.uni <- rbind(cov.uni,coverage.u[1])
	
		# gen mod
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.gee <- rbind(bias.gee,bias[1])
		mse.gee <- rbind(mse.gee,mse[1])
	
		# poisson
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.pois <- rbind(bias.pois,bias[1])
		mse.pois <- rbind(mse.pois,mse[1])
	}
}

mse.eff <- matrix(mse.eff,nrow=18,byrow=T)
mse.uni <- matrix(mse.uni,nrow=18,byrow=T)
mse.gee <- matrix(mse.gee,nrow=18,byrow=T)
mse.pois <- matrix(mse.pois,nrow=18,byrow=T)
mse.pois.eff <- matrix(mse.pois.eff,nrow=18,byrow=T)

bias.eff <- matrix(bias.eff,nrow=18,byrow=T)
bias.uni <- matrix(bias.uni,nrow=18,byrow=T)
bias.gee <- matrix(bias.gee,nrow=18,byrow=T)
bias.pois <- matrix(bias.pois,nrow=18,byrow=T)
bias.pois.eff <- matrix(bias.pois.eff,nrow=18,byrow=T)


header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - GEE'
f <- paste(getwd(),'/bias-gee.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.gee,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')

header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - Efficient'
f <- paste(getwd(),'/bias-eff.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.eff,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')

header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - Assume independence'
f <- paste(getwd(),'/bias-uni.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.uni,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')

header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - poisson'
f <- paste(getwd(),'/bias-pois.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.pois,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')




(mse.eff-bias.eff^2)<(mse.pois-bias.pois^2)







###########################################
###########################################






for (n in 5:76){
	load(paste(getwd(),'/output',n,'.RData',sep=''))
	# efficient
	bias <- colMeans(betas.eff) - b
	v    <- apply(betas.eff,2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	cov.out <- rbind(cov.out,coverage)
	
	# independent
	bias <- colMeans(betas.uni) - b
	v    <- apply(betas.uni,2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	cov.out <- rbind(cov.out,coverage.u)
	
	# gen mod
	fail.out <- c(fail.out , sum(betas.sas==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
	v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	#cov.out <- rbind(cov.out,coverage)
	
	# poisson
	fail.out <- c(fail.out , sum(betas.pois==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
	v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
}







#################################
################################
###############################
###############################
####### new ones ###############

# eff results
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results\\results\\p0=.2varre=.2')
source('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\my_tex.r')

fail.gee <- fail.pois <- NULL
fail.row <- c(79,83,77,81)
for (ii in 0:8){
	for( r in fail.row ){
		if (file.exists(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))){
			load(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))
			betas.sas <- betas.sas[betas.sas[,1]!=0,]
			fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
			betas.pois <- betas.pois[betas.pois[,1]!=0,]
			fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
		}
	}
}
fail.gee <- matrix(fail.gee,nrow=9,byrow=T)
fail.pois <- matrix(fail.pois,nrow=9,byrow=T)


header <- c('Covariate RR','Interventional','','Observational','')
caption <- 'Convergence rate'
f <- paste(getwd(),'/convergence.tex',sep='')
label <- NULL
rowname <- c('','1','','','1.25','','','2','','')
my.tex(file=f,table=rbind(c('','20 Clusters' , '50 Clusters','20 Clusters' , '50 Clusters'),cbind(rep(c(1,1.25,2),3),1-fail.gee)),sideways=F,center=T,align='|l|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(2,5,8,11),caption=caption,label=label,rownames=rowname,row.header='Treatment RR')





bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
rows <- c(79,83,77,81)
temp <- NULL
for (ii in 0:8){
	for ( r in rows ){
		if (file.exists(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))){
		temp <- c(temp,r+8*ii)
		load(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))
		# efficient
		bias <- colMeans(betas.eff) - b
		v    <- apply(betas.eff,2,var)
		mse  <- bias^2 + v 
		bias.eff <- c(bias.eff,bias[1])
		mse.eff <- c(mse.eff,mse[1])
		cov.eff <- c(cov.eff,coverage[1])
		
		# independent
		bias <- colMeans(betas.uni) - b
		v    <- apply(betas.uni,2,var)
		mse  <- bias^2 + v 
		bias.uni <- rbind(bias.uni,bias[1])
		mse.uni <- rbind(mse.uni,mse[1])
		cov.uni <- rbind(cov.uni,coverage.u[1])
	
		# gen mod
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.gee <- rbind(bias.gee,bias[1])
		mse.gee <- rbind(mse.gee,mse[1])
	
		# poisson
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.pois <- rbind(bias.pois,bias[1])
		mse.pois <- rbind(mse.pois,mse[1])
		}
	}
}
for (ii in 0:8){
	for ( r in rows ){
		if (file.exists(paste(getwd(),'/output',r+8*ii,'.RData',sep=''))){
		load(paste(getwd(),'/output',r+8*ii+4,'.RData',sep=''))
		# efficient
		bias <- colMeans(betas.eff) - b
		v    <- apply(betas.eff,2,var)
		mse  <- bias^2 + v 
		bias.eff <- c(bias.eff,bias[1])
		mse.eff <- c(mse.eff,mse[1])
		cov.eff <- c(cov.eff,coverage[1])
		
		# independent
		bias <- colMeans(betas.uni) - b
		v    <- apply(betas.uni,2,var)
		mse  <- bias^2 + v 
		bias.uni <- rbind(bias.uni,bias[1])
		mse.uni <- rbind(mse.uni,mse[1])
		cov.uni <- rbind(cov.uni,coverage.u[1])
	
		# gen mod
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.gee <- rbind(bias.gee,bias[1])
		mse.gee <- rbind(mse.gee,mse[1])
	
		# poisson
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
		v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
		mse  <- bias^2 + v 
		bias.pois <- rbind(bias.pois,bias[1])
		mse.pois <- rbind(mse.pois,mse[1])
		}
	}
}

mse.eff <- matrix(mse.eff,nrow=18,byrow=T)
mse.uni <- matrix(mse.uni,nrow=18,byrow=T)
mse.gee <- matrix(mse.gee,nrow=18,byrow=T)
mse.pois <- matrix(mse.pois,nrow=18,byrow=T)

bias.eff <- matrix(bias.eff,nrow=18,byrow=T)
bias.uni <- matrix(bias.uni,nrow=18,byrow=T)
bias.gee <- matrix(bias.gee,nrow=18,byrow=T)
bias.pois <- matrix(bias.pois,nrow=18,byrow=T)


header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - GEE'
f <- paste(getwd(),'/bias-gee.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.gee,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')

header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - Efficient'
f <- paste(getwd(),'/bias-eff.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.eff,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')

header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - Assume independence'
f <- paste(getwd(),'/bias-uni.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.uni,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')

header <- c('Trt RR','Cov RR','Interventional','','Observational','')
caption <- 'Bias - poisson'
f <- paste(getwd(),'/bias-pois.tex',sep='')
label <- NULL
rowname <- c('','20',rep('',8),'50',rep('',8))
my.tex(file=f,table=rbind(c('','','Binary','Continuous','Binary','Continuous'),cbind(c(1,'','',1.25,'','',2,'',''),rep(c(1,1.25,2),3),round(1000*bias.pois,2))),sideways=F,center=T,align='|l|c|c|cc|cc|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,2,11),caption=caption,label=label,rownames=rowname,row.header='No of Cluster')




(mse.eff-bias.eff^2)<(mse.pois-bias.pois^2)







###########################################
###########################################






for (n in 5:76){
	load(paste(getwd(),'/output',n,'.RData',sep=''))
	# efficient
	bias <- colMeans(betas.eff) - b
	v    <- apply(betas.eff,2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	cov.out <- rbind(cov.out,coverage)
	
	# independent
	bias <- colMeans(betas.uni) - b
	v    <- apply(betas.uni,2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	cov.out <- rbind(cov.out,coverage.u)
	
	# gen mod
	fail.out <- c(fail.out , sum(betas.sas==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.sas),ncol=2)) - b
	v    <- apply(matrix(as.numeric(betas.sas),ncol=2),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	#cov.out <- rbind(cov.out,coverage)
	
	# poisson
	fail.out <- c(fail.out , sum(betas.pois==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.pois),ncol=2)) - b
	v    <- apply(matrix(as.numeric(betas.pois),ncol=2),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
}








#################################
################################
###############################
###############################
####### 1-4 ###############

# eff results
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results\\results')
source('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\my_tex.r')

fail.gee <- fail.pois <- NULL
index <- 1:3
for( r in index ){
		if (file.exists(paste(getwd(),'/output',r,'.RData',sep=''))){
			load(paste(getwd(),'/output',r,'.RData',sep=''))
			betas.sas <- betas.sas[betas.sas[,1]!=0,]
			fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
			betas.pois <- betas.pois[betas.pois[,1]!=0,]
			fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
		}
}

bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
rows <- 1:3
temp <- NULL
for ( r in rows ){
		load(paste(getwd(),'/output',r,'.RData',sep=''))
		# efficient
		bias <- colMeans(betas.eff) - b
		v    <- apply(betas.eff,2,var)
		mse  <- bias^2 + v 
		bias.eff <- rbind(bias.eff,bias)
		mse.eff <- rbind(mse.eff,mse)
		cov.eff <- rbind(cov.eff,coverage)
		
		# independent
		bias <- colMeans(betas.uni) - b
		v    <- apply(betas.uni,2,var)
		mse  <- bias^2 + v 
		bias.uni <- rbind(bias.uni,bias)
		mse.uni <- rbind(mse.uni,mse)
		cov.uni <- rbind(cov.uni,coverage.u)
	
		# gen mod
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
		v    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
		mse  <- bias^2 + v 
		bias.gee <- rbind(bias.gee,bias)
		mse.gee <- rbind(mse.gee,mse)
	
		# poisson
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
		v    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
		mse  <- bias^2 + v 
		bias.pois <- rbind(bias.pois,bias)
		mse.pois <- rbind(mse.pois,mse)
}



#################################
################################
###############################
###############################
####### 500's and 1000's ###############

# eff results
rm(list=ls())
setwd("/Users/matt/Documents/Harvard/Research/Correlated Binary/results/clustersize5/1000/continuous/exchangeable")
source('/Users/matt/Documents/Harvard/Research/my_tex.r')
fail.gee <- fail.pois <- NULL
index <- 1:20
for( r in index ){
		if (file.exists(paste(getwd(),'/output',500+r,'.RData',sep=''))){
			load(paste(getwd(),'/output',500+r,'.RData',sep=''))
			betas.sas <- betas.sas[betas.sas[,1]!=0,]
			fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
			betas.pois <- betas.pois[betas.pois[,1]!=0,]
			fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
		}
}
bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
bias.pois.eff <- mse.pois.eff <- cov.pois.eff <- NULL
rows <- 1:20
temp <- NULL
my.beta <- NULL
for ( r in rows ){
		load(paste(getwd(),'/output',500+r,'.RData',sep=''))
		my.beta <- rbind(my.beta,b)
		# efficient
		bias <- colMeans(betas.eff) - b
		v    <- apply(betas.eff,2,var)
		mse  <- bias^2 + v 
		bias.eff <- rbind(bias.eff,bias)
		mse.eff <- rbind(mse.eff,mse)
		cov.eff <- rbind(cov.eff,coverage)
		
		# update from poisson
		bias <- colMeans(betas.pois.eff) - b
		v    <- apply(betas.pois.eff,2,var)
		mse  <- bias^2 + v 
		bias.pois.eff <- rbind(bias.pois.eff,bias)
		mse.pois.eff <- rbind(mse.pois.eff,mse)
		coverage <- numeric(q)
		for (j in 1:q){
			coverage[j] <- sum((ci.pois.eff[[j]][,1]<b[j]) & (ci.pois.eff[[j]][,2]>b[j]),na.rm=T)/M
		}
		cov.pois.eff <- rbind(cov.pois.eff,coverage)
		
		# independent
		bias <- colMeans(betas.uni) - b
		v    <- apply(betas.uni,2,var)
		mse  <- bias^2 + v 
		bias.uni <- rbind(bias.uni,bias)
		mse.uni <- rbind(mse.uni,mse)
		cov.uni <- rbind(cov.uni,coverage.u)
	
		# gen mod
		betas.sas <- betas.sas[betas.sas[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
		v    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
		mse  <- bias^2 + v 
		bias.gee <- rbind(bias.gee,bias)
		mse.gee <- rbind(mse.gee,mse)
	
		# poisson
		betas.pois <- betas.pois[betas.pois[,1]!=0,]
		bias <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
		v    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
		mse  <- bias^2 + v 
		bias.pois <- rbind(bias.pois,bias)
		mse.pois <- rbind(mse.pois,mse)
		coverage <- numeric(q)
		for (j in 1:q){
			coverage[j] <- sum((as.numeric(ci.pois[[j]][,1])<b[j]) & (as.numeric(ci.pois[[j]][,2])>b[j]))/M
		}
		cov.pois <- rbind(cov.pois,coverage)
		
}
tab.true <- round(cbind(exp(my.beta[c(1:20),1]),mse.uni[c(1:20),1]*1000,mse.pois[c(1:20),1]*1000,mse.pois.eff[c(1:20),1]*1000,exp(my.beta[c(1:20),5]),mse.uni[c(1:20),5]*1000,mse.pois[c(1:20),5]*1000,mse.pois.eff[c(1:20),5]*1000,fail.gee[c(1:20)]*100),3)
tab.true <- cbind(tab.true,cov.pois.eff[,c(1,5)]*100)

setwd("/Users/matt/Documents/Harvard/Research/Correlated Binary/results/clustersize5/1000/continuous/household")
fail.gee <- fail.pois <- NULL
index <- 1:20
for( r in index ){
   if (file.exists(paste(getwd(),'/output',500+r,'.RData',sep=''))){
      load(paste(getwd(),'/output',500+r,'.RData',sep=''))
      betas.sas <- betas.sas[betas.sas[,1]!=0,]
      fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
      betas.pois <- betas.pois[betas.pois[,1]!=0,]
      fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
   }
}
bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
bias.pois.eff <- mse.pois.eff <- cov.pois.eff <- NULL
rows <- 1:20
temp <- NULL
my.beta <- NULL
for ( r in rows ){
   load(paste(getwd(),'/output',500+r,'.RData',sep=''))
   my.beta <- rbind(my.beta,b)
   # efficient
   bias <- colMeans(betas.eff) - b
   v    <- apply(betas.eff,2,var)
   mse  <- bias^2 + v 
   bias.eff <- rbind(bias.eff,bias)
   mse.eff <- rbind(mse.eff,mse)
   cov.eff <- rbind(cov.eff,coverage)
   
   # update from poisson
   bias <- colMeans(betas.pois.eff) - b
   v    <- apply(betas.pois.eff,2,var)
   mse  <- bias^2 + v 
   bias.pois.eff <- rbind(bias.pois.eff,bias)
   mse.pois.eff <- rbind(mse.pois.eff,mse)
   coverage <- numeric(q)
   for (j in 1:q){
      coverage[j] <- sum((ci.pois.eff[[j]][,1]<b[j]) & (ci.pois.eff[[j]][,2]>b[j]),na.rm=T)/M
   }
   cov.pois.eff <- rbind(cov.pois.eff,coverage)
   
   # independent
   bias <- colMeans(betas.uni) - b
   v    <- apply(betas.uni,2,var)
   mse  <- bias^2 + v 
   bias.uni <- rbind(bias.uni,bias)
   mse.uni <- rbind(mse.uni,mse)
   cov.uni <- rbind(cov.uni,coverage.u)
   
   # gen mod
   betas.sas <- betas.sas[betas.sas[,1]!=0,]
   bias <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
   v    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
   mse  <- bias^2 + v 
   bias.gee <- rbind(bias.gee,bias)
   mse.gee <- rbind(mse.gee,mse)
   
   # poisson
   betas.pois <- betas.pois[betas.pois[,1]!=0,]
   bias <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
   v    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
   mse  <- bias^2 + v 
   bias.pois <- rbind(bias.pois,bias)
   mse.pois <- rbind(mse.pois,mse)
   coverage <- numeric(q)
   for (j in 1:q){
      coverage[j] <- sum((as.numeric(ci.pois[[j]][,1])<b[j]) & (as.numeric(ci.pois[[j]][,2])>b[j]))/M
   }
   cov.pois <- rbind(cov.pois,coverage)
   
}
tab.miss <- round(cbind(mse.uni[c(1:20),1]*1000,mse.pois[c(1:20),1]*1000,mse.pois.eff[c(1:20),1]*1000,mse.uni[c(1:20),5]*1000,mse.pois[c(1:20),5]*1000,mse.pois.eff[c(1:20),5]*1000,fail.gee[c(1:20)]*100),3)
tab.miss <- cbind('',tab.miss[,1:3],'',tab.miss[,4:7],cov.pois.eff[,c(1,5)]*100)

setwd("/Users/matt/Documents/Harvard/Research/Correlated Binary/results/clustersize5/1000/binary/exchangeable")
fail.gee <- fail.pois <- NULL
index <- 1:20
for( r in index ){
   if (file.exists(paste(getwd(),'/output',500+r,'.RData',sep=''))){
      load(paste(getwd(),'/output',500+r,'.RData',sep=''))
      betas.sas <- betas.sas[betas.sas[,1]!=0,]
      fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
      betas.pois <- betas.pois[betas.pois[,1]!=0,]
      fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
   }
}
bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
bias.pois.eff <- mse.pois.eff <- cov.pois.eff <- NULL
rows <- 1:20
temp <- NULL
my.beta <- NULL
for ( r in rows ){
   load(paste(getwd(),'/output',500+r,'.RData',sep=''))
   my.beta <- rbind(my.beta,b)
   # efficient
   bias <- colMeans(betas.eff) - b
   v    <- apply(betas.eff,2,var)
   mse  <- bias^2 + v 
   bias.eff <- rbind(bias.eff,bias)
   mse.eff <- rbind(mse.eff,mse)
   cov.eff <- rbind(cov.eff,coverage)
   
   # update from poisson
   bias <- colMeans(betas.pois.eff) - b
   v    <- apply(betas.pois.eff,2,var)
   mse  <- bias^2 + v 
   bias.pois.eff <- rbind(bias.pois.eff,bias)
   mse.pois.eff <- rbind(mse.pois.eff,mse)
   coverage <- numeric(q)
   for (j in 1:q){
      coverage[j] <- sum((ci.pois.eff[[j]][,1]<b[j]) & (ci.pois.eff[[j]][,2]>b[j]),na.rm=T)/M
   }
   cov.pois.eff <- rbind(cov.pois.eff,coverage)
   
   # independent
   bias <- colMeans(betas.uni) - b
   v    <- apply(betas.uni,2,var)
   mse  <- bias^2 + v 
   bias.uni <- rbind(bias.uni,bias)
   mse.uni <- rbind(mse.uni,mse)
   cov.uni <- rbind(cov.uni,coverage.u)
   
   # gen mod
   betas.sas <- betas.sas[betas.sas[,1]!=0,]
   bias <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
   v    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
   mse  <- bias^2 + v 
   bias.gee <- rbind(bias.gee,bias)
   mse.gee <- rbind(mse.gee,mse)
   
   # poisson
   betas.pois <- betas.pois[betas.pois[,1]!=0,]
   bias <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
   v    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
   mse  <- bias^2 + v 
   bias.pois <- rbind(bias.pois,bias)
   mse.pois <- rbind(mse.pois,mse)
   coverage <- numeric(q)
   for (j in 1:q){
      coverage[j] <- sum((as.numeric(ci.pois[[j]][,1])<b[j]) & (as.numeric(ci.pois[[j]][,2])>b[j]))/M
   }
   cov.pois <- rbind(cov.pois,coverage)
   
}
tab.true <- round(cbind(exp(my.beta[c(1:20),1]),mse.uni[c(1:20),1]*1000,mse.pois[c(1:20),1]*1000,mse.pois.eff[c(1:20),1]*1000,exp(my.beta[c(1:20),5]),mse.uni[c(1:20),5]*1000,mse.pois[c(1:20),5]*1000,mse.pois.eff[c(1:20),5]*1000,fail.gee[c(1:20)]*100),3)
tab.true <- cbind(tab.true,cov.pois.eff[,c(1,5)]*100)
tab.true.b <- round(cbind(bias.uni[c(1:20),1]*1000,bias.pois[c(1:20),1]*1000,bias.pois.eff[c(1:20),1]*1000,bias.uni[c(1:20),5]*1000,bias.pois[c(1:20),5]*1000,bias.pois.eff[c(1:20),5]*1000,fail.gee[c(1:20)]*100),3)
tab.true.b <- cbind('',tab.true.b[,1:3],'',tab.true.b[,4:7],cov.pois.eff[,c(1,5)]*100)


setwd("/Users/matt/Documents/Harvard/Research/Correlated Binary/results/clustersize5/1000/binary/household")
fail.gee <- fail.pois <- NULL
index <- 1:20
for( r in index ){
   if (file.exists(paste(getwd(),'/output',500+r,'.RData',sep=''))){
      load(paste(getwd(),'/output',500+r,'.RData',sep=''))
      betas.sas <- betas.sas[betas.sas[,1]!=0,]
      fail.gee <- c(fail.gee , 1-dim(betas.sas)[1]/M)
      betas.pois <- betas.pois[betas.pois[,1]!=0,]
      fail.pois <- c(fail.pois , 1-dim(betas.pois)[1]/M)
   }
}
bias.eff <- mse.eff <- cov.eff <- NULL
bias.uni <- mse.uni <- cov.uni <- NULL
bias.gee <- mse.gee <- cov.gee <- NULL
bias.pois <- mse.pois <- cov.pois <- NULL
bias.pois.eff <- mse.pois.eff <- cov.pois.eff <- NULL
rows <- 1:20
temp <- NULL
my.beta <- NULL
for ( r in rows ){
   load(paste(getwd(),'/output',500+r,'.RData',sep=''))
   my.beta <- rbind(my.beta,b)
   # efficient
   bias <- colMeans(betas.eff) - b
   v    <- apply(betas.eff,2,var)
   mse  <- bias^2 + v 
   bias.eff <- rbind(bias.eff,bias)
   mse.eff <- rbind(mse.eff,mse)
   cov.eff <- rbind(cov.eff,coverage)
   
   # update from poisson
   bias <- colMeans(betas.pois.eff) - b
   v    <- apply(betas.pois.eff,2,var)
   mse  <- bias^2 + v 
   bias.pois.eff <- rbind(bias.pois.eff,bias)
   mse.pois.eff <- rbind(mse.pois.eff,mse)
   coverage <- numeric(q)
   for (j in 1:q){
      coverage[j] <- sum((ci.pois.eff[[j]][,1]<b[j]) & (ci.pois.eff[[j]][,2]>b[j]),na.rm=T)/M
   }
   cov.pois.eff <- rbind(cov.pois.eff,coverage)
   
   # independent
   bias <- colMeans(betas.uni) - b
   v    <- apply(betas.uni,2,var)
   mse  <- bias^2 + v 
   bias.uni <- rbind(bias.uni,bias)
   mse.uni <- rbind(mse.uni,mse)
   cov.uni <- rbind(cov.uni,coverage.u)
   
   # gen mod
   betas.sas <- betas.sas[betas.sas[,1]!=0,]
   bias <- colMeans(matrix(as.numeric(betas.sas),ncol=5)) - b
   v    <- apply(matrix(as.numeric(betas.sas),ncol=5),2,var)
   mse  <- bias^2 + v 
   bias.gee <- rbind(bias.gee,bias)
   mse.gee <- rbind(mse.gee,mse)
   
   # poisson
   betas.pois <- betas.pois[betas.pois[,1]!=0,]
   bias <- colMeans(matrix(as.numeric(betas.pois),ncol=5)) - b
   v    <- apply(matrix(as.numeric(betas.pois),ncol=5),2,var)
   mse  <- bias^2 + v 
   bias.pois <- rbind(bias.pois,bias)
   mse.pois <- rbind(mse.pois,mse)
   coverage <- numeric(q)
   for (j in 1:q){
      coverage[j] <- sum((as.numeric(ci.pois[[j]][,1])<b[j]) & (as.numeric(ci.pois[[j]][,2])>b[j]))/M
   }
   cov.pois <- rbind(cov.pois,coverage)
   
}
tab.miss <- round(cbind(mse.uni[c(1:20),1]*1000,mse.pois[c(1:20),1]*1000,mse.pois.eff[c(1:20),1]*1000,mse.uni[c(1:20),5]*1000,mse.pois[c(1:20),5]*1000,mse.pois.eff[c(1:20),5]*1000,fail.gee[c(1:20)]*100),3)
tab.miss <- cbind('',tab.miss[,1:3],'',tab.miss[,4:7],cov.pois.eff[,c(1,5)]*100)
tab.miss.b <- round(cbind(bias.uni[c(1:20),1]*1000,bias.pois[c(1:20),1]*1000,bias.pois.eff[c(1:20),1]*1000,bias.uni[c(1:20),5]*1000,bias.pois[c(1:20),5]*1000,bias.pois.eff[c(1:20),5]*1000,fail.gee[c(1:20)]*100),3)
tab.miss.b <- cbind('',tab.miss.b[,1:3],'',tab.miss.b[,4:7],cov.pois.eff[,c(1,5)]*100)





################################################

# exp(my.beta[c(1:3,5:20),1]),

tab <- NULL
for (i in 1:dim(tab.miss)[1]){
   tab <- rbind(tab,tab.true[i,],tab.miss[i,])
}

tab.b <- NULL
for (i in 1:dim(tab.miss)[1]){
   tab.b <- rbind(tab.b,tab.true.b[i,],tab.miss.b[i,])
}


header <- c('Relative risk','Modified Poisson','Efficient','Relative risk','Modified Poisson','Efficient','Failure rate')
caption <- 'MSE ($10^{-3}$)'
f <- './mse.tex'
label <- NULL
rowname <- rep(c('\\rule{0pt}{1.5ex}Exch.','\\rule[-1.5ex]{0pt}{0pt}Other'),100)
my.tex(file=f,table=tab,sideways=F,center=T,align='|l|l|cc|l|cc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,7,15,23,31),caption=caption,label=label,rownames=rowname,row.header='True CS')



## smaller tables 
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results\\final\\binary')

tab.binary <- tab[1:8,c(5:8,11)]
header <- c('Relative risk','Independent','Modified Poisson','Efficient','Coverage')
caption <- 'Mean square error ($10^{-3}$) of the modified Poisson approach and the efficient approach for estimating the relative risk of a binary covariate when there are 500 clusters of size 5 under both correct specification and incorrect specification of the working correlation structure. In all scenarios, the working correlation structure is assumed to be exchangeable, while the true correlation structure is either exchangeable with $\\rho=.3$ or the household structure given by Equation~\\ref{eqn:house}.'
f <- './mse-binary.tex'
label <- 'tab:binary'
rowname <- rep(c('\\rule{0pt}{1.7ex}Exch.','\\rule[-1.7ex]{0pt}{0pt}Household'),100)
my.tex(file=f,table=tab.binary,sideways=F,center=T,align='|l|l|ccc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='True CS')


tab.cont <- tab[sort(rep(0:3,2)*8+13)+c(0,1),c(1:4,10)]
header <- c('Relative risk','Independent','Modified Poisson','Efficient','Coverage')
caption <- 'Mean square error ($10^{-3}$) of the modified Poisson approach and the efficient approach for estimating the relative risk of a continuous covariate when there are 500 clusters of size 5 under both correct specification and incorrect specification of the working correlation structure. In all scenarios, the working correlation structure is assumed to be exchangeable, while the true correlation structure is either exchangeable with $\\rho=.3$ or the household structure given by Equation~\\ref{eqn:house}'
f <- './mse-cont.tex'
label <- 'tab:cont'
rowname <- rep(c('\\rule{0pt}{1.7ex}Exch.','\\rule[-1.7ex]{0pt}{0pt}Household'),100)
my.tex(file=f,table=tab.cont,sideways=F,center=T,align='|l|l|ccc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='True CS')


tab.exch.binary <- tab.binary
tab.binary.b <- tab.b[1:8,c(5:8,11)]
tab.exch.binary[c(2,4,6,8),] = tab.binary.b[c(1,3,5,7),]

tab.house.binary <- tab.binary
tab.house.binary[c(1,3,5,7),] = tab.house.binary[c(2,4,6,8),]
tab.house.binary[c(2,4,6,8),] = tab.binary.b[c(2,4,6,8),]
tab.house.binary[,1] = c("1","","1.05","","1.5","","2","")

setwd("/Users/matt/Documents/Harvard/Research/Correlated Binary/results/")

header <- c('Relative risk','Independent','Modified Poisson','Efficient','Coverage')
caption <- 'Bias ($10^{-3}$) and mean square error ($10^{-3}$) of the modified Poisson approach and the efficient approach for estimating the relative risk of a binary covariate when there are 1000 clusters of size 5 under correct specification of the working correlation structure. In all scenarios, the true correlation structure is exchangeable with $\\rho=0.3$.'
f <- './mse-exch-binary.tex'
label <- 'tab:exchbinary'
rowname <- rep(c('\\rule{0pt}{1.7ex}MSE','\\rule[-1.7ex]{0pt}{0pt}Bias'),100)
my.tex(file=f,table=tab.exch.binary,sideways=F,center=T,align='|l|l|ccc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='True CS')

header <- c('Relative risk','Independent','Modified Poisson','Efficient','Coverage')
caption <- 'Bias ($10^{-3}$) and mean square error ($10^{-3}$) of the modified Poisson approach and the efficient approach for estimating the relative risk of a binary covariate when there are 1000 clusters of size 5 under incorrect specification of the working correlation structure. In all scenarios, the working correlation structure is assumed to be exchangeable, while the true correlation structure is the household structure given by Equation~\\ref{eqn:house}.'
f <- './mse-house-binary.tex'
label <- 'tab:housebinary'
rowname <- rep(c('\\rule{0pt}{1.7ex}MSE','\\rule[-1.7ex]{0pt}{0pt}Bias'),100)
my.tex(file=f,table=tab.house.binary,sideways=F,center=T,align='|l|l|ccc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='True CS')



tab.exch.cont <- tab.cont
tab.cont.b <- tab.b[sort(rep(0:3,2)*8+13)+c(0,1),c(1:4,10)]
tab.exch.cont[c(2,4,6,8),] = tab.cont.b[c(1,3,5,7),]

tab.house.cont <- tab.cont
tab.house.cont[c(1,3,5,7),] = tab.house.cont[c(2,4,6,8),]
tab.house.cont[c(2,4,6,8),] = tab.cont.b[c(2,4,6,8),]
tab.house.cont[,1] = c("1","","1.05","","1.5","","2","")

setwd("/Users/matt/Documents/Harvard/Research/Correlated Binary/results/")

header <- c('Relative risk','Independent','Modified Poisson','Efficient','Coverage')
caption <- 'Bias ($10^{-3}$) and mean square error ($10^{-3}$) of the modified Poisson approach and the efficient approach for estimating the relative risk of a continuous covariate when there are 1000 clusters of size 5 under correct specification of the working correlation structure. In all scenarios, the true correlation structure is exchangeable with $\\rho=0.3$.'
f <- './mse-exch-cont.tex'
label <- 'tab:exchcont'
rowname <- rep(c('\\rule{0pt}{1.7ex}MSE','\\rule[-1.7ex]{0pt}{0pt}Bias'),100)
my.tex(file=f,table=tab.exch.cont,sideways=F,center=T,align='|l|l|ccc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='True CS')

header <- c('Relative risk','Independent','Modified Poisson','Efficient','Coverage')
caption <- 'Bias ($10^{-3}$) and mean square error ($10^{-3}$) of the modified Poisson approach and the efficient approach for estimating the relative risk of a contnuous covariate when there are 1000 clusters of size 5 under incorrect specification of the working correlation structure. In all scenarios, the working correlation structure is assumed to be exchangeable, while the true correlation structure is the household structure given by Equation~\\ref{eqn:house}.'
f <- './mse-house-cont.tex'
label <- 'tab:housecont'
rowname <- rep(c('\\rule{0pt}{1.7ex}MSE','\\rule[-1.7ex]{0pt}{0pt}Bias'),100)
my.tex(file=f,table=tab.house.cont,sideways=F,center=T,align='|l|l|ccc|c|',rowcolors=c('',''),header=header,header.color='ltgray',hlines=c(0,1,3,5,7,9),caption=caption,label=label,rownames=rowname,row.header='True CS')


