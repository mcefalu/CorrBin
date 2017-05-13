### new results ###

# eff results
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\results\\results')

bias.out <- NULL
mse.out <- NULL
cov.out <- NULL
fail.out <- NULL

for (i in 5:12){
	load(paste(getwd(),'/output',i,'.RData',sep=''))
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
	fail.out <- c(fail.out , sum(betas.sas[,1]==0)/M)
	bias <- colMeans(matrix(as.numeric(betas.sas),ncol=dim(bias.out)[2])) - b
	v    <- apply(matrix(as.numeric(betas.sas),ncol=dim(bias.out)[2]),2,var)
	mse  <- bias^2 + v 
	bias.out <- rbind(bias.out,bias)
	mse.out <- rbind(mse.out,mse)
	#cov.out <- rbind(cov.out,coverage)
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



#####################
#####################
### Failure table ###

int20 <- 0:8*8+7
int50 <- 0:8*8+11
obs20 <- 0:8*8+5
obs50 <- 0:8*8+9

fail.out <- matrix(0,length(obs20),4)

for (ii in 1:length(int20)){
	if (file.exists(paste(getwd(),'/output',int20[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',int20[ii],'.RData',sep=''))
		# gen mod
		fail.out[ii,1] <- (length(unique(betas.sas[,2]))-1)/M*100
	}else{
		fail.out[ii,1] <- 888
	}
}
for (ii in 1:length(int50)){
	if (file.exists(paste(getwd(),'/output',int50[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',int50[ii],'.RData',sep=''))
		# gen mod
		fail.out[ii,2] <- (length(unique(betas.sas[,2]))-1)/M*100
	}else{
		fail.out[ii,2] <- 888
	}
}
for (ii in 1:length(obs20)){
	if (file.exists(paste(getwd(),'/output',obs20[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',obs20[ii],'.RData',sep=''))
		# gen mod
		fail.out[ii,3] <- (length(unique(betas.sas[,2]))-1)/M*100
	}else{
		fail.out[ii,3] <- 888
	}
}
for (ii in 1:length(obs50)){
	if (file.exists(paste(getwd(),'/output',obs50[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',obs50[ii],'.RData',sep=''))
		# gen mod
		fail.out[ii,4] <- (length(unique(betas.sas[,2]))-1)/M*100
	}else{
		fail.out[ii,4] <- 888
	}
}


#################################
#################################
### bias table -- binary only ###

int20 <- 0:8*8+8
int50 <- 0:8*8+12
obs20 <- 0:8*8+6
obs50 <- 0:8*8+10

bias.b <- matrix(0,2*3*3,7)
bias.b[c(1,9),1] <- c(20,50)
colnames(bias.b) <- c('clusters','trt RR','cov RR','ind','eff','GEE','poisson')
bias.b[c(1,4,7,10,13,16),2] <- rep(c(1,1.25,2),2)
bias.b[,3] <- rep(c(1,1.25,2),6)

n1 <- length(int20)
for (ii in 1:length(int20)){
	if (file.exists(paste(getwd(),'/output',int20[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',int20[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.b[ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.b[ii,4] <- mean((exp(temp) - exp(b[1])))
		# poisson
		bias.b[ii,7] <- mean((unique(exp(matrix(as.numeric(betas.pois),ncol=2)[,1])) - exp(b[1])))
		# gen mod
		bias.b[ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1])) - exp(b[1])))
	}else{
		bias.b[ii,4] <- 888
		bias.b[ii,5] <- 888
		bias.b[ii,6] <- 888
	}
}
for (ii in 1:length(int50)){
	if (file.exists(paste(getwd(),'/output',int50[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',int50[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.b[n1+ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.b[n1+ii,4] <- mean((exp(temp) - exp(b[1])))
		# poisson
		bias.b[n1+ii,7] <- mean((unique(exp(matrix(as.numeric(betas.pois),ncol=2)[,1])) - exp(b[1])))
		# gen mod
		bias.b[n1+ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1]))- exp(b[1])))
	}else{
		bias.b[n1+ii,4] <- 888
		bias.b[n1+ii,5] <- 888
		bias.b[n1+ii,6] <- 888
	}
}


###########
#########################
#####################################

int20 <- 0:8*8+8
int50 <- 0:8*8+12
obs20 <- 0:8*8+6
obs50 <- 0:8*8+10

bias.o <- matrix(0,2*3*3,7)
bias.o[c(1,9),1] <- c(20,50)
colnames(bias.o) <- c('clusters','trt RR','cov RR','ind','eff','GEE','pois')
bias.o[c(1,4,7,10,13,16),2] <- rep(c(1,1.25,2),2)
bias.o[,3] <- rep(c(1,1.25,2),6)

n1 <- length(obs20)
for (ii in 1:length(obs20)){
	if (file.exists(paste(getwd(),'/output',obs20[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',obs20[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.o[ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.o[ii,4] <- mean((exp(temp) - exp(b[1])))
		# poisson
		bias.o[ii,7] <- mean((unique(exp(matrix(as.numeric(betas.pois),ncol=2)[,1])) - exp(b[1])))
		# gen mod
		bias.o[ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1])) - exp(b[1])))
	}else{
		bias.o[ii,4] <- 888
		bias.o[ii,5] <- 888
		bias.o[ii,6] <- 888
	}
}
for (ii in 1:length(obs50)){
	if (file.exists(paste(getwd(),'/output',obs50[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',obs50[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.o[n1+ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.o[n1+ii,4] <- mean((exp(temp) - exp(b[1])))
		# poisson
		bias.o[n1+ii,7] <- mean((unique(exp(matrix(as.numeric(betas.pois),ncol=2)[,1])) - exp(b[1])))
		# gen mod
		bias.o[n1+ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1]))- exp(b[1])))
	}else{
		bias.o[n1+ii,4] <- 888
		bias.o[n1+ii,5] <- 888
		bias.o[n1+ii,6] <- 888
	}
}




#################################
#################################
### mse table -- binary only ###

int20 <- 0:8*8+8
int50 <- 0:8*8+12
obs20 <- 0:8*8+6
obs50 <- 0:8*8+10

bias.b <- matrix(0,2*3*3,6)
bias.b[c(1,9),1] <- c(20,50)
colnames(bias.b) <- c('clusters','trt RR','cov RR','ind','eff','GEE')
bias.b[c(1,4,7,10,13,16),2] <- rep(c(1,1.25,2),2)
bias.b[,3] <- rep(c(1,1.25,2),6)

n1 <- length(int20)
for (ii in 1:length(int20)){
	if (file.exists(paste(getwd(),'/output',int20[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',int20[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.b[ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.b[ii,4] <- mean((exp(temp) - exp(b[1])))
		# gen mod
		bias.b[ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1])) - exp(b[1])))
	}else{
		bias.b[ii,4] <- 888
		bias.b[ii,5] <- 888
		bias.b[ii,6] <- 888
	}
}
for (ii in 1:length(int50)){
	if (file.exists(paste(getwd(),'/output',int50[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',int50[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.b[n1+ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.b[n1+ii,4] <- mean((exp(temp) - exp(b[1])))
		# gen mod
		bias.b[n1+ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1]))- exp(b[1])))
	}else{
		bias.b[n1+ii,4] <- 888
		bias.b[n1+ii,5] <- 888
		bias.b[n1+ii,6] <- 888
	}
}


###########
#########################
#####################################


int20 <- 0:8*8+8
int50 <- 0:8*8+12
obs20 <- 0:8*8+6
obs50 <- 0:8*8+10

bias.o <- matrix(0,2*3*3,6)
bias.o[c(1,9),1] <- c(20,50)
colnames(bias.o) <- c('clusters','trt RR','cov RR','ind','eff','GEE')
bias.o[c(1,4,7,10,13,16),2] <- rep(c(1,1.25,2),2)
bias.o[,3] <- rep(c(1,1.25,2),6)

n1 <- length(obs20)
for (ii in 1:length(obs20)){
	if (file.exists(paste(getwd(),'/output',obs20[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',obs20[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.o[ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.o[ii,4] <- mean((exp(temp) - exp(b[1])))
		# gen mod
		bias.o[ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1])) - exp(b[1])))
	}else{
		bias.o[ii,4] <- 888
		bias.o[ii,5] <- 888
		bias.o[ii,6] <- 888
	}
}
for (ii in 1:length(obs50)){
	if (file.exists(paste(getwd(),'/output',obs50[ii],'.RData',sep=''))){
		load(paste(getwd(),'/output',obs50[ii],'.RData',sep=''))
		# gen mod
		temp <- unique(betas.eff[,1])
		temp <- temp[abs(temp)<20]
		bias.o[n1+ii,5] <- mean((exp(temp) - exp(b[1])))
		# independent
		temp <- unique(betas.uni[,1])
		temp <- temp[abs(temp)<20]
		bias.o[n1+ii,4] <- mean((exp(temp) - exp(b[1])))
		# gen mod
		bias.o[n1+ii,6] <- mean((unique(exp(matrix(as.numeric(betas.sas),ncol=2)[,1]))- exp(b[1])))
	}else{
		bias.o[n1+ii,4] <- 888
		bias.o[n1+ii,5] <- 888
		bias.o[n1+ii,6] <- 888
	}
}



