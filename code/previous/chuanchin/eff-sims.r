# solve score eqn
# generate data
# solve score eqn
# generate data
rm(list=ls())

source('./code/score.r')
source('./code/my_estimator.r')
source('./code/gen_data_my.r')
source('./code/gen_data_other.r')
source('./code/variance.r')

###########################
###########################

load('./input/input.RData')
num <- Sys.getenv('LSB_JOBINDEX')
num = as.numeric(num)
set.seed(num)

n <- input[[num]][['n']]
k <- input[[num]][['k']] #5 -- number per cluster
q <- input[[num]][['q']] #5 -- number of covariates
b <- input[[num]][['b']] #c(.1,-.1,.05,-.05,0) -- betas 
a <- input[[num]][['a']] #-1.5 -- intercept
M <- input[[num]][['M']] #10000 -- replications
rho <- input[[num]][['rho']]
who <- input[[num]][['who']] # my sims or yelland sims
type <- input[[num]][['type']] # obs or intervention
cont <- input[[num]][['cont']] # T or F
M <- 1000
#n <- 1000

# Unstructured correlation matrix.
R = matrix(rho , k ,k)
for (i in 1:k){
	R[i,i] = 1
}

#R <- diag(5)
#R <- R + (1-diag(5))*.1
#R[1,2] <- R[2,1] <- .05
#R[3,4] <- R[3,5] <- .2
#R[4,3] <- R[4,5] <- .2
#R[5,3] <- R[5,4] <- .2


betas.uni  <- matrix(0,M,q)
betas.eff  <- matrix(0,M,q)
betas.sas  <- matrix(0,M,q)
betas.pois <- matrix(0,M,q)
betas.pois.eff <- matrix(0,M,q)
ci.uni <- list()
ci.eff <- list()
ci.sas <- list()
ci.pois <- list()
ci.pois.eff <- list()
for (i in 1:q){
	ci.uni[[i]] <- matrix(0,M,2)
	ci.eff[[i]] <- matrix(0,M,2)
	ci.sas[[i]] <- matrix(0,M,2)
	ci.pois[[i]] <- matrix(0,M,2)
	ci.pois.eff[[i]] <- matrix(0,M,2)
}


if (who!='me'){
	if ((type!='obs')&(cont)){
		p0 <- .2
		rrtrt <- exp(b[1])
		varre <- .2
		rrcov <- exp(b[2])
		varnc <- .25
		rhonc <- .05
		meannc <- .5
		params <- matrix(c(M,n,k,p0,rrtrt,varre,rrcov,varnc,rhonc,meannc),nrow=1)
		colnames(params) <- c('M','n','k','p0','rrtrt','varre','rrcov','varnc','rhonc','meannc')
		write.csv(params,file=paste('./temp/params',num,'.csv',sep=''),row.names=F,quote=F)
		
		system(paste("sas -nolog ./code/fwdajepaper/get_params_int_cont.sas -sysparm ", num ," -print ./out/temp-",num,".lst" ,sep=''),wait=TRUE)
	}
	if ((type!='obs')&(!cont)){
		p0 <- .2
		rrtrt <- exp(b[1])
		varre <- .2
		rrcov <- exp(b[2])
		pbc <- .5
		rhocov <- .05
		params <- matrix(c(M,n,k,p0,rrtrt,varre,rrcov,pbc,rhocov),nrow=1)
		colnames(params) <- c('M','n','k','p0','rrtrt','varre','rrcov','pbc','rhocov')
		write.csv(params,file=paste('./temp/params',num,'.csv',sep=''),row.names=F,quote=F)
		
		system(paste("sas -nolog ./code/fwdajepaper/get_params_int_binary.sas -sysparm ", num ," -print ./out/temp-", num ,".lst"  ,sep=''),wait=TRUE)
	}
	if ((type=='obs')&(cont)){
		p0 <- .2
		rrtrt <- exp(b[1])
		varre <- .2
		rrcov <- exp(b[2])
		varcov <- .25
		rhocov <- .05
		meanunexp <- .4
		meanexp <- .6
		params <- matrix(c(M,n,k,p0,rrtrt,varre,rrcov,varcov,rhocov,meanunexp,meanexp),nrow=1)
		colnames(params) <- c('M','n','k','p0','rrtrt','varre','rrcov','varcov','rhocov','meanunexp','meanexp')
		write.csv(params,file=paste('./temp/params',num,'.csv',sep=''),row.names=F,quote=F)

		system(paste("sas -nolog ./code/fwdajepaper/get_params_obs_cont.sas -sysparm ", num ," -print ./out/temp-",num,".lst"  ,sep=''),wait=TRUE)
	}
	if ((type=='obs')&(!cont)){
		p0 <- .2
		rrtrt <- exp(b[1])
		varre <- .2
		pbcu <- .4
		pbce <- .6
		rhocov <- .05
		rrcov <- exp(b[2])	
		params <- matrix(c(M,n,k,p0,rrtrt,varre,pbcu,pbce,rhocov,rrcov),nrow=1)
		colnames(params) <- c('M','n','k','p0','rrtrt','varre','pbcu','pbce','rhocov','rrcov')
		write.csv(params,file=paste('./temp/params',num,'.csv',sep=''),row.names=F,quote=F)
		
		system(paste("sas -nolog ./code/fwdajepaper/get_params_obs_binary.sas -sysparm ", num ," -print ./out/temp-",num,".lst"  ,sep=''),wait=TRUE)
	}


	fullDTA <-  as.matrix(read.csv(paste("./temp/fullDTA",num,'.csv',sep='')))
	colnames(fullDTA)[9] <- 'Y'
	colnames(fullDTA)[2] <- 'id'
}


for (i in 1:M){
	print(i)
	dta <- 'temp'
	count <- 0
	while (class(dta)=="character"){
		if (who=='me'){
			dta <- gen.data.my(n=n , k=k , q=q , b=b , a=a , R=R)
		}
		else{ 
			#dta <- gen.data.other(n=n , k=k , b.trt=b[1] , b.cov=b[2] , b0=.1 , rho.cov=0.05 , R=R , type=type , continuous=cont) 
			#system()
			dta <- fullDTA[fullDTA[,'simulation']==i,]	
			dta <- dta[,c(2,6,8,9)]
		}
		count <- count+1
		if (count>1000){
			print(dta)
			break
		}
	}

	out <- my.estimator(Y=dta[,'Y'], X=dta[,2:(1+q)] ,ID=dta[,1])
	
	betas.uni[i,] <- out$independent$beta
	betas.eff[i,] <- out$efficient$beta
	for (j in 1:q){
		ci.uni[[j]][i,]    <- out$independent$ci[j,]
		ci.eff[[j]][i,]    <- out$efficient$ci[j,]
	}
	
	colnames(dta)[2:(1+q)] <- paste('c',1:q,sep='') 	
	write.csv(dta,file=paste('./temp/dta',num,'.csv',sep=''),row.names=F)
	
	if (who=='me'){
		system(paste("sas -nolog ./code/genmod.sas -sysparm ", num, " -print ./out/genmod-",num,".lst" , sep=''),wait=TRUE)
	}
	else{
		system(paste("sas -nolog ./code/genmod-other.sas -sysparm ", num ," -print ./out/genmod-other-",num,".lst",sep=''),wait=TRUE)
	}

	if (file.exists(paste("./temp/sasout",num,'.csv',sep=''))){	
		sas <- as.matrix(read.csv(paste("./temp/sasout",num,'.csv',sep='')))
		betas.sas[i,] <- sas[2:(1+q),2]
        	for (j in 1:q){
                	ci.sas[[j]][i,]    <- sas[(j+1),4:5]
        	}
		system(paste("rm ./temp/sasout", num , '.csv',sep=''),wait=TRUE)
			
	}
	if (file.exists(paste("./temp/poisson",num,'.csv',sep=''))){	
		sas <- as.matrix(read.csv(paste("./temp/poisson",num,'.csv',sep='')))
		betas.pois[i,] <- sas[2:(1+q),2]
        	for (j in 1:q){
                	ci.pois[[j]][i,]    <- sas[(j+1),4:5]
        	}	
		system(paste("rm ./temp/poisson", num , '.csv',sep=''),wait=TRUE)
	}
	system(paste('rm ./temp/dta',num,'.csv',sep=''),wait=TRUE)
	
	# use poisson regression as the starting point	
	out <- my.estimator(Y=dta[,'Y'], X=dta[,2:(1+q)] ,ID=dta[,1] , beta.pois=as.numeric(betas.sas[i,]))
	betas.pois.eff[i,] <- out$efficient$beta
	for (j in 1:q){
		ci.pois.eff[[j]][i,]    <- out$efficient$ci[j,]
	}
}
warnings()
#apply(betas.uni,2,mean)
#apply(betas.eff,2,mean)

#apply(betas.uni,2,var)
#apply(betas.eff,2,var)

coverage <- numeric(q)
for (j in 1:q){
	coverage[j] <- sum((ci.eff[[j]][,1]<b[j]) & (ci.eff[[j]][,2]>b[j]))/M
}

coverage.u <- numeric(q)
for (j in 1:q){
	coverage.u[j] <- sum((ci.uni[[j]][,1]<b[j]) & (ci.uni[[j]][,2]>b[j]))/M
}

coverage.sas <- numeric(q)
for (j in 1:q){
	coverage.sas[j] <- sum((ci.sas[[j]][,1]<b[j]) & (ci.sas[[j]][,2]>b[j]))/M
}

save.image(paste('./results/output',num,'.RData',sep=''))
