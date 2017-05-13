# solve score eqn
# generate data
# solve score eqn
# generate data
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric')

source('./code/score.r')
source('./code/my_estimator.r')
source('./code/gen_data_my.r')
source('./code/gen_data_other.r')


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

# Unstructured correlation matrix.
R = matrix(rho , k ,k)
#r2 = c(0.1, 1, 0.3)
#r3 = c(0.2, 0.3, 1)
for (i in 1:k){
	R[i,i] = 1
}

betas.uni <- matrix(0,M,q)
betas.eff <- matrix(0,M,q)
betas.sas <- matrix(0,M,q)
ci.uni <- list()
ci.eff <- list()
ci.sas <- list()
for (i in 1:q){
	ci.uni[[i]] <- matrix(0,M,2)
	ci.eff[[i]] <- matrix(0,M,2)
	ci.sas[[i]] <- matrix(0,M,2)
}

for (i in 1:M){
	dta <- 'temp'
	count <- 0
	while (class(dta)=="character"){
		if (who=='me'){
			dta <- gen.data.my(n=n , k=k , q=q , b=b , a=a , R=R)
		}
		else{ dta <- gen.data.other(n=n , k=k , b.trt=b[1] , b.cov=b[2] , b0=.1 , rho.cov=0.05 , R=R , type=type , continuous=cont) }
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
	write.csv(dta,file=paste('./results/dta',num,'.csv',sep=''),row.names=F)
	system(paste("sas genmod",num,'.sas',sep=''))
	
	if (file.exists(paste("./results/sasout",num,'.csv',sep=''))){	
		sas <- as.matrix(read.csv(paste("./results/sasout",num,'.csv',sep='')))
		betas.sas[i,] <- sas[2:6,2]
        	for (j in 1:q){
                	ci.sas[[j]][i,]    <- sas[(j+1),4:5]
        	}	
	}
}

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
