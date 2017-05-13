# solve score eqn
# generate data
# solve score eqn
# generate data
rm(list=ls())
library(mvtnorm)
library(binarySimCLF)
logit <- function(x){
	log(x/(1-x))
}
expit <- function(x){
	exp(x)/(1+exp(x))
}

gen.data <- function(n , k , q , b , a , R){
#n <- 500 # number of individuals
#k <- 3   # number of obs per cluster
#q <- 2	 # number of covariates

#b <- c(.3,-.5) # beta for q's
#a <- 0 # intercepts for individuals
#R <- correlation matrix
	dta <- NULL
	for (i in 1:n){
		dta <- rbind(dta , cbind(i,rmvnorm(k,rep(0,4))))
	}
	dta <- cbind(dta , rbinom(n*k,1,expit(dta[,2:5]%*%c(-1,-1,.3,1.3))))
	colnames(dta)[1] <- 'id'
	y <- numeric(n*k)
	for (i in 1:n){
		id <- which(dta[,'id']==i)
		mu <- as.vector(exp(dta[id,-1]%*%b+a))
		V = cor2var(R,mu)
		B = allReg(V);
		# Checks CLF compatibility.
		clf.compat = blrchk1(mu,B)
		if (clf.compat){
			y[id] = mbsclf(1,mu,B)$y
		} 
		else{ print("Not CLF compatible"); return("Not CLF compatible")}
	}
	dta <- cbind(dta,as.vector(y))
	colnames(dta)[q+2] <- 'Y'
	return(dta)
}

#################
score.uni <- function(Y,X,beta){
	exp.a <- t(Y)%*%exp(-X%*%beta)/length(Y)
	p.hat <- exp(X%*%beta)*as.numeric(exp.a)
	w.opt <- as.vector((1-p.hat)^{-1})*( X - matrix( (t(X)%*%(p.hat/(1-p.hat))) / as.numeric(t(p.hat)%*%(1-p.hat)) , dim(X)[1] , dim(X)[2] , byrow=T) )
	beta.opt <- beta + solve(t(w.opt)%*%(Y*X))%*%t(Y%*%w.opt)
	return(beta.opt)
}

score.u <- function(Y,X,beta){
	exp.a <- t(Y)%*%exp(-X*beta)/length(Y)
	p.hat <- exp(X*beta)*as.numeric(exp.a)
	w.opt <- as.vector((1-p.hat)^{-1})*( X - (X%*%(p.hat/(1-p.hat))) / as.numeric(t(p.hat)%*%(1-p.hat)) )
	beta.opt <- beta + solve(t(w.opt)%*%(Y*X))%*%t(Y%*%w.opt)
	return(beta.opt)
}

library(gregmisc)
score <- function(id,Y,X,beta,alpha,time=NULL){

	# some constants
	#n <- length(unique(id))
	#k <- length(Y)/n
	N <- length(Y)
	
	if (!is.null(time)){
		# allows for different intercepts
		exp.a <- Y*exp(-X%*%beta)	
		exp.a <- tapply(exp.a,time,mean)	
		mu    <- as.vector(exp.a[time])*exp(X%*%beta)
		eps   <- (Y-mu)
	}
	else{
		# single intercept
		if (alpha=='unspec'){
			exp.a <- t(Y)%*%exp(-X%*%beta)/(N)
		}else{ 
			if (alpha=='other'){
				exp.a <- Y*exp(-X%*%beta)
				exp.a <- tapply(exp.a,time,mean)
			} else{
				exp.a <- exp(alpha) }
		}
	
		mu    <- rep(exp.a,N)*exp(X%*%beta)
		eps   <- (Y-mu)
	}
	
	#check that all mu's are between 0 and 1...if not smooth!!
	if (any(mu>1)){
		mu.t <- mu
		mu <- glm(Y~mu.t+mu.t^2+mu.t^3,family='binomial')$fitted.values
	}
	
	# estimate rho from exchangeble correlation
	if ( N==length(unique(id))){
		rho <- 1
	}else{
		# create all pairwise combos of Y
		Y.t <- NULL
		for (i in id){
			Y.id <- Y[id==i]
			if (length(Y.id)>1){
				MM <- combinations(n=length(Y.id),r=2)
				Y.t <- rbind(Y.t , matrix(Y.id[MM],ncol=2))
			}
		}
		rho <- cor(Y.t)[1,2]
		#Y.t <- matrix(Y,n,byrow=T)
		#rho <- sum(cor(Y.t)-diag(rep(1,k)))/(k*(k-1))
	}
	# exchangeable correlation for each size of cluster
	R <- list()
	for (k in 1:max(tapply(rep(1,length(id)),id,sum))){
		R[[k]] <- matrix(rho,k,k)
		for (i in 1:k){
			R[[k]][i,i] <- 1
		}
	}
	
	
	
	# calculate working vcov
	V <- list()
	V.inv <- list()
	for (i in id){
		mu.V <- mu[id==i]
		k <- length(mu.V)
		V[[i]] <- diag(sqrt(mu.V*(1-mu.V)),nrow=k)%*%R[[k]]%*%diag(sqrt(mu.V*(1-mu.V)),nrow=k)
		V.inv[[i]] <- solve(V[[i]])
	}
	
	#V 	  <- tapply(eps,id,function(x){x%*%t(x)})
	#V 	  <- Reduce('+',V)/n
	#V <- diag(sum(eps^2)/n,k)
	#V.inv <- solve(V)
	MV <- list()
	MVM <- list()
	MVMX <- list()
	XMVMX <- list()
	MVY <- list()
	XMVY <- list()
	for (i in id){
		MV[[i]]    <- diag(mu[id==i],nrow=length(mu[id==i]))%*%V.inv[[i]]
		MVM[[i]]   <- MV[[i]]%*%diag(mu[id==i],nrow=length(mu[id==i]))
		MVMX[[i]]  <- MVM[[i]]%*%X[id==i,]
		if (sum(id==i)==1){
			XMVMX[[i]] <- X[id==i,]%*%MVMX[[i]]
			MVY[[i]]  <-  MV[[i]]%*%Y[id==i]
			XMVY[[i]]  <- X[id==i,]%*%MVY[[i]]
		}else{
			XMVMX[[i]] <- t(X[id==i,])%*%MVMX[[i]]
			MVY[[i]]  <-  MV[[i]]%*%Y[id==i]
			XMVY[[i]]  <- t(X[id==i,])%*%MVY[[i]]
		}
	}
	# semi-parametric variance bound
	k <- tapply(rep(1,length(id)),id,sum)
	
	if (!is.null(time)){
		VB1 <- Reduce('+',XMVMX)
		VB2 <- list()
		VB3 <- list()
		#U1 <- Reduce('+',XMVY)
		U1 <- list()
		U2 <- list()
		U3 <- list()
		for (j in k){
			VB2[[j]] <- 0
			VB3[[j]] <- 0
			U1[[j]] <- 0
			U2[[j]] <- 0
			U3[[j]] <- 0
		}
		for (i in id){
			VB2[[k[i]]] <- VB2[[k[i]]] + t(MVMX[[i]])%*%rep(1,k[[i]])/k[[i]]
			VB3[[k[i]]] <- VB3[[k[i]]] + sum(MVM[[i]])/k[[i]]
			#VB <- Reduce('+',XMVMX)-t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVMX)/sum(Reduce('+',MVM))
			# sum U_i
			#U <- Reduce('+',XMVY) - t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVY)/sum(Reduce('+',MVM))
			U1[[k[i]]] <- U1[[k[i]]] + XMVY[[i]]/k[[i]]
			U2[[k[i]]] <- U2[[k[i]]] + t(MVMX[[i]])%*%rep(1,k[[i]])/k[[i]]
			U3[[k[i]]] <- U3[[k[i]]] + t(rep(1,k[[i]]))%*%MVY[[i]]/k[[i]]
		}
		VB <- 0
		U <- 0
		for (j in unique(k)){
			VB <- VB + (VB1 - VB2[[j]]%*%t(VB2[[j]])/VB3[[j]])*sum(k==j)/length(unique(id))
			U <- U + (U1[[j]] - U2[[j]]%*%U3[[j]]/VB3[[j]])*sum(k==j)/length(unique(id))
		}
	}
	else{
		VB1 <- Reduce('+',XMVMX)
		VB2 <- 0
		VB3 <- 0
		U1 <- Reduce('+',XMVY)
		U2 <- 0
		U3 <- 0
		for (i in id){
			VB2 <- VB2 + t(MVMX[[i]])%*%rep(1,k[[i]])/k[[i]]
			VB3 <- VB3 + sum(MVM[[i]])/k[[i]]
			#VB <- Reduce('+',XMVMX)-t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVMX)/sum(Reduce('+',MVM))
			# sum U_i
			#U <- Reduce('+',XMVY) - t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVY)/sum(Reduce('+',MVM))
			U2 <- U2 + t(MVMX[[i]])%*%rep(1,k[[i]])/k[[i]]
			U3 <- U3 + t(rep(1,k[[i]]))%*%MVY[[i]]/k[[i]]
		}
		VB <- VB1 - VB2%*%t(VB2)/VB3
		U  <- U1 - U2%*%U3/VB3
	}
	# updated beta
	bb <- beta + solve(VB)%*%U
	
	# conf int
	ci <- cbind(bb,bb) + t(c(-1.96,1.96)%*%t(sqrt(diag(solve(VB)))))
	
	return(list(beta=bb,ci=ci))
}

my.estimator <- function(Y,X,ID,time=NULL){
	X <- as.matrix(X)
	# initial fit 
	if (is.null(time)){
		cox.model <- coxph(Surv(rep(1,length(Y)),Y)~X)
	}
	else{
		cox.model <- coxph(Surv(rep(1,length(Y)),Y)~X+strata(time))
	}
	# excludes any NA is outcome
	Y.a  <- Y[!is.na(Y)]
	X.a  <- as.matrix(X[!is.na(Y),])
	ID.a <- ID[!is.na(Y)]
	t.a <- time[!is.na(Y)]
	# assumes independence
	model.uni <- score(id=1:length(Y.a),Y=Y.a,X=as.matrix(X.a),beta=cox.model$coef,alpha='unspec',time=t.a)
	# exchangeable correlation structure
	model.eff <- score(id=ID.a,Y=Y.a,X=X.a,beta=model.uni$beta,alpha='unspec',time=t.a)
	# more updates
	model2 <- score(id=ID.a,Y=Y.a,X=X.a,beta=model.eff$beta,alpha='unspec',time=t.a)
	model3 <- score(id=ID.a,Y=Y.a,X=X.a,beta=model2$beta,alpha='unspec',time=t.a)
	model4 <- score(id=ID.a,Y=Y.a,X=X.a,beta=model3$beta,alpha='unspec',time=t.a)
	return(list(efficient=model.eff , independent=model.uni , step2=model2 , step3=model3 ,step4=model4))
}
library(survival)

###########################
###########################

library(gee)


num <- Sys.getenv('LSB_JOBINDEX')
num = as.numeric(num)
set.seed(num)

n <- num
k <- 5
q <- 5
b <- c(.1,-.1,.05,-.05,0)
a <- -1.5
M <- 10000

# Unstructured correlation matrix.
R = matrix(0.3 , k ,k)
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

library(survival)

for (i in 1:M){
	dta <- 'temp'
	count <- 0
	while (class(dta)=="character"){
		dta <- gen.data(n=n , k=k , q=q , b=b , a=a , R=R)
		count <- count+1
		if (count>1000){
			print(dta)
			break
		}
	}
	# initial guess
	#t.c <- rep(1,dim(dta)[1])
	#event <- dta[,'Y']
	#cox.model <- coxph(Surv(t.c,event)~dta[,2:(1+q)])
	#model.glm <- glm(dta[,'Y']~dta[,2:(1+q)] , family=binomial(link='log'))
	
	# update for independent observations
	#model.uni <- score(id=1:dim(dta)[1],Y=dta[,'Y'],X=dta[,2:(1+q)],beta=cox.model$coef[1:(q)],alpha='unspec')
	#model.uni <- score.uni(Y=dta[,'Y'],X=dta[,2:(1+q)],beta=model.glm$coef[2:(q+1)])
	#model.u1  <- score.u(Y=Y,X=X[,1],beta=model.glm$coef[2])
	#model.gee <- gee(dta[,'Y']~dta[,2:(1+q)], id=dta[,1], R = NULL, b = NULL, tol = 0.001, maxiter = 25, family = 'binomial', corstr = "independence", silent=T)
	
	
	# update for our estimator
	#model.eff <- score(id=1:600,Y=dta[,'Y'],X=dta[,2:(1+q)],beta=model.glm$coef[2:(q+1)],alpha=exp(model.glm$coef[1]))
	#model.eff <- score(id=dta[,1],Y=dta[,'Y'],X=dta[,2:(1+q)],beta=model.uni$beta,alpha='unspec')
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


coverage <- numeric(q)
for (j in 1:q){
	coverage[j] <- sum((ci.eff[[j]][,1]<b[j]) & (ci.eff[[j]][,2]>b[j]))/M
}

coverage.sas <- numeric(q)
for (j in 1:q){
	coverage.sas[j] <- sum((ci.sas[[j]][,1]<b[j]) & (ci.sas[[j]][,2]>b[j]))/M
}

save.image(paste('./results/output',num,'.RData',sep=''))
