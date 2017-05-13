
# function to find updated estimator
update.RR <- function(id,Y,X,beta,alpha,time=NULL,smooth='bspline',corr='exchangeable'){
	# some constants
	N <- length(Y)
	if (!is.null(time)){
		# allows for different intercepts
		exp.a <- Y*exp(-X%*%beta)	
		exp.a <- tapply(exp.a,time,mean)	
		mu    <- as.vector(exp.a[time])*exp(X%*%beta)
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
	}
	mu.true <- mu	
	#check that all mu's are between 0 and 1...if not smooth!!
	if (any(mu>1)){
		# just use 3rd order polynomial
		mu.t <- X%*%beta
		mu.t2 <- mu.t^2
		mu.t3 <- mu.t^3
		mu <- glm(Y~mu.t+mu.t2+mu.t3,family='binomial')$fitted.values
	}
	# calculate dispersion parameter
	phi <- max(1,sum(((Y-mu)/sqrt(mu*(1-mu)))^2)/(length(Y)-length(beta)))
	# estimate rho from exchangeble correlation or AR-1
	if ( N==length(unique(id))){
		rho <- 1
	}else{
      if (corr=='exchangeable'){
	   	eps   <- (Y-mu)/sqrt(mu*(1-mu))
		   rho <- 0
		   N.star <- 0
		   for (i in unique(id)){
			   id.temp <- which(id==i)
			   if ( length(id.temp)>1 ){
   				for (k in 2:length(id.temp)){
   					for (j in 1:(k-1)){
   						rho <- rho + eps[id.temp[j]]*eps[id.temp[k]]	
   					}
   				}
   			}
		   	N.star <- N.star + .5*length(id.temp)*(length(id.temp)-1)
	   	}
	   	rho <- rho / (N.star-length(beta))/phi
      }
      if (corr=='AR-1'){
         eps   <- (Y-mu)/sqrt(mu*(1-mu))
         rho <- 0
         N.star <- 0
         for (i in unique(id)){
            id.temp <- which(id==i)
            if ( length(id.temp)>1 ){
               #for (k in 2:length(id.temp)){
               for (j in 1:(length(id.temp)-1)){
                  rho <- rho + eps[id.temp[j]]*eps[id.temp[j+1]]	
               }
               #}
            }
            N.star <- N.star + (length(id.temp)-1)
         }
         rho <- rho / (N.star-length(beta))/phi
      }
	}
	#  correlation for each size of cluster
   if (corr=='exchangeable'){
   	R <- list()
   	for (k in 1:max(tapply(rep(1,length(id)),id,sum))){
   		R[[k]] <- matrix(rho,k,k)
   		for (i in 1:k){
   			R[[k]][i,i] <- 1
   		}
   	}
   }
	if (corr=='AR-1'){
	   R <- list()
	   for (k in 1:max(tapply(rep(1,length(id)),id,sum))){
	      R[[k]] <- abs(outer(1:k,1:k,'-'))
         R[[k]] <- rho^R[[k]]
	   }
	}
   
	# calculate working vcov
	V <- list()
	V.inv <- list()
	for (i in unique(id)){
		mu.V <- mu[id==i]
		k <- length(mu.V)
		V[[i]] <- diag(sqrt(mu.V*(1-mu.V)),nrow=k)%*%R[[k]]%*%diag(sqrt(mu.V*(1-mu.V)),nrow=k)
		V.inv[[i]] <- solve(V[[i]]) 
	}
	# reset mean vector
	mu <- exp(X%*%beta)*as.numeric(exp.a)
	# Y is now epsilon
	Y <- Y-mu
	# calculate some stuff thats needed later
	MV <- list()
	MVM <- list()
	MVMX <- list()
	XMVMX <- list()
	MVY <- list()
	XMVY <- list()
	for (i in unique(id)){
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
	VB1 <- Reduce('+',XMVMX)
	VB2 <- 0
	VB3 <- 0
	U1 <- Reduce('+',XMVY)
	U2 <- 0
	U3 <- 0
	for (i in unique(id)){
		VB3 <- VB3 + sum(MVM[[i]]) 
		U2 <- U2 + t(MVMX[[i]])%*%rep(1,k[[i]]) 
		U3 <- U3 + t(rep(1,k[[i]]))%*%MVY[[i]]  
	}
	VB <- VB1 - U2%*%t(U2)/VB3
	U  <- U1 - U2%*%U3/VB3
	# updated beta
	bb <- beta + solve(VB)%*%U
	# reset Y
	Y <- Y + mu
	# find variance estimate
	sigma.hat <- var.eff(id=id,Y=Y,X=X,beta=bb,exp.a=as.numeric(exp.a),time=NULL,phi=phi,V.inv=V.inv,EXMVm=U2,EmVm=VB3)$sigma 	
	# conf int
	ci <- cbind(bb,bb) + t(c(-1.96,1.96)%*%t(sqrt(diag(sigma.hat))))
	# test stats
	z <- bb/sqrt(diag(sigma.hat))
	p.value <- 2*(1-pnorm(abs(z)))
	return(list(beta=bb,ci=ci,sigma.hat=sigma.hat,test.stat=z,p.value=p.value))
}

# function to calculate variance estimate based on sandwich estimator
var.eff <- function(id,Y,X,beta,exp.a,time=NULL,phi,V.inv,EXMVm,EmVm){
	# some constants
	N <- length(Y)
	if (!is.null(time)){
		# allows for different intercepts
		exp.a <- Y*exp(-X%*%beta)	
		exp.a <- tapply(exp.a,time,mean)	
		mu    <- as.vector(exp.a[time])*exp(X%*%beta)
	}
	else{
		mu    <- rep(exp.a,N)*exp(X%*%beta)
	}
	# Y is now epsilon
	Y <- Y-mu	
	# calculate various terms that are needed
	MV <- list()
	MVM <- list()
	MVMX <- list()
	XMVMX <- list()
	MVY <- list()
	XMVY <- list()
	for (i in unique(id)){
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
	VB1 <- Reduce('+',XMVMX)
	VB2 <- 0
	# empirical estimate of variance
	I1 <- matrix(0 , length(beta) , length(beta)) 
	for (i in unique(id)){
		VB2 <- VB2 + t(MVMX[[i]])%*%rep(1,k[[i]])
		temp <- XMVY[[i]]-EXMVm%*%(t(rep(1,k[[i]]))%*%MVY[[i]])/EmVm
		I1 <- I1 + temp%*%t(temp) 
	}
	VB <- VB1 - EXMVm%*%t(VB2)/EmVm
	# sandwich variance estimate
	sigma.hat <- phi*solve(VB)%*%I1%*%t(solve(VB)) 	
	return(list(sigma=sigma.hat))
}

binary.gee <- function(Y,X,ID,time=NULL,beta.init=NULL,corr='exchangeable'){
		ID.t <- numeric(length(ID))
		i <- 1
		for (id in unique(ID)){
			ID.t[ID==id] <- i
			i <- i+1
		}
		X <- as.matrix(X)
		if (is.null(beta.init)){
			library(gee)
			# modified poisson as initial estimate!
         #if (corr=='AR-1'){
         #   capture.output(beta.init <- gee(Y~X,family='poisson',corstr='AR-M',Mv=1,id=ID.t)$coef[-1])
         #}else{
			   capture.output(beta.init <- gee(Y~X,family='poisson',corstr='exchangeable',id=ID.t)$coef[-1])
         #}
      }
		model.eff <- update.RR(id=ID.t,Y=Y,X=X,beta=beta.init,alpha='unspec',time=time,corr=corr)
		class(model.eff) <- 'binary.gee'
		print(model.eff)
		return(model.eff)
}

`print.binary.gee` <-
function(x, ...){
		# print results to screen
		cat("\nCoefficients:\n")
		out <- cbind(Estimate=as.numeric(format(as.vector(x$beta),digits=5)) , StdError=as.numeric(format(sqrt(diag(x$sigma.hat)),digits=5)) , TestStat=as.numeric(format(as.vector(x$test.stat),digits=3)), Pvalue=as.numeric(format(as.vector(x$p.value),scientific=F)))
		rownames(out) <- rownames(x$beta)
		print(out)
}

###########################
####### Example ###########

# this function requires the library gee
# so install it if it is not already in library
# library(gee)
# install.packages('gee')

# Generate some binary data
x <- matrix(rnorm(200),100,2) # the covariates
y <- rbinom(100,1,.5) # the outcome -- the model is Y~X
id <- rep(1:50,2) # id for clusters -- needs to be indexed from 1 to # of clusters

model <- binary.gee(Y=y,X=x,ID=id)

#############################
# what is saved in the model?

# the estimated RR
model$beta
# the std errors
sqrt(diag(model$sigma.hat))
# the conf int
model$ci
# the test stats
model$test.stat
# the p-values
model$p.value


# Generate random cluster sizes
k <- 1+rpois(50,5)
id <- paste('A',rep(1:50,k),sep='')

# generate the data
x <- matrix(rnorm(2*length(id)),length(id),2) # the covariates
y <- rbinom(length(id),1,.5) # the outcome -- the model is Y~X

model <- binary.gee(Y=y,X=x,ID=id)


#############################
# what is saved in the model?

# the estimated RR
model$beta
# the std errors
sqrt(diag(model$sigma.hat))
# the conf int
model$ci
# the test stats
model$test.stat
# the p-values
model$p.value



################
## likelihood ##
library(mvtnorm)
log.L <- function(Y,X,beta,rho,alpha,ID){
   mu <- exp(x%*%t(beta)+alpha)
   if (all(mu<1)){
      R <- (Y-mu)/sqrt(mu*(1-mu))
      l.star <- Y%*%log(mu)+(1-Y)%*%log(1-mu)
      b.star <- 0
      for (id in ID){
         temp <- R[ID==id]
         temp <- sapply(temp,function(x) rho*x*temp)
         for (i in 1:(dim(temp)[1]-1)){
            temp1 <- temp[i:dim(temp)[1],i:dim(temp)[1]]
            if (sum(temp1[lower.tri(temp1)])<(-1)){
               return(-10000000)
            }else{
               b.star <- b.star + log(1+sum(temp1[lower.tri(temp1)]))
            }
         }
      }
      return(l.star+b.star)
   }
   else{
      return(-10000000)
   }
}

beta. = model$beta
beta.var = model$sigma.hat
beta.var[1,2]=beta.var[2,1]

i=1
beta = matrix(0,1000,2)
rho = numeric(1000)
alpha = numeric(1000)
while (i <= 1000){
   beta.prop = rmvnorm(1,beta.,beta.var)
   rho.prop = runif(1,-1,1)
   alpha.prop = rnorm(1)
   num = log.L(y,x,beta.prop,rho.prop,alpha.prop,ID=id)
   den = log.L(y,x,beta.prop,rho.prop,alpha.prop,ID=id)
   if (exp(num-den)>runif(1)){
      beta[i,] = beta.prop
      alpha[i] = alpha.prop
      rho[i] = rho.prop
      i = i+1
   }
}





#########################

x <- read.table('x.txt')
y <- read.table('y.txt')[[1]]
id <- read.table('id.txt')[[1]]

scale.x <- scale(x[,c(2,3,4,15,16)],scale=F)

x[,c(2,3,4,15,16)] <- scale.x

model <- binary.gee(Y=y,X=x,ID=id)
