
var.eff <- function(id,Y,X,beta,exp.a,time=NULL,phi,V.inv,EXMVm,EmVm){

	# some constants
	#n <- length(unique(id))
	#k <- length(Y)/n
	N <- length(Y)
	
	if (!is.null(time)){
		# allows for different intercepts
		exp.a <- Y*exp(-X%*%beta)	
		exp.a <- tapply(exp.a,time,mean)	
		mu    <- as.vector(exp.a[time])*exp(X%*%beta)
		#eps   <- (Y-mu)/mu/(1-mu)
	}
	else{
		# single intercept
		#if (alpha=='unspec'){
		# exp.a <- t(Y)%*%exp(-X%*%beta)/(N)
		#}else{ 
		#	if (alpha=='other'){
		#		exp.a <- Y*exp(-X%*%beta)
		#		exp.a <- tapply(exp.a,time,mean)
		#	} else{
		#		exp.a <- exp(alpha) }
		#}
	
		mu    <- rep(exp.a,N)*exp(X%*%beta)
		#eps   <- (Y-mu)/mu/(1-mu)
	}
	
	Y <- Y-mu	
	
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
		#U1 <- Reduce('+',XMVY)
		#U2 <- 0
		#U3 <- 0
		# empirical estimate of variance
		I1 <- matrix(0 , length(beta) , length(beta)) 
		for (i in unique(id)){
			#VB1 <- VB1 + XMVMX[[i]]/k[[i]]
			VB2 <- VB2 + t(MVMX[[i]])%*%rep(1,k[[i]]) #/k[[i]]
			#VB3 <- VB3 + sum(MVM[[i]]) #/k[[i]]
			#VB <- Reduce('+',XMVMX)-t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVMX)/sum(Reduce('+',MVM))
			# sum U_i
			#U <- Reduce('+',XMVY) - t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVY)/sum(Reduce('+',MVM))
			#U2 <- U2 + t(MVMX[[i]])%*%rep(1,k[[i]]) # /k[[i]]
			#U3 <- U3 + t(rep(1,k[[i]]))%*%MVY[[i]]  # /k[[i]]
			
			temp <- XMVY[[i]]-EXMVm%*%(t(rep(1,k[[i]]))%*%MVY[[i]])/EmVm
			I1 <- I1 + temp%*%t(temp) 
		}
		VB <- VB1 - EXMVm%*%t(VB2)/EmVm
		#VB <- VB1 - VB2%*%t(VB2)/VB3
		#U  <- U1 - U2%*%U3/VB3
		# empirical estimate of variance
		#I1 <- matrix(0 , length(beta) , length(beta)) 
		#for (i in unique(id)){
		#	temp <- XMVY[[i]]-VB2%*%(t(rep(1,k[[i]]))%*%MVY[[i]])/VB3
		#	I1 <- I1 + temp%*%t(temp) 
		#}

	# sandwich variance estimate
	sigma.hat <- phi*solve(VB)%*%I1%*%solve(VB) 	
	print(phi)
	return(list(sigma=sigma.hat))
}
