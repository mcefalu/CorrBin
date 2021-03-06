
library(gregmisc)
library(gee)
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
	}else{
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
		#mu.t <- exp(X%*%beta)
		#mu.t <- mu
		#mu.t2 <- mu.t^2
		#mu.t3 <- mu.t^3
		#mut.t4 <- mu.t^4
		#mut.t5 <- mu.t^5
		#mu <- glm(Y~mu.t+mu.t2+mu.t3 ,family='binomial')$fitted.values
		#mu <- gee(Y~mu.t+mu.t2+mu.t3,family='binomial',id=id,corstr='exchangeable')$fitted.values
		#mu[mu>1] <- .5 # mu.t[mu>1]
		mu <- mu/max(1,mu)/2
	}
	
	# estimate rho from exchangeble correlation
	if ( N==length(unique(id))){
		rho <- 1
	}else{
		# create all pairwise combos of Y
		Y.t <- NULL
		for (i in unique(id)){
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
	for (i in unique(id)){
		mu.V <- mu[id==i]
		k <- length(mu.V)
		V[[i]] <- diag(sqrt(mu.V*(1-mu.V)),nrow=k)%*%R[[k]]%*%diag(sqrt(mu.V*(1-mu.V)),nrow=k)
		V.inv[[i]] <- solve(V[[i]])
	}
	# reset mean vector
	#mu <- mu.t
	
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
	
	#if (!is.null(time)){
		#VB1 <- Reduce('+',XMVMX)
		#VB2 <- list()
		#VB3 <- list()
		#U1 <- Reduce('+',XMVY)
		#U1 <- list()
		#U2 <- list()
		#U3 <- list()
		#for (j in k){
		#	VB2[[j]] <- 0
		#	VB3[[j]] <- 0
		#	U1[[j]] <- 0
		#	U2[[j]] <- 0
		#	U3[[j]] <- 0
		#}
		#for (i in id){
		#	VB2[[k[i]]] <- VB2[[k[i]]] + t(MVMX[[i]])%*%rep(1,k[[i]])/k[[i]]
		#	VB3[[k[i]]] <- VB3[[k[i]]] + sum(MVM[[i]])/k[[i]]
		#	#VB <- Reduce('+',XMVMX)-t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVMX)/sum(Reduce('+',MVM))
		#	# sum U_i
		#	#U <- Reduce('+',XMVY) - t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVY)/sum(Reduce('+',MVM))
		#	U1[[k[i]]] <- U1[[k[i]]] + XMVY[[i]]/k[[i]]
		#	U2[[k[i]]] <- U2[[k[i]]] + t(MVMX[[i]])%*%rep(1,k[[i]])/k[[i]]
		#	U3[[k[i]]] <- U3[[k[i]]] + t(rep(1,k[[i]]))%*%MVY[[i]] /k[[i]]
		#}
		#VB <- 0
		#U <- 0
		#for (j in unique(k)){
		#	VB <- VB + (VB1 - VB2[[j]]%*%t(VB2[[j]])/VB3[[j]])*sum(k==j)/length(unique(id))
		#	U <- U + (U1[[j]] - U2[[j]]%*%U3[[j]]/VB3[[j]])*sum(k==j)/length(unique(id))
		#}
	#}
	#else{
		VB1 <- Reduce('+',XMVMX)
		VB2 <- 0
		VB3 <- 0
		U1 <- Reduce('+',XMVY)
		U2 <- 0
		U3 <- 0
		for (i in unique(id)){
			#VB1 <- VB1 + XMVMX[[i]]/k[[i]]
			VB2 <- VB2 + t(MVMX[[i]])%*%rep(1,k[[i]]) #/k[[i]]
			VB3 <- VB3 + sum(MVM[[i]]) #/k[[i]]
			#VB <- Reduce('+',XMVMX)-t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVMX)/sum(Reduce('+',MVM))
			# sum U_i
			#U <- Reduce('+',XMVY) - t(Reduce('+',MVMX))%*%rep(1,k)%*%t(rep(1,k))%*%Reduce('+',MVY)/sum(Reduce('+',MVM))
			U2 <- U2 + t(MVMX[[i]])%*%rep(1,k[[i]]) # /k[[i]]
			U3 <- U3 + t(rep(1,k[[i]]))%*%MVY[[i]]  # /k[[i]]
		}
		VB <- VB1 - VB2%*%t(VB2)/VB3
		U  <- U1 - U2%*%U3/VB3
	#}
	# updated beta
	bb <- beta + solve(VB)%*%U
	
	# conf int
	ci <- cbind(bb,bb) + t(c(-1.96,1.96)%*%t(sqrt(diag(solve(VB)))))
	
	return(list(beta=bb,ci=ci))
}
