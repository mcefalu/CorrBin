#####################
# binary_gee.r
# 03/22/15
# written by m. cefalu
# no warranty or copyright -- free to modify or distribute

library(gee)

#####################
# the main function is binary.gee()
# it does a one step update of either a user supplied beta.init or a Poisson GEE model fit

#####################
# Y is the outcome -- must be binary
# X are the covariates
# ID is the cluster ID
# beta.init is an optional initial coefficient vector
# corr is the correlation structure -- exchangeable or AR-1 -- exchangeable is the default
# data needs to be in long format

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

# ###########################
# ####### Example ###########
# 
# # this function requires the library gee to find initial model fit
# # so install it if it is not already available
# # library(gee)
# # install.packages('gee')
# 
# # Generate some clustered binary data in long format
# x <- matrix(rnorm(200),100,2) # the covariates
# y <- rbinom(100,1,.5) # the outcome -- the model is Y~X
# id <- rep(1:50,2) # id for clusters -- needs to be indexed from 1 to # of clusters
# 
# model <- binary.gee(Y=y,X=x,ID=id)
# 
# #############################
# # what is saved in the model?
# 
# # the estimated RR
# model$beta
# # the std errors
# sqrt(diag(model$sigma.hat))
# # the conf int
# model$ci
# # the test stats
# model$test.stat
# # the p-values
# model$p.value
# 
# 
# #######################################################
# ####### Example -- independent observations ###########
# 
# # with independent observations, there are no clusters
# # so set ID=1:n , where n is the # of observations
# n <- length(y)
# model.ind <- binary.gee(Y=y,X=x,ID=1:n)
# 
# # if data was actually clustered -- we could use this as an initial fit
# model <- binary.gee(Y=y,X=x,ID=id,beta.init=model.ind$beta)




