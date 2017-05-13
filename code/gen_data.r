logit <- function(x){
   log(x/(1-x))
}
expit <- function(x){
   exp(x)/(1+exp(x))
}
library(mvtnorm)
library(binarySimCLF)

gen_data <- function(n , k , q , b , a , type="Exchangeable" , rho=0){
   #n <- 500 # number of individuals
   #k <- 3   # number of obs per cluster
   #q <- 2	 # number of covariates
   
   #b <- c(.3,-.5) # beta for q's
   #a <- 0 # intercepts for individuals
   #R <- correlation matrix
   
   if (type=="AR-1"){
      # AR1
      R = matrix(rho , k ,k)
      for (i in 1:k){
         R[i,i] = 1
      }
      
      H <- abs(outer(1:k, 1:k, "-")) 
      R <-  rho^H 
   }
   if (type=="Exchangeable"){
      # exchangeable correlation matrix.
      R = matrix(rho , k ,k)
      for (i in 1:k){
         R[i,i] = 1
      }
   }
   if (type=="Household"){
      # household correlation matrix.
      R = matrix(0.1 , k ,k)
      R <- R + diag(5)*.9
      
      R[1,2] <- R[2,1] <- .05
      R[3,4] <- R[3,5] <- .3
      R[4,3] <- R[4,5] <- .3
      R[5,3] <- R[5,4] <- .3
   }
   
   dta <- NULL
   for (i in 1:n){
      # k. <- min(max(1,rpois(1,k)),10)
      # R = matrix(rho , k. ,k.)
      # for (ii in 1:k.){
      #    R[ii,ii] = 1
      # }
      #k. <- k
      while (TRUE){
         #temp <- cbind(i,rmvnorm(k,rep(0,q-1)),rbinom(k,1,.5))
         #k. <- max(1,rpois(1,k))
         #print('cluster size:')
         #print(k.)
         # R = matrix(rho , k. ,k.)
         # for (j in 1:k.){
         #    R[j,j] = 1
         # }
         temp <- cbind(i,rep(rbinom(1,1,.5),k),rmvnorm(k,rep(0,4)))
         #colnames(temp)[1] <- 'id'
         #id <- which(temp[,'id']==i)
         mu <- as.vector(exp(temp[,-1]%*%b+a))
         if (all((mu>0.001)&(mu<0.999))){
            if (k>1){
               V = cor2var(R,mu)
               B = allReg(V);
               # Checks CLF compatibility.
               clf.compat = blrchk1(mu,B)
               if (clf.compat){
                  y = mbsclf(1,mu,B)$y
                  dta <- rbind(dta , cbind(temp,t(y)))
                  colnames(dta)[1] <- 'id'
                  break
               }
            }else{
               y <- rbinom(1,1,mu)
               dta <- rbind(dta,cbind(temp,t(y)))
               break
            }
            
         }
      } 
   }
   #dta <- cbind(dta , rbinom(n*k,1,expit(dta[,2:5]%*%c(-1,-1,.3,1.3))))
   #for (i in 1:n){
   #	id <- which(dta[,'id']==i)
   #	mu <- as.vector(exp(dta[id,-1]%*%b+a))
   #	V = cor2var(R,mu)
   #	B = allReg(V);
   #	# Checks CLF compatibility.
   #	clf.compat = blrchk1(mu,B)
   #	if (clf.compat){
   #		y[id] = mbsclf(1,mu,B)$y
   #	} 
   #	else{ print("Not CLF compatible"); return("Not CLF compatible")}
   #}
   #dta <- cbind(dta,as.vector(y))
   colnames(dta)[q+2] <- 'Y'
   colnames(dta)[2:(q+1)] <- paste0("X",1:k)
   return(data.frame(dta))
}
