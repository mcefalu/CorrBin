## bayesian correlated binary ##

rm(list=ls())

loglik <- function(alpha,beta,rho,...){
   mu <- apply(X,1,function(x) exp(x%*%beta+alpha))
   if (all(mu<1)){
      R <- (Y-mu)/sqrt(mu*(1-mu))
      l.star <- Y%*%log(mu)+(1-Y)%*%log(1-mu)
      b.star <- 0
      for (id in ID){
         temp <- R[ID==id]
         temp <- sapply(temp,function(x) rho*x*temp)
         for (i in 2:(dim(temp)[1])){
            #temp1 <- temp[1:i,1:i]
            if (sum(temp[i,1:(i-1)])<(-1)){
               return(-Inf)
            }else{
               b.star <- b.star + log(1+sum(temp[i,1:(i-1)]))
            }
         }
      }
      return(l.star+b.star)
   }
   else{
      return(-Inf)
   }
}

logit <- function(x){
   log(x/(1-x))
}
expit <- function(x){
   exp(x)/(1+exp(x))
}

gen.data.my <- function(n , k , q , b , a , R){
   library(mvtnorm)
   library(binarySimCLF)
   #n <- 500 # number of individuals
   #k <- 3   # number of obs per cluster
   #q <- 2	 # number of covariates
   
   #b <- c(.3,-.5) # beta for q's
   #a <- 0 # intercepts for individuals
   #R <- correlation matrix
   dta <- NULL
   y <- numeric(n*k)
   for (i in 1:n){
      while (TRUE){
         temp <- cbind(i,rmvnorm(k,rep(0,q-1)),rbinom(k,1,.5))
         #temp <- cbind(i,rmvnorm(k,rep(0,q)))
         mu <- as.vector(exp(temp[,-1]%*%b+a))
         if (all((mu>0.05)&(mu<0.95))){
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
         }
      } 
   }
   colnames(dta)[q+2] <- 'Y'
   return(dta)
}

R = matrix(.3 , 5 ,5)
for (i in 1:5){
   R[i,i] = 1
}

rr1 = 1.05
rr2 = 1.05
dta <- gen.data.my(100,5,5,c(log(rr1),log(rr2),log(1.05),log(.95),log(1.05)),log(exp(-1)),R)


Y <- dta[,'Y']
X <- dta[,2:6]
ID <- dta[,'id']



model <- binary.gee(Y=Y,X=X,ID=ID,beta.init=c(log(rr1),log(rr2),log(1.05),log(.95),log(1.05)))
beta. = model$beta
beta.var = model$sigma.hat
beta.var[1,2]=beta.var[2,1]
beta.var = diag(diag(beta.var))

i=1
beta = matrix(0,1000,5)
rho = numeric(1000)
alpha = numeric(1000)
while (i <= 1000){
   beta.prop = t(rmvnorm(1,beta.,beta.var))
   rho.prop = runif(1,.2,.4)
   alpha.prop = rnorm(1,log(.4))
   num = loglik(Y=Y,X=X,beta=beta.prop,rho=rho.prop,alpha=alpha.prop,ID=ID)
   den = loglik(Y=Y,X=X,beta=beta.prop,rho=rho.prop,alpha=alpha.prop,ID=ID)
   if( (num!=-Inf)&(den!=-Inf)){
   if (exp(num-den)>runif(1)){
      beta[i,] = beta.prop
      alpha[i] = alpha.prop
      rho[i] = rho.prop
      i = i+1
   }
   }
}

apply(beta,2,var)/diag(beta.var)


loglik(alpha=log(.4),beta=c(log(1.1),log(1.05),log(1.5),log(1.5),log(.9)),Y=Y,X=X,ID=ID,rho=.3)
