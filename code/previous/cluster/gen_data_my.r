logit <- function(x){
	log(x/(1-x))
}
expit <- function(x){
	exp(x)/(1+exp(x))
}
library(mvtnorm)
library(binarySimCLF)

gen.data.my <- function(n , k , q , b , a , R){
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
			#colnames(temp)[1] <- 'id'
			#id <- which(temp[,'id']==i)
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
	return(dta)
}
