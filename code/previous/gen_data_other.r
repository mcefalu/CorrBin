logit <- function(x){
	log(x/(1-x))
}
expit <- function(x){
	exp(x)/(1+exp(x))
}
library(mvtnorm)
library(binarySimCLF)

# n is number of clusters
# k is number per cluster

gen.data.other <- function(n=20 , k=25 , b.trt=log(1) , b.cov=log(1) , b0=.1 , rho.cov=0.1 , R , type='obs' , continuous=T){
	# generate id
	dta <- NULL
	for (i in 1:n){
		dta <- c(dta , rep(i,k))
	}
	dta <- as.matrix(dta)
	colnames(dta)[1] <- 'id'
	
	R.c = matrix(rho.cov , k ,k)
	#r2 = c(0.1, 1, 0.3)
	#r3 = c(0.2, 0.3, 1)
	for (i in 1:k){
		R.c[i,i] = 1
	}
	
	if (type=='intervention'){
		# generate binary covariate
		X <- numeric(n*k)
		x1 <- numeric(n*k)
		a <- rnorm(n , 0 , sqrt(.25/2)) # cluster specific effect for covariate
		y <- numeric(n*k)
		u <- rnorm(n,0,sqrt(.1)) # cluster specific effect for 
		for (i in 1:n){
			while (T){
				# cluster size
				# k  <- rpois(1 , mu.k)
				a[i] <- rnorm(1 , 0 , sqrt(.05))
				id <- which(dta[,'id']==i)
				X[id] <- (i<=(n/2))*1
				if (!continuous){
					mu <- rep(.5,k)
					V = cor2var(R.c,mu)
					B = allReg(V)
					x1[id] = mbsclf(1,mu,B)$y
				}
				else{
					# continuous cov
					x1[id] <- .5 + a[i] + rnorm(k,0,sqrt(.25/2))
				}
				# outcome
				mu <- as.vector(exp(log(b0) + b.trt*X[id] + b.cov*x1[id] + u[i]))
				if (all(mu<1)){
					V = cor2var(R,mu)
					B = allReg(V)
					# Checks CLF compatibility.
					clf.compat = blrchk1(mu,B)
					if (clf.compat){
						y[id] = mbsclf(1,mu,B)$y
						break
					}
				}
			}
		}
	}
	else{
		# generate binary covariate
		X <- numeric(n*k)
		x1 <- numeric(n*k)
		a <- rnorm(n , 0 , sqrt(.25/2)) # cluster specific effect for covariate
		y <- numeric(n*k)
		u <- rnorm(n,0,sqrt(.1)) # cluster specific effect for 
		for (i in 1:n){
			while (T){
				# cluster size
				# k  <- rpois(1 , mu.k)
				a[i] <- rnorm(1 , 0 , sqrt(.05))
				id <- which(dta[,'id']==i)
				X[id] <- (i<=(n/2))*1
				if (!continuous){
					if (i<=(n/2)){
						mu <- rep(.6,k)
					}else{ mu <- rep(.4,k) }
					V = cor2var(R.c,mu)
					B = allReg(V)
					x1[id] = mbsclf(1,mu,B)$y
				}else{
					# continuous cov
					if (i<=(n/2)){
						x1[id] <- .6 + a[i] + rnorm(k,0,sqrt(.25/2))
					} else{ x1[id] <- .4 + a[i] + rnorm(k,0,sqrt(.25/2)) }
				}
				# outcome
				mu <- as.vector(exp(log(b0) + b.trt*X[id] + b.cov*x1[id] + u[i]))
				if (all(mu<1)){
					V = cor2var(R,mu)
					B = allReg(V)
					# Checks CLF compatibility.
					clf.compat = blrchk1(mu,B)
					if (clf.compat){
						y[id] = mbsclf(1,mu,B)$y
						break
					}
				}
			}
		}
	}
	
	dta <- cbind(dta,as.vector(X),as.vector(x1),as.vector(y))
	colnames(dta)[4] <- 'Y'
	colnames(dta)[3] <- 'C'
	colnames(dta)[2] <- 'X'
	return(dta)
}