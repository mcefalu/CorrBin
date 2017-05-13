rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric')


input <- list()

input[[1]] <- list(n=500,k=5,q=5,b=c(.1,-.1,.05,-.05,0),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[2]] <- list(n=1000,k=2,q=5,b=c(.1,-.1,.05,-.05,0),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[3]] <- list(n=200,k=5,q=5,b=c(.1,-.1,.05,-.05,0),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[4]] <- list(n=500,k=5,q=5,b=c(.1,-.1,.05,-.05,0),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)


i <- 4
for (b in c(log(1),log(1.25),log(2))){
	for (bb in c(log(1),log(1.25),log(2))){
		input[[i+1]] <- list(n=20,k=25,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=T)
		input[[i+2]] <- list(n=20,k=25,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=F)
		input[[i+3]] <- list(n=20,k=25,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=T)
		input[[i+4]] <- list(n=20,k=25,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=F)
		input[[i+5]] <- list(n=50,k=10,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=T)
		input[[i+6]] <- list(n=50,k=10,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=F)
		input[[i+7]] <- list(n=50,k=10,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=T)
		input[[i+8]] <- list(n=50,k=10,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=F)
		i <- i+8
	}
}

i <- 76
for (b in c(log(1),log(1.25),log(2))){
	for (bb in c(log(1),log(1.25),log(2))){
		input[[i+1]] <- list(n=500,k=5,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=T)
		input[[i+2]] <- list(n=500,k=5,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=F)
		input[[i+3]] <- list(n=500,k=5,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=T)
		input[[i+4]] <- list(n=500,k=5,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=F)
		input[[i+5]] <- list(n=1000,k=2,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=T)
		input[[i+6]] <- list(n=1000,k=2,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='obs',cont=F)
		input[[i+7]] <- list(n=1000,k=2,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=T)
		input[[i+8]] <- list(n=1000,k=2,q=2,b=c(b,bb),a=NULL,M=5000,rho=.3,who='other',type='intervention',cont=F)
		i <- i+8
	}
}


### RR for binary treatment
input[[501]] <- list(n=500,k=5,q=5,b=c(log(1.5),log(1.5),log(1.05),log(.95),log(1)),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[502]] <- list(n=500,k=5,q=5,b=c(log(1.5),log(1.5),log(1.05),log(.95),log(1.05)),a=-1.55,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[503]] <- list(n=500,k=5,q=5,b=c(log(1.5),log(1.5),log(1.05),log(.95),log(1.5)),a=-1.9,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[504]] <- list(n=500,k=5,q=5,b=c(log(1.5),log(1.5),log(1.05),log(.95),log(2)),a=-2.2,M=5000,rho=.3,who='me',type=NULL,cont=NULL)

input[[1001]] <- list(n=1000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1)),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[1002]] <- list(n=1000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1.05)),a=-1.55,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[1003]] <- list(n=1000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1.5)),a=-1.9,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[1004]] <- list(n=1000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(2)),a=-2.2,M=5000,rho=.3,who='me',type=NULL,cont=NULL)

input[[5001]] <- list(n=5000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1)),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[5002]] <- list(n=5000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1.05)),a=-1.55,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[5003]] <- list(n=5000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1.5)),a=-1.9,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[5004]] <- list(n=5000,k=2,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(2)),a=-2.2,M=5000,rho=.3,who='me',type=NULL,cont=NULL)

input[[5501]] <- list(n=5000,k=5,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1)),a=-1.5,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[5502]] <- list(n=5000,k=5,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1.05)),a=-1.55,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[5503]] <- list(n=5000,k=5,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(1.5)),a=-1.9,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
input[[5504]] <- list(n=5000,k=5,q=5,b=c(log(1.1),log(.9),log(1.05),log(.95),log(2)),a=-2.2,M=5000,rho=.3,who='me',type=NULL,cont=NULL)


### RR for two continuous trt
i <- 5
for (rr1 in c(1,1.1,1.5,2)){
	for (rr2 in c(1,1.1,1.5,2)){
		input[[500+i]] <- list(n=500,k=5,q=5,b=c(log(rr1),log(rr2),log(1.05),log(.95),log(1.05)),a=-1,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
		input[[1000+i]] <- list(n=1000,k=2,q=5,b=c(log(rr1),log(rr2),log(1.05),log(.95),log(1.05)),a=-1,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
		input[[5000+i]] <- list(n=5000,k=2,q=5,b=c(log(rr1),log(rr2),log(1.05),log(.95),log(1.05)),a=-1,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
		input[[5500+i]] <- list(n=5000,k=5,q=5,b=c(log(rr1),log(rr2),log(1.05),log(.95),log(1.05)),a=-1,M=5000,rho=.3,who='me',type=NULL,cont=NULL)
		i <- i+1
	}
}



save.image('./input/input.RData')