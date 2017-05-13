#####
rm(list=ls())
data <- read.csv('./data/YC/dataDec10.csv',header=T,as.is=T)

# y=1, completely true; y=2, somewhat true / dont know ; y=3,4 worse
#y <- as.matrix(data[,'Q70_hivworse'])
#y <- (y==1)*1 + (y==2)*1

# new outcome
y = data[,'Kids_HIV_Teach_STD']
y <- (y>3)*1


# covariates
other <- (data$ethnic%in%c(5,3,9))*1
chagga <- (data$ethnic==1)*1
pare <- (data$ethnic==2)*1 

#weatlh -- 0 is baseline
wealth1 <- (data$wealth==1)*1
wealth2 <- (data$wealth==2)*1
wealth3 <- (data$wealth==3)*1
wealth4 <- (data$wealth==4)*1
wealth5 <- (data$wealth==5)*1

#religion
religion1 <- ((data$religion==1)|(data$religion==4))*1
religion2 <- (data$religion==2)*1
religion3 <- (data$religion==3)*1
religion4 <- (data$religion==5)*1


#x = cbind(as.matrix(data[,c('group','mtaa_mean','knowmtaa','goodfloor')]),pare,other)
x = cbind(as.matrix(data[,c('group','urb_rur')]),pare,other,religion1,religion2,religion3,wealth1,wealth2,wealth3,wealth4,wealth5)
#x = cbind(as.matrix(data[,c('group','urb_rur')]),pare,other,religion1,religion2,religion3)

# 
id = data[,'com_id']

model <- binary.gee(Y=y,X=x,ID=id)
model <- binary.gee(Y=y,X=x,ID=id,beta.init=model$beta)


coef(gee(y~x,family='poisson',corstr='exchangeable',id=id))[2]+c(-1.96,1.96)*sqrt(gee(y~x,family='poisson',corstr='exchangeable',id=id)$robust[2,2])

sqrt(gee(y~x,family='binomial',corstr='exchangeable',id=id)$robust)


tab <- cbind(c(model$beta[1],))



# bootstrap 
M=50
beta.eff = numeric(M)
beta.pois = numeric(M)
beta.ind = numeric(M)

for (i in 1:M){
   index = sample(unique(id),length(unique(id)),replace=T)
   x.boot = x[unlist(sapply(index,function(x) which(id==x))),]
   y.boot = unlist(sapply(index,function(x) y[id==x]))
   id.boot=NULL
   for (k in 1:length(index)){
      id.boot = c(id.boot, rep(k,sum(id==index[k])))
   }
   beta.pois[i]=coef(gee(y.boot~x.boot,family='poisson',corstr='exchangeable',id=id.boot))[2]
   model <- binary.gee(Y=y.boot,X=x.boot,ID=id.boot)
   model <- binary.gee(Y=y.boot,X=x.boot,ID=id.boot,beta.init=model$beta)
   beta.eff[i] = model$beta[1]
}

sd(beta.eff)
sd(beta.pois)

               
               
               


#######################
#######################

# kids a big problem versus not
y1 <- (data$Q51==1)*1
# robbery often versus not
y2 <- (data$Q65==1)*1
# violent arguement often versus not
y3 <- (data$Q67==1)*1
# rape yes vs no
y4 <- (data$Q68==1)|(data$Q68==2)*1
# AIDS better or worse
y3 <- (data$Q70==1)*1+(data$Q70==2)*1


#######################
#######################
# new analysis 

data <- read.csv('./data/YC/dataDec10.csv',header=T,as.is=T)

# y=1, completely true; y=2, somewhat true / dont know ; y=3,4 worse
y <- as.matrix(data[,'Q70_hivworse'])
y <- (y==1)*1 + (y==2)*1






