rm(list=ls())
dta <- read.csv('./data/arrhythmia.csv',header=T,as.is=T)

ID = dta[,'ID']
Y = apply(dta[,c('ALLSVE','ALLVE')]>0,1,any)*1
#Y=(dta[,c('ALLSVE')]>0)*1
X = dta[,c('c_pm_d5','age','bmi','c_so_d5')]#,'apptemp_l0')]
#X = dta[,c('c_pm_d5','age','bmi')]#,'c_so_d5')]#,'apptemp_l0')]

#X = dta[,c(paste('c_pm_d',1:5,sep=''),'age','bmi',paste('c_so_d',1:5,sep=''))]#,'apptemp_l0')]

index = !apply(t(apply(X,1,function(x) is.na(x))),1,any)
#index = index & sapply(ID,function(x){!(x%in%names(table(ID[index])[which(table(ID[index])==1)])})
index = index & (ID!=30) & (ID!=31)


model <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index])
model <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index],beta.init=model$beta)

model2 <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index],corr='AR-1')
model2 <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index],beta.init=model2$beta,corr='AR-1')


md = gee(Y[index]~as.matrix(X[index,]),family='binomial',corstr='AR-M',id=ID[index])

lb = 5*coef(md)-1.96*5*sqrt(diag(md$robust))
ub = 5*coef(md)+1.96*5*sqrt(diag(md$robust))

se = sqrt(diag(md$robust))


cbind(coef(md),se,2*pnorm(abs(coef(md)/se),lower.tail=FALSE),lb,ub)[c(2,5),]
se*5model

mdo = gee(Y[index]~as.matrix(X[index,]),family='binomial',corstr='exchangeable',id=ID[index])
lbo = coef(mdo)-1.96*sqrt(diag(mdo$robust))
ubo = coef(mdo)+1.96*sqrt(diag(mdo$robust))
seo = sqrt(diag(mdo$robust))


cbind(coef(mdo),seo,2*pnorm(abs(coef(mdo)/seo),lower.tail=FALSE),lbo,ubo)[c(2,5),]

# independent
model <- binary.gee(Y=Y[index],X=X[index,],ID=1:length(ID[index]))
model <- binary.gee(Y=Y[index],X=X[index,],ID=1:length(ID[index]),beta.init=model$beta)





## SVE pm only
Coefficients:
   Estimate  StdError TestStat      Pvalue
c_pm_d5 0.0019162 0.0017573    1.090 0.275526122
age     0.0230236 0.0070579    3.262 0.001105965
bmi     0.0114430 0.0166844    0.686 0.492808232












## look at all lags
b.sulfur = matrix(0,2,5)
b.pm = matrix(0,2,5)
for (i in 1:5){
   X = dta[,c(paste('c_pm_d',i,sep=''),'age','bmi')]#,'apptemp_l0')]
   index = !apply(t(apply(X,1,function(x) is.na(x))),1,any)
   model <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index])
   for (j in 1:20){
      model <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index],beta.init=model$beta)
   }
   b.pm[1,i]=model$beta[1]
   md = gee(Y[index]~as.matrix(X[index,]),family='poisson',corstr='exchangeable',id=ID[index])
   b.pm[2,i] = coef(md)[2]
   lb = coef(md)-1.96*sqrt(diag(md$robust))
   ub = coef(md)+1.96*sqrt(diag(md$robust))
   se = sqrt(diag(md$robust))   
   
   X = dta[,c(paste('c_so_d',i,sep=''),'age','bmi')]#,'apptemp_l0')]
   index = !apply(t(apply(X,1,function(x) is.na(x))),1,any)
   model <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index])
   for (j in 1:20){
      model <- binary.gee(Y=Y[index],X=X[index,],ID=ID[index],beta.init=model$beta)
   }
   b.sulfur[1,i]=model$beta[1]
   md = gee(Y[index]~as.matrix(X[index,]),family='poisson',corstr='exchangeable',id=ID[index])
   b.sulfur[2,i] = coef(md)[2]
   lb = coef(md)-1.96*sqrt(diag(md$robust))
   ub = coef(md)+1.96*sqrt(diag(md$robust))
   se = sqrt(diag(md$robust))   
}







