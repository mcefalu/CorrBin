
rm(list=ls())

source('./code/binary_gee.r')
source('./code/gen_data.r')

###########################
###########################

num <- 1
rho = .3
type= "Household"

set.seed(num+341)

n <- 500
k <- 5  #5 -- number per cluster
q <- 5 #5 -- number of covariates
b <- log(c(1,1.1,1.2,1.3,1.4)) #c(.1,-.1,.05,-.05,0) -- betas 
a <- -0.5  #-1.5 -- intercept
M <- 100 #10000 -- replications
b[1] = log(seq(1,1.1,by=0.0025))[num]

results <- data.table()

start.time <- Sys.time()
for (i in 1:M){
   print(i)
   dta <- 'temp'
   count <- 0
   while (class(dta)=="character"){
      dta <- gen_data(n=n , k=k , q=q , b=b , a=a , rho=rho , type=type)
      count <- count+1
      if (count>1000){
         print(dta)
         break
      }
   }
   
   # fit model with independent structure
   model.ind <- binary.gee(Y=dta$Y,X=dta[,grep("X",colnames(dta),value=T)],ID=1:nrow(dta))
   
   # fit modified poisson w/exchangeable or AR-1
   pois.exch <- gee(as.formula(paste0("Y~",paste0(grep("X",colnames(dta),value=T),collapse="+"))),family='poisson',corstr='exchangeable',id=id,data=dta)
   pois.ar1 <- gee(as.formula(paste0("Y~",paste0(grep("X",colnames(dta),value=T),collapse="+"))),family='poisson',corstr='AR-M',id=id,data=dta,Mv=1)
   
   # fit efficient estimator with exchangeable or ar1
   eff.exch <- binary.gee(Y=dta$Y,X=dta[,grep("X",colnames(dta),value=T)],ID=dta$id,beta.init=coef(pois.exch)[-1])
   eff.ar1 <- binary.gee(Y=dta$Y,X=dta[,grep("X",colnames(dta),value=T)],ID=dta$id,beta.init=coef(pois.ar1)[-1],corr='AR-1')
   
   # save results
   ind = data.table(n=n,rr=exp(b),covariate=c("binary",rep("continuous",k-1)),true.corr=type,estimator="Efficient",corr="Independent",coef=model.ind$beta[,1],var=diag(model.ind$sigma.hat),lb=model.ind$ci[,1],ub=model.ind$ci[,2])
   eff_exch = data.table(n=n,rr=exp(b),covariate=c("binary",rep("continuous",k-1)),true.corr=type,estimator="Efficient",corr="Exchangeable",coef=eff.exch$beta[,1],var=diag(eff.exch$sigma.hat),lb=eff.exch$ci[,1],ub=eff.exch$ci[,2])
   eff_ar1 = data.table(n=n,rr=exp(b),covariate=c("binary",rep("continuous",k-1)),true.corr=type,estimator="Efficient",corr="AR-1",coef=eff.ar1$beta[,1],var=diag(eff.ar1$sigma.hat),lb=eff.ar1$ci[,1],ub=eff.ar1$ci[,2])
   
   pois_exch = data.table(n=n,rr=exp(b),covariate=c("binary",rep("continuous",k-1)),true.corr=type,estimator="Poisson",corr="Exchangeable",coef=coef(pois.exch)[-1],var=diag(pois.exch$robust.variance)[-1])
   pois_exch[ , lb := coef - 1.96*sqrt(var) ]
   pois_exch[ , ub := coef + 1.96*sqrt(var) ]
   
   pois_ar1 = data.table(n=n,rr=exp(b),covariate=c("binary",rep("continuous",k-1)),true.corr=type,estimator="Poisson",corr="AR-1",coef=coef(pois.ar1)[-1],var=diag(pois.ar1$robust.variance)[-1])
   pois_ar1[ , lb := coef - 1.96*sqrt(var) ]
   pois_ar1[ , ub := coef + 1.96*sqrt(var) ]
   
   # store
   results <- rbind( results , ind , eff_exch , eff_ar1 , pois_exch , pois_ar1)
}
Sys.time()-start.time

save.image(paste('./results/',type,'/output',num,'.RData',sep=''))



results[covariate=="binary",var(coef),by=.(n,covariate,rr,true.corr,estimator,corr)]
