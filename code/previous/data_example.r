## data example

## load score function

dta <- read.csv("C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric\\data.csv")

Y <- dta$Y
ID <- dta$ID
time <- dta$time
X <- dta$dose
m1 <- my.estimator(Y=Y,X=X,ID=ID,time=time)

## same model as paper
X <- cbind(dta$dose*time,dta$dose*time^2)
m2 <- my.estimator(Y=Y,X=X,ID=ID,time=time)
