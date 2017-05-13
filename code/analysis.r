#################GBC######################
rm(list=ls())
setwd('C:\\Users\\Matt\\Documents\\School\\Harvard\\Research\\eric')
source('./code/score.r')
source('./code/my_estimator.r')

data <- read.csv('./data/GBM-master.csv')
Y <- data$obs30
X <- cbind( data[,c("sex",'age','white','Black','Other')] , data[,46:79])

# treatment variables
Z <- data[,c("Excision","Biopsy_and_Excision","Biopsy_only")]
Z[,'Excision'] <- Z[,'Excision']*(1-Z[,'Biopsy_and_Excision'])
names(Z)[1] <- 'Excision_only'

# try some things
X <- as.matrix(X)
Z <- as.matrix(Z)

# rm cancers and colinear things
X <- X[,c(-5,-6,-19,-27,-32)]

subs <- 1:1000

glm(Y~X+Z,family=binomial(link='logit'))
glm(Y~X+Z,family=binomial(link='log'))
glm(Y~X+Z,family=poisson(link='log'))

# fit my estimator
aa <- unique(data$prov_num)
ID <- numeric(length(Y))
for (i in 1:length(Y)){
	ID[i] = which(data$prov_num[i]==aa)
}

a <- my.estimator(Y=Y[subs],X=cbind(X,Z)[subs,],ID=ID[subs])