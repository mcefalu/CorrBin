library(survival)
my.estimator <- function(Y,X,ID,time=NULL,beta.pois=NULL){
	if (is.null(beta.pois)){
		X <- as.matrix(X)
		# initial fit 
		if (is.null(time)){
			cox.model <- coxph(Surv(rep(1,length(Y)),Y)~X)
		}
		else{
			cox.model <- coxph(Surv(rep(1,length(Y)),Y)~X+strata(time))
		}
		# excludes any NA is outcome
		Y.a  <- Y[!is.na(Y)]
		X.a  <- as.matrix(X[!is.na(Y),])
		ID.a <- ID[!is.na(Y)]
		t.a <- time[!is.na(Y)]
		# assumes independence
		model.uni <- score(id=1:length(Y.a),Y=Y.a,X=as.matrix(X.a),beta=cox.model$coef,alpha='unspec',time=t.a)
		# exchangeable correlation structure
		model.eff <- score(id=ID.a,Y=Y.a,X=X.a,beta=model.uni$beta,alpha='unspec',time=t.a)
		# more updates
		#model2 <- score(id=ID.a,Y=Y.a,X=X.a,beta=model.eff$beta,alpha='unspec',time=t.a)
		#model3 <- score(id=ID.a,Y=Y.a,X=X.a,beta=model2$beta,alpha='unspec',time=t.a)
		#model4 <- score(id=ID.a,Y=Y.a,X=X.a,beta=model3$beta,alpha='unspec',time=t.a)
		return(list(efficient=model.eff , independent=model.uni ))
	}
	else{
		X <- as.matrix(X)
		# excludes any NA is outcome
		Y.a  <- Y[!is.na(Y)]
		X.a  <- as.matrix(X[!is.na(Y),])
		ID.a <- ID[!is.na(Y)]
		t.a <- time[!is.na(Y)]
		# assumes independence
		model <- score(id=ID.a,Y=Y.a,X=X.a,beta=beta.pois,alpha='unspec',time=t.a)
		return(list(efficient=model))
	}
}
