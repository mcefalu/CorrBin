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
		#var.eff <-  score(id=ID.a,Y=Y.a,X=X.a,beta=model.eff$beta,alpha='unspec',time=t.a) 	
		#ci <- cbind(model.eff$beta,model.eff$beta) + t(c(-1.96,1.96)%*%t(sqrt(diag(model.eff$sigma))))
		#model.eff$ci <- ci	

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
		#temp1 <- beta.pois
		beta.pois <- gee(Y.a~X.a,family='poisson',corstr='exchangeable',id=ID.a)$coef[-1]
		model.eff <- score(id=ID.a,Y=Y.a,X=X.a,beta=beta.pois,alpha='unspec',time=t.a)
		#var.eff <-  score(id=ID.a,Y=Y.a,X=X.a,beta=model.eff$beta,alpha='unspec',time=t.a) 
		
		#ci <- cbind(model.eff$beta,model.eff$beta) + t(c(-1.96,1.96)%*%t(sqrt(diag(model.eff$sigma))))
		#model.eff$ci <- ci	

	
		#temp1 <- gee(Y~X,family='poisson',corstr='exchangeable',id=ID)$coef[-1]
		#while (sum(abs(temp1-temp2))>.00001){
		#model.eff <- score(id=ID.a,Y=Y.a,X=X.a,beta=temp1[-1],alpha=temp1[1],time=t.a)
		#	temp2 <- temp1
		#	temp1 <- model.eff$beta
		#}


		return(list(efficient=model.eff))
	}
}
