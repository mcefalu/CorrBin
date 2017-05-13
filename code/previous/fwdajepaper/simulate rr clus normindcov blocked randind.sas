**********************************************************************************************************************;
*																						 	 					  	 *;
* Macro to Simulate Clustered Data with Blocked Randomisation of Individuals and a Normal Individual Level Covariate *;
*																						 	 					  	 *;
**********************************************************************************************************************;



/* data = name of dataset
   nsim = number of simulated datasets contained in data
   nclus = number of clusters
   fixed = whether cluster size is fixed or varies according to a poisson distribution
   nobsclus = (mean) number of observations per cluster
   p0 = expected probability of success in control group (not the same as exp(beta_0) since addition of random effects on the 
		log scale changes the mean on the original scale from exp(beta_0) to exp(beta_0+varre/2) => beta_0=log(p0/exp(varre/2))
   rrt = relative risk (treatment)
   varre = variance of random effects
   meannc = mean of covariate
   varnc = variance of covariate
   rhonc = intracluster correlation of covariate
   rrnc = relative risk (covariate)

   the macro:
   1) determines the size of each cluster (fixed or randomly generated from a poisson distribution)
   2) assigns half of all patients to intervention group and half to control group (this would occur at least 
      approximately in practice if randomly permuted blocks were used)
   3) generates a random effect for each cluster
   4) generates a covariate from the model x_ij = meannc + a_i + b_ij where a_i~N(0,sigma1) and b_ij~N(0,sigma2)
   5) resamples from b_ij distribution until probability < 1 if necessary (will reduce sigma2 and hence increase icc)
   6) randomly allocates patients to have the outcome with probability (p0/exp(varre/2))*exp(log(rrt)*treatment)*exp(log(rrnc)*normcov)*exp(random)*

   note: ignore pbc, rhobc, rrbc (these are used to simulate more complicated datasets) */

%macro simulate(data, nsim, nclus, fixed, nobsclus, p0, rrt, varre, pbc, rhobc, rrbc, meannc, varnc, rhonc, rrnc);

	data &data.;
		length simulation cluster_id cluster_size subject_id intercept treatment random normcov outcome 3.;
	  	do simulation = 1 to &nsim.;
	   		do cluster_id = 1 to &nclus.;
				if &fixed. = 1 then cluster_size = &nobsclus.;
				else cluster_size = ranpoi(cluster_id*simulation, &nobsclus.);
				intercept = 1;
				random = sqrt(&varre.)*rannor(2*cluster_id*simulation);
				%if &rrnc. = 1 %then %do; 	* added 8/6/10 to avoid invalid values of random effect when covariate has no effect;
					prob = (&p0./exp(&varre./2))*&rrt.*exp(random);
					j = 3;
					do while (prob >= 1);
						j = j + 1;
						random = sqrt(&varre.)*rannor(j*cluster_id*simulation);
						prob = (&p0./exp(&varre./2))*&rrt.*exp(random);
					end;
					drop j;
				%end;
				ai = sqrt(&rhonc.*&varnc.)*rannor(3*cluster_id*simulation);
				do subject_id = 1 to cluster_size;
					treatment = .;
					normcov = .;	
					outcome = .;
					uniform = uniform(cluster_id*subject_id*simulation);		* use this to allocate treatment;
					* select value of bij that makes probability ok if patient assigned to control group;
					i = 1;
					do until (prob1 < 1);
	    				i = i + 1; 	* seed starts at 2*cluster_id*subject_id*simulation;
						bij1 = sqrt(&varnc.*(1-&rhonc.))*rannor(i*cluster_id*subject_id*simulation);
						normcov1 = &meannc. + ai + bij1;	
						prob1 = (&p0./exp(&varre./2))*exp(log(&rrnc.)*normcov1)*exp(random);
					end;
					* select value of bij that makes probability ok if patient assigned to invervention group;
					normcov2 = normcov1;
					prob2 = prob1*&rrt.;
					do while (prob2 >= 1);
	    				i = i + 1; 	
						bij2 = sqrt(&varnc.*(1-&rhonc.))*rannor(i*cluster_id*subject_id*simulation);
						normcov2 = &meannc. + ai + bij2;	
						prob2 = (&p0./exp(&varre./2))*&rrt.*exp(log(&rrnc.)*normcov2)*exp(random);
					end;
	    			output;
				end;	
	 		end;
		end;    
	run;

	proc means data=&data. noprint;
		by simulation;
		var intercept;
		output out=temp sum(intercept)=nobs; 				* calculate total number of patients;
	run;

	data temp;
		set temp;
		rand = ranbin(simulation+1000000, 1, 0.5);			* treatment allocation of first half of patients (excluding extra patient if odd number);
	run;

	data &data.;
		merge &data. temp(drop=_type_ _freq_);
		by simulation;
	run;

	proc sort data=&data.;
		by simulation uniform;
	data &data.;
		set &data.;
		by simulation;
		retain n;
		if first.simulation then n = 1;
		else n = n + 1;
		if n <= int(nobs/2) then treatment = rand;
		else treatment = 1 - rand;
		if treatment = 0 then normcov = normcov1;
		else normcov = normcov2;
		outcome = ranbin(cluster_id*subject_id*simulation, 1, (&p0./exp(&varre./2))*exp(log(&rrt.)*treatment)*exp(log(&rrnc.)*normcov)*exp(random));
		drop i prob1 prob2 normcov1 normcov2 uniform nobs rand n ai bij1 bij2;
	run;

	proc sort data=&data.;
		by simulation cluster_id subject_id;
	run;

%mend;


/* test macro;

%simulate(test, 1000, 50, 0, 10, 0.1, 2, 0.1, , , , 0.5, 0.25, 0.2, 2);

* check correct proportion of successes by treatment group;

proc freq data=test;
	tables treatment*outcome / nocol nopercent;
run;

* check correct proportion of successes by covariate category;

data test;
	set test;
	if normcov < 0.5 - sqrt(0.25) then catcov = 1;
	if 0.5 - sqrt(0.25) <= normcov < 0.5 then catcov = 2;
	if 0.5 <= normcov < 0.5 + sqrt(0.25) then catcov = 3;
	if normcov >= 0.5 + sqrt(0.25) then catcov = 4;
run;

proc freq data=test;
	tables catcov*treatment*outcome / nocol nopercent;
run;

* check correct number of treatment allocations by simulation;

proc freq data=test;
	tables simulation*treatment / norow nocol nopercent;
run;

* check distribution of covariate, random effects and cluster size;

proc univariate data=test;
	var normcov; 		
	histogram normcov;
run;

data temp;
	set test;
	by simulation cluster_id;
	if first.cluster_id;
run;

proc univariate data=temp;
	var random cluster_size; 		
	histogram random cluster_size;
run;

* check distribution of proportion of successes in each cluster (variance shows overdispersion due to clustering as expected);

proc means data=test noprint;
	by simulation cluster_id;
	var outcome;
	output out=temp mean(outcome)=;
run;

proc univariate data=temp;
	var outcome;
run;

* check for zero cluster sizes (temp will have nsim*nclus observations if there are no clusters of size zero);

data temp;
	set test;
	by simulation cluster_id;
	if first.cluster_id;
run; */





