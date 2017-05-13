********************************************************************************************;
*																						   *;
* Macro to Simulate Clustered Data with Confounding by a Normal Individual Level Covariate *;
*																						   *;
********************************************************************************************;



/* data = name of dataset
   nsim = number of simulated datasets contained in data
   nclus = number of clusters
   fixed = whether cluster size is fixed or varies according to a poisson distribution
   nobsclus = (mean) number of observations per cluster
   p0 = expected probability of success in unexposed group (not the same as exp(beta_0) since addition of random effects on the 
		log scale changes the mean on the original scale from exp(beta_0) to exp(beta_0+varre/2) => beta_0=log(p0/exp(varre/2))
   rre = relative risk (exposure)
   varre = variance of random effects
   meanncu = mean of covariate (unexposed)
   meannce = mean of covariate (exposed)
   varnc = variance of covariate
   rhonc = intracluster correlation of covariate
   rrnc = relative risk (covariate)

   the macro:
   1) determines the size of each cluster (fixed or randomly generated from a poisson distribution)
   2) assigns half of all patients to exposed group and half to unexposed group
   3) generates a random effect for each cluster
   4) generates a covariate from the model x_ij = meannc_ij + a_i + b_ij where a_i~N(0,sigma1), b_ij~N(0,sigma2) and meannc_ij = meanncu or meannce
   5) resamples from b_ij distribution until probability < 1 if necessary (will reduce sigma2 and hence increase icc)
   6) randomly allocates patients to have the outcome with probability (p0/exp(varre/2))*exp(log(rre)*exposure)*exp(log(rrnc)*normcov)*exp(random)*

   note: ignore pbcu, pbce, rhobc, rrbc (these are used to simulate more complicated datasets) */

%macro simulate(data, nsim, nclus, fixed, nobsclus, p0, rre, varre, pbcu, pbce, rhobc, rrbc, meanncu, meannce, varnc, rhonc, rrnc);

	data &data.;
		length simulation cluster_id cluster_size subject_id intercept exposed random normcov outcome 3.;
	  	do simulation = 1 to &nsim.;
	   		do cluster_id = 1 to &nclus.;
				if &fixed. = 1 then cluster_size = &nobsclus.;
				else cluster_size = ranpoi(cluster_id*simulation, &nobsclus.);
				intercept = 1;
				random = sqrt(&varre.)*rannor(2*cluster_id*simulation);
				%if &rrnc. = 1 %then %do; 	* added 9/6/10 to avoid invalid values of random effect when covariate has no effect;
					prob = (&p0./exp(&varre./2))*&rre.*exp(random);
					j = 3;
					do while (prob >= 1);
						j = j + 1;
						random = sqrt(&varre.)*rannor(j*cluster_id*simulation);
						prob = (&p0./exp(&varre./2))*&rre.*exp(random);
					end;
					drop j;
				%end;
				ai = sqrt(&rhonc.*&varnc.)*rannor(3*cluster_id*simulation);
				do subject_id = 1 to cluster_size;
					exposed = .;
					normcov = .;	
					outcome = .;
					uniform = uniform(cluster_id*subject_id*simulation);		* use this to allocate exposure;
					* select value of bij that makes probability ok if patient unexposed;
					i = 1;
					do until (prob1 < 1);
	    				i = i + 1; 	* seed starts at 2*cluster_id*subject_id*simulation;
						bij1 = sqrt(&varnc.*(1-&rhonc.))*rannor(i*cluster_id*subject_id*simulation);
						normcov1 = &meanncu. + ai + bij1;	
						prob1 = (&p0./exp(&varre./2))*exp(log(&rrnc.)*normcov1)*exp(random);
					end;
					* select value of bij that makes probability ok if patient exposed;
					normcov2 = &meannce. + ai + bij1;
					prob2 = (&p0./exp(&varre./2))*&rre.*exp(log(&rrnc.)*normcov2)*exp(random);
					do while (prob2 >= 1);
	    				i = i + 1; 	
						bij2 = sqrt(&varnc.*(1-&rhonc.))*rannor(i*cluster_id*subject_id*simulation);
						normcov2 = &meannce. + ai + bij2;	
						prob2 = (&p0./exp(&varre./2))*&rre.*exp(log(&rrnc.)*normcov2)*exp(random);
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
		rand = ranbin(simulation+1000000, 1, 0.5);			* exposure allocation of first half of patients (excluding extra patient if odd number);
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
		if n <= int(nobs/2) then exposed = rand;
		else exposed = 1 - rand;
		if exposed = 0 then normcov = normcov1;
		else normcov = normcov2;
		prob = (&p0./exp(&varre./2))*exp(log(&rre.)*exposed)*exp(log(&rrnc.)*normcov)*exp(random);
		outcome = ranbin(cluster_id*subject_id*simulation, 1, (&p0./exp(&varre./2))*exp(log(&rre.)*exposed)*exp(log(&rrnc.)*normcov)*exp(random));
		drop i prob1 prob2 normcov1 normcov2 uniform nobs rand n ai bij1 bij2;
	run;

	proc sort data=&data.;
		by simulation cluster_id subject_id;
	run;

%mend;


*simulate(data, nsim, nclus, fixed, nobsclus, p0, rre, varre, pbcu, pbce, rhobc, rrbc, meanncu, meannce, varnc, rhonc, rrnc);
%simulate(test, 1000, 50, 1,       10,      0.1, 2,   0.1,       ,     ,     ,      , 0.5     , 0.5    , 0.25, 0.05, 2);

/* test macro;
data = name of dataset
   nsim = number of simulated datasets contained in data
   nclus = number of clusters
   fixed = whether cluster size is fixed or varies according to a poisson distribution
   nobsclus = (mean) number of observations per cluster
   p0 = expected probability of success in unexposed group (not the same as exp(beta_0) since addition of random effects on the 
		log scale changes the mean on the original scale from exp(beta_0) to exp(beta_0+varre/2) => beta_0=log(p0/exp(varre/2))
   rre = relative risk (exposure)
   varre = variance of random effects
   meanncu = mean of covariate (unexposed)
   meannce = mean of covariate (exposed)
   varnc = variance of covariate
   rhonc = intracluster correlation of covariate
   rrnc = relative risk (covariate)


* check correct proportion of successes by exposure;

proc freq data=test;
	tables exposed*outcome / nocol nopercent;
run;

* check correct proportion of successes by covariate category;

data test;
	set test;
	if exposed = 0 then do;
		if normcov < 0.4 - sqrt(0.25) then catcov = 1;
		if 0.4 - sqrt(0.25) <= normcov < 0.4 then catcov = 2;
		if 0.4 <= normcov < 0.4 + sqrt(0.25) then catcov = 3;
		if normcov >= 0.4 + sqrt(0.25) then catcov = 4;
	end;
	if exposed = 1 then do;
		if normcov < 0.6 - sqrt(0.25) then catcov = 1;
		if 0.6 - sqrt(0.25) <= normcov < 0.6 then catcov = 2;
		if 0.6 <= normcov < 0.6 + sqrt(0.25) then catcov = 3;
		if normcov >= 0.6 + sqrt(0.25) then catcov = 4;
	end;
run;

proc freq data=test;
	tables catcov*exposed*outcome / nocol nopercent;
run;

* check correct number of exposure allocations by simulation;

proc freq data=test;
	tables simulation*exposed / norow nocol nopercent;
run;

* check distribution of covariate, random effects and cluster size;

proc sort data=test out=temp;
	by exposed;
proc univariate data=temp;
	by exposed;
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





