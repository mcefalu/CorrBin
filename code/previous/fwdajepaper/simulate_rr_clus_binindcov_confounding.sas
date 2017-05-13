********************************************************************************************;
*																						   *;
* Macro to Simulate Clustered Data with Confounding by a Binary Individual Level Covariate *;
*																						   *;
********************************************************************************************;



/* data = name of dataset
   nsim = number of simulated datasets contained in data
   nclus = number of clusters
   fixed = whether cluster size is fixed or varies according to a poisson distribution
   nobsclus = (mean) number of observations per cluster
   p0 = expected probability of success in control group (not the same as exp(beta_0) since addition of random effects on the 
		log scale changes the mean on the original scale from exp(beta_0) to exp(beta_0+varre/2) => beta_0=log(p0/exp(varre/2))
   rre = relative risk (exposure)
   varre = variance of random effects
   pbcu = probability of having covariate (unexposed)
   pbce = probability of having covariate (exposed)
   rhobc = intracluster correlation for binary covariate
   rrbc = relative risk (covariate)

   the macro:
   1) determines the size of each cluster (fixed or randomly generated from a poisson distribution)
   2) assigns half of all patients to exposed group and half to unexposed group

   note: ignore meanncu, meannce, varnc, rhonc, rrnc (these are used to simulate more complicated datasets) */

%macro simulate(data, nsim, nclus, fixed, nobsclus, p0, rre, varre, pbcu, pbce, rhobc, rrbc, meanncu, meannce, varnc, rhonc, rrnc);

	data &data.;
		length id 8. simulation cluster_id cluster_size subject_id intercept exposed 3. random 8. bincov outcome 3.;
		id = .;
	  	do simulation = 1 to &nsim.;
	   		do cluster_id = 1 to &nclus.;
				if &fixed. = 1 then cluster_size = &nobsclus.;
				else cluster_size = ranpoi(cluster_id*simulation, &nobsclus.);
				if cluster_size = 0 then cluster_size = 1; 						* need this since clfsim falls over if have any clusters of size 0 when it tries to loop;
				intercept = 1;
				i = 1;
				do until (prob < 1);
	    			i = i + 1; 	* seed starts at 2*cluster_id*simulation;
					random = sqrt(&varre.)*rannor(i*cluster_id*simulation);	
					prob = (&p0./exp(&varre./2))*&rre.*&rrbc.*exp(random);		* exposure and covariate could be 1 for some patients in cluster so make sure prob ok for these patients;
				end;
				do subject_id = 1 to cluster_size;
					exposed = .;
					outcome = .;
					uniform = uniform(cluster_id*subject_id*simulation);		* use this to allocate exposure;
					bincov = .;	
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
		drop prob i uniform nobs rand n;
	run;

	proc sort data=&data.;
		by simulation cluster_id subject_id;
	run;

	data &data.;
		set &data.;
		id = _n_;
		if exposed = 0 then prob_success = &pbcu.;
		if exposed = 1 then prob_success = &pbce.;
		compatible = .;
	run;

	options nonotes;

	proc iml;

		%include "./code/fwdajepaper/clfsim.sas";
		*%include "\\Dmac-fs02\dmc\Stats\Lisa\SAS Macros\clfsim.sas";

	  	rho = &rhobc.;  											
		use &data.;								
		read all var {id simulation cluster_id cluster_size subject_id intercept exposed random bincov outcome prob_success compatible} into &data.;

		do sim = 1 to &nsim.;
			do clus = 1 to &nclus.;
				use &data. where(simulation=sim & cluster_id=clus);								
				read all var {id cluster_size prob_success} into temp;

				first_id = temp[1,1];
			  	cluster_size = temp[1,2];   										* cluster size;

			  	mean = temp[,3]; 													* mean vector;
				if cluster_size = 1 then corr = 1;
			  	else corr = xch(cluster_size, rho);                      			* correlation matrix;
			  	cov = cor2var(corr, mean);         									* covariance matrix;

			  	B = allreg(cov); 													* check for clf compatability;
			  	err = blrchk1 (mean, B);
				if (err) then &data.[first_id:first_id+cluster_size-1,12] = 0;
				else &data.[first_id:first_id+cluster_size-1,12] = 1;

				&data.[first_id:first_id+cluster_size-1,9] = mbsclf1 (mean, B);     * simulate data;
			end;
		end;

	  	create _output from &data.; 												* output data;
	    append from &data.;

	quit;

	options notes;

	data &data.;
		set _output;
		rename col2 = simulation 
			   col3 = cluster_id 
			   col4 = cluster_size 
			   col5 = subject_id 
			   col6 = intercept 
			   col7 = exposed
			   col8 = random 
			   col9 = bincov
			   col10 = outcome
			   col12 = compatible;
		drop col1 col11;
	run;

	data &data.;
		set &data.;
		outcome = ranbin(cluster_id*subject_id*simulation, 1, (&p0./exp(&varre./2))*exp(log(&rre.)*exposed)*exp(log(&rrbc.)*bincov)*exp(random));	
	run;

	proc means data=test noprint;
		var compatible;
		output out=temp min(compatible)=;
	run;

	data temp;
		set temp;
		call symput('continue', compatible);
	run;

	%if &continue. = 0 %then %do; 					* this will delete dataset so no analysis will be done and can just try a different seed;

		proc datasets lib=work nolist;
		  	delete &data.;
		run;
		quit;

	%end;

%mend;

/* test macro;

%let seed = 4966;

%simulate(test, 10, 50, 0, 10, 0.1, 3, 0.1, 0.1, 0.9, 0.05, 2, , , , , );				* compatible = 0 - deletes test as it is supposed to;

%simulate(test, 50, 50, 0, 10, 0.1, 3, 0.1, 0.4, 0.6, 0.05, 2, , , , , );				* compatible = 1 - creates test;

* check correct proportion of successes by exposure and covariate;

proc freq data=test;
	tables bincov*exposed*outcome / nocol nopercent;
run;

* check correct covariate prevalence by exposure;

proc freq data=test;
	tables exposed*bincov / nocol nopercent;
run;

* check correct number of treatment allocations by simulation;

proc freq data=test;
	tables simulation*exposed / norow nocol nopercent;
run;

* check distribution of random effects and cluster size;

data temp;
	set test;
	by simulation cluster_id;
	if first.cluster_id;
run;

proc univariate data=temp;
	var random cluster_size; 		
	histogram random cluster_size;
run;

* check distribution of proportion of successes in each cluster (mean is as expected, variance shows overdispersion due to clustering as expected);

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
run; 

* check correlations for covariate;

proc mixed data=test;
	by simulation;
	class exposed cluster_id;
	model outcome = exposed;
	random cluster_id;
	ods output CovParms = _covparms;
run;

proc transpose data=_covparms out=_corr;
	by simulation;
	var estimate;
	id covparm;
run;

data _corr;
	set _corr;
	icc = cluster_id / (cluster_id + residual);
run;

proc univariate data=_corr;
	var icc;
	histogram icc;
run; */
