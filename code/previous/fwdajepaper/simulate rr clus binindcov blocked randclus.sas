*******************************************************************************************************************;
*																						 	 				   	  *;
* Macro to Simulate Clustered Data with Blocked Randomisation of Clusters and a Binary Individual Level Covariate *;
*																						 	 				   	  *;
*******************************************************************************************************************;



/* data = name of dataset
   nsim = number of simulated datasets contained in data
   nclus = number of clusters
   fixed = whether cluster size is fixed or varies according to a poisson distribution
   nobsclus = (mean) number of observations per cluster
   p0 = expected probability of success in control group (not the same as exp(beta_0) since addition of random effects on the 
		log scale changes the mean on the original scale from exp(beta_0) to exp(beta_0+varre/2) => beta_0=log(p0/exp(varre/2))
   rrt = relative risk (treatment)
   varre = variance of random effects
   pbc = probability of having covariate
   rhobc = intracluster correlation for binary covariate
   rrbc = relative risk (covariate)

   the macro:
   1) determines the size of each cluster (fixed or randomly generated from a poisson distribution)
   2) assigns half of all clusters to intervention group and half to control group (this would occur at least 
      approximately in practice if randomly permuted blocks were used)
   3) generates a probability of having the covariate for each cluster (pbci) from a beta distribution with mean pbc
   4) randomly allocates patients to have covariate with probability pbci
   5) generates a random effect for each cluster 
   6) resamples from random effect distribution until probability < 1 if necessary
   7) randomly allocates patients to have the outcome with probability (p0/exp(varre/2))*exp(log(rrt)*treatment)*exp(log(rrbc)*bincov)*exp(random)*

   note: ignore meannc, varnc, rhonc, rrnc (these are used to simulate more complicated datasets) */

%macro simulate(data, nsim, nclus, fixed, nobsclus, p0, rrt, varre, pbc, rhobc, rrbc, meannc, varnc, rhonc, rrnc);

	data &data.;
		length simulation cluster_id cluster_size subject_id intercept treatment random bincov outcome 3.;
	  	do simulation = 1 to &nsim.;
	   		do cluster_id = 1 to &nclus.;
				if &fixed. = 1 then cluster_size = &nobsclus.;
				else cluster_size = ranpoi(cluster_id*simulation, &nobsclus.);
				intercept = 1;
				if cluster_id <= &nclus./2 then treatment = 0;
				else treatment = 1;
   				call streaminit(2*cluster_id*simulation);
				pbci = rand('beta', ((1/&rhobc.)-1)*&pbc., ((1/&rhobc.)-1)*(1-&pbc.));	
				i = 2;
				do until (prob < 1);
	    			i = i + 1; 	* seed starts at 3*cluster_id*simulation;
					random = sqrt(&varre.)*rannor(i*cluster_id*simulation);	
					prob = (&p0./exp(&varre./2))*exp(log(&rrt.)*treatment)*&rrbc.*exp(random); 	* bincov could be 1 for some patients in cluster so make sure prob ok for these patients;
				end;
				do subject_id = 1 to cluster_size;
					bincov = ranbin(2*cluster_id*subject_id*simulation, 1, pbci);
					outcome = ranbin(cluster_id*subject_id*simulation, 1, (&p0./exp(&varre./2))*exp(log(&rrt.)*treatment)*exp(log(&rrbc.)*bincov)*exp(random));
	    			output;
				end;	
	 		end;
		end; 
		drop pbci prob i; 
	run;

%mend;


/* test macro;

%simulate(test, 1000, 50, 0, 10, 0.1, 2, 0.1, 0.75, 0.2, 2, , , , );

* check correct proportion of successes by treatment group;

proc freq data=test;
	tables bincov*treatment*outcome / nocol nopercent;
run;

* check correct number of treatment allocations by simulation;

data temp;
	set test;
	by simulation cluster_id;
	if first.cluster_id;
run;

proc freq data=temp;
	tables simulation*treatment / norow nocol nopercent;
run;

* check distribution of random effects, cluster size and cluster specific covariate probabilities;

data temp;
	set test;
	by simulation cluster_id;
	if first.cluster_id;
run;

proc univariate data=temp;
	var random cluster_size pbci; 		
	histogram random cluster_size pbci;
run;

data random;
   	do i=1 to 1000;
   		call streaminit(2*i);
      	*x1=rand('beta',9.5,9.5); 			* icc=0.05, pbc=0.5;
      	*x1=rand('beta',14.25,4.75);		* icc=0.05, pbc=0.75;
      	*x1=rand('beta',2,2); 				* icc=0.2, pbc=0.5;
      	x1=rand('beta',3,1); 				* icc=0.2, pbc=0.75;
      	output; 
   end; 
run;

proc univariate data=random noprint;
	histogram x1;
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
run; */



/* check worst case truncation point for random effects;

%simulate(test, 1000, 50, 0, 10, 0.1, 2, 0.2, 0.5, 0.2, 2, , , , );

data temp;
	set test;
	if treatment = 1 and bincov = 1;
run;

proc univariate data=temp;
	var random; 		
	histogram random;
run; */
