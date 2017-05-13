
proc import datafile="./temp/dta&SYSPARM..csv"
  out=Data&SYSPARM.
  dbms=csv
  replace;
  getnames=yes;
run;


ods trace on;
proc genmod data=Data&SYSPARM. descending ;
  class id;
  model Y= c1 c2 / dist=poisson link=log;
  repeated subject=id / type=cs;
  ods output GEEEmpPEst=covariates1&SYSPARM. ;
run;
ods trace off;

proc export data=covariates1&SYSPARM. outfile="./temp/poisson&SYSPARM..csv" dbms=csv replace; run;


ods trace on;
proc genmod data=Data&SYSPARM. descending ;
  class id;
  model Y= c1 c2 / dist=binomial link=log;
  repeated subject=id / type=cs;
  ods output GEEEmpPEst=covariates&SYSPARM. ;
run;
ods trace off;

proc export data=covariates&SYSPARM. outfile="./temp/sasout&SYSPARM..csv" dbms=csv replace; run;

