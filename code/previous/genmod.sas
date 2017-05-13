proc import datafile="./temp/dta&SYSPARM..csv"
  out=Data
  dbms=csv
  replace;
  getnames=yes;
run;


ods trace on;
proc genmod data=Data descending ;
  class id;
  model Y= c1 c2 c3 c4 c5 / dist=poisson link=log;
  repeated subject=id / type=cs;
  ods output GEEEmpPEst=covariates1 ;
run;
ods trace off;

proc export data=covariates1 outfile="./temp/poisson&SYSPARM..csv" dbms=csv replace; run;
