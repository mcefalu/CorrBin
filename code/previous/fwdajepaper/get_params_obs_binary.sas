PROC DATASETS LIB=WORK MEMTYPE=DATA
   kill;
RUN;

%let seed=1;


* import parameter values for this simulation;
proc import datafile="./temp/params&SYSPARM..csv"
  out=mydata
  dbms=csv
  replace;
  getnames=yes;
run;

* sets parameter values as macro variables;
%MACRO Get_data(myDataset,myLine,myColumn,myMVar);
%GLOBAL &myMVar.;
data _null_;
set &myDataset.;
if _N_ = &myLine. then do;
call symput(symget('myMVar'),&myColumn.);
end;
run;
%MEND Get_data;

%Get_data(myDataset=mydata,myLine=1,myColumn=M,myMVar=M);
%Get_data(myDataset=mydata,myLine=1,myColumn=n,myMVar=n);
%Get_data(myDataset=mydata,myLine=1,myColumn=k,myMVar=k);
%Get_data(myDataset=mydata,myLine=1,myColumn=p0,myMVar=p0);
%Get_data(myDataset=mydata,myLine=1,myColumn=rrtrt,myMVar=rrtrt);
%Get_data(myDataset=mydata,myLine=1,myColumn=varre,myMVar=varre);
%Get_data(myDataset=mydata,myLine=1,myColumn=pbce,myMVar=pbce);
%Get_data(myDataset=mydata,myLine=1,myColumn=pbcu,myMVar=pbcu);
%Get_data(myDataset=mydata,myLine=1,myColumn=rhocov,myMVar=rhocov);
%Get_data(myDataset=mydata,myLine=1,myColumn=rrcov,myMVar=rrcov);

%include "./code/fwdajepaper/simulate_rr_clus_binindcov_confounding.sas";

*simulate(data, nsim, nclus, fixed, nobsclus, p0, rre       , varre, pbcu,      pbce, rhobc,           rrbc , meanncu     , meannce , varnc   , rhonc, rrnc);
%simulate(test, &M., &n., 1,       &k.,   &p0., &rrtrt. , &varre., &pbcu.   , &pbce. , &rhocov.    , &rrcov.     ,  ,    , , , );

proc export data=WORK.test outfile="./temp/fullDTA&SYSPARM..csv" dbms=csv replace; run;
