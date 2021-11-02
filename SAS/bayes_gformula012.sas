/*DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; */
/**********************************************************************************************************************
* Author: Alex Keil
* Program: bayes_gformula012.sas
* System notes: SAS 9.4 (SAS/STAT 13.1), Windows 8.1 [does not work in lower versions of SAS due to MCMC changes]
* Date: Wednesday, July 2, 2014 at 11:41:42 AM
* Project: miscellaneous
* Tasks: Implement bayesian g-formula algorithm in PROC MCMC
* Data in: NA
* Data out: NA
* Description: program simulates data for a small longitudinal analysis
   This takes M cohorts (1000 seems to work fine) of N simulated individuals and performs
      1) frequentist g-formula
         uses bootstrap analysis using percentile based confidence intervals from 1,000 samples
      2) bayesian g-formula under vague null centered priors for every coefficient
         uses MCMC methods
         using percentile based confidence intervals from posterior samples 
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
*UPDATES: calculates non-bootstrap RD for frequentist g-formula for comparison
* Changes from previous version: Includes an L1 loss function
* Next steps:

**********************************************************************************************************************/
*clear the log window and the output window;
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120  PAGESIZE=80 SKIP = 2 FORMDLIM = '-' MPRINT NOCENTER;
OPTIONS FORMCHAR = '|----|+|---+=|-/\<>*';



%LET n=40; *sample size - does Bayes g-formula outperform frequentist g-formula in small samples?;
%LET m=1000; *number of different cohorts (500 seems to be the minimum needed);
%LET MC=2000; *number of replicates for MC approximation of g-formula;
%LET nboots = 1000; *number of bootstrap samples;
%LET trueRD = .2;

OPTIONS NONOTES;
*0) simulate data;
DATA d (keep = cohort id X: L: Y: py00);
LENGTH cohort id x1 l2 x2 y 3;
*CALL STREAMINIT(1192887);
DO cohort = 1 TO &m;
*observed data;
 DO id = 1 TO &N;
  x1 = RAND("bernoulli", 0.5);
  py00 = RAND("uniform")*0.1 + 0.4;  
  l2 = RAND("bernoulli", 1/(1+exp(-1 + x1 + py00)));
  x2 = RAND("bernoulli", 1/(1+exp(-1 + x1 + l2)));
  py = py00 + &trueRD*((x1 + x2)/2); *true risk difference per unit exposure;
  y = RAND("bernoulli", py);
  OUTPUT;
 END;
 END;
RUN;

*0.1) Example data;
PROC PRINT DATA = d (OBS=10); TITLE "Example data";
PROC MEANS DATA = d MAXDEC=3 FW=5; TITLE "Observed covariate distributions"; VAR x1 x2 l2 y py00; RUN;
************************************************************************ 
 1) frequentist g-formula;
 ************************************************************************;
*bootstrap sample (1000 samples);
PROC SURVEYSELECT DATA=d OUT=d2(DROP=NumberHits ExpectedHits SamplingWeight) METHOD=URS N=%SYSEVALF(&n*&nboots) NOPRINT OUTHITS; STRATA cohort;RUN;
DATA d2;
 SET d2;
  CALL STREAMINIT(12132);
  randord = RAND('uniform');

PROC SORT DATA=d2; BY cohort randord;
DATA d2; LENGTH cohort set id count 3; SET D2; BY cohort;
 RETAIN set count;
 IF _n_=1 THEN set=0;
 IF first.cohort THEN DO; set=1; count=1; END;
  IF count>&n THEN DO;
   set=set+1;
   count=1;
  END;
  OUTPUT;
  count=count+1;
run;

PROC MEANS DATA = d2 MEAN NOPRINT;
 CLASS cohort set;
 VAR id;
 OUTPUT OUT = meanids MEAN=meanid;
PROC MEANS DATA = meanids;
 TITLE "Checking that IDs are more or less evenly distributed b/w bootstrap samples";
 VAR meanid;
run;

*original (non bootstrap) parameter estimates;
PROC LOGISTIC NOPRINT DATA = d DESCENDING OUTEST=lm (KEEP=cohort intercept x1 RENAME=(intercept=a0 x1=a1));
 BY cohort;
 MODEL l2 = x1;
PROC LOGISTIC NOPRINT DATA = d DESCENDING OUTEST=ym (KEEP=cohort intercept x1 x2 l2 RENAME=(intercept=b0 x1=b1 x2=b2 l2=b3));
 BY cohort;
 MODEL y = x1 x2 l2;
RUN;
DATA e;
 MERGE d ym lm;
 BY cohort ;
RUN;
*mc sample for g-formula (non bootstrap);
PROC SURVEYSELECT DATA = e OUT=interv(KEEP=cohort a0 a1 b0 b1 b2 b3) N=&mc SEED=1223 METHOD=URS OUTHITS NOPRINT;  STRATA cohort; run;
PROC DATASETS LIBRARY=work; DELETE e MEANIDS;QUIT;
DATA interv(KEEP=cohort rdi spe bias);
 *CALL STREAMINIT(121232);
  SET interv;
    *py0 = 1/(1+exp(-(b0 + (b1+b2)*0 + b3*1/(1+exp(-(a0 + a1*0))))));
    *rdi = 1/(1+exp(-(b0 + (b1+b2)*1 + b3*1/(1+exp(-(a0 + a1*1))))))-py0;
    rdi = 1/(1+exp(-(b0 + (b1+b2)*1 + b3*1/(1+exp(-(a0 + a1*1))))))-
          1/(1+exp(-(b0 + (b1+b2)*0 + b3*1/(1+exp(-(a0 + a1*0))))));
    bias = rdi-&truerd;
    spe = (bias)**2;
  LABEL rdi = "Risk difference (frequentist g-formula)"
       /* py0 = "Potential risk given no exposure"*/;
/*PROC SORT DATA = interv; BY cohort; */
PROC MEANS DATA = interv MEAN NOPRINT;
  BY cohort;
  VAR  /*py0*/ rdi spe bias;
  OUTPUT OUT=fr0 MEAN= /*py0*/ rd mspe bias;
RUN;
   
 PROC DATASETS LIBRARY=work; DELETE interv;QUIT;


*bootstrap estimates;
PROC LOGISTIC NOPRINT DATA = d2 DESCENDING OUTEST=lm (KEEP=cohort set intercept x1 RENAME=(intercept=a0 x1=a1));
 BY cohort set;
 MODEL l2 = x1;
PROC LOGISTIC NOPRINT DATA = d2 DESCENDING OUTEST=ym (KEEP=cohort set intercept x1 x2 l2 RENAME=(intercept=b0 x1=b1 x2=b2 l2=b3));
 BY cohort set;
 MODEL y = x1 x2 l2;
RUN;
DATA e;
 MERGE d2 lm ym;
 BY cohort set;
RUN;
PROC DATASETS LIBRARY=work; DELETE lm ym d2;QUIT;

*mc sample for g-formula (bootstrap samples);
PROC SURVEYSELECT DATA = e OUT=interv(KEEP=cohort set a0 a1 b0 b1 b2 b3) N=&mc SEED=1223 METHOD=URS OUTHITS NOPRINT;  STRATA cohort set; run;
PROC DATASETS LIBRARY=work; DELETE e ;QUIT;
/********** trouble is here *************/
DATA interv(KEEP=cohort set rdi );
 *CALL STREAMINIT(121232);
  SET interv;
    *py0 = 1/(1+exp(-(b0 + (b1+b2)*0 + b3*1/(1+exp(-(a0 + a1*0))))));
    *rdi = 1/(1+exp(-(b0 + (b1+b2)*1 + b3*1/(1+exp(-(a0 + a1*1))))))-py0;
    rdi = 1/(1+exp(-(b0 + (b1+b2)*1 + b3*1/(1+exp(-(a0 + a1*1))))))-
          1/(1+exp(-(b0 + (b1+b2)*0 + b3*1/(1+exp(-(a0 + a1*0))))));
RUN;
DATA interv(KEEP=cohort set rdi ape spe bias);
 *CALL STREAMINIT(121232);
  SET interv;
    bias = rdi-&truerd;
    spe = (bias)**2;
    ape = abs(bias);
   LABEL rdi = "Risk difference (frequentist g-formula)"
        /*py0 = "Potential risk given no exposure"*/;
/*PROC SORT DATA = interv; BY cohort; */
PROC MEANS DATA = interv MEAN NOPRINT;
  BY cohort set;
  VAR /*py0*/ rdi spe ape bias;
  OUTPUT OUT=fr(DROP=_:) MEAN=  /*py0*/ rd mspe mape bias;
RUN;
PROC DATASETS LIBRARY=work; DELETE interv ;QUIT;
PROC MEANS DATA = fr MEAN STD NOPRINT;
  BY cohort;
  VAR rd /*py0*/ bias mape mspe;
  OUTPUT out = fr2(DROP=_:) MEAN=rd /*py0*/ bias mape mspe STD=mcse /*spy0*/ sbias smape smspe;
RUN;


PROC UNIVARIATE DATA = fr NOPRINT;
  BY cohort;
  VAR rd;
  OUTPUT out = fr3 PCTLPRE=prd_  PCTLPTS=2.5 97.5;
RUN;

DATA fr2;
 MERGE fr2 fr3;
 BY cohort;
 mse = MCSE*MCSE + BIAS*BIAS;
 mcvar = mcse*mcse;
 cover = rd-mcse*1.96 < &truerd < rd+mcse*1.96;
 cover_p = prd_2_5 < &truerd < prd_97_5;
RUN;
 
************************************************************************
2) Bayes g-formula
************************************************************************;
DATA d;
 SET d;
 BY cohort ;
 IF first.cohort THEN ind=0;
 ind+1;
RUN;
ODS LISTING CLOSE;
PROC MCMC DATA = d/* SEED = 12123*/ NMC = 10000 OUTPOST=dpost 

MONITOR=(rd bias mape mspe);

 BY cohort;
 TITLE "Bayesian g-formula risk difference (truth=&truerd)";
 *ODS SELECT PostSumInt ESS TADPanel;
*priors;
 PARMS a0 a1 b0 b1 b2 b3;
 PRIOR a1 b1 b2 b3 ~ NORMAL(0, var=3); PRIOR a0 b0~ NORMAL(LOG(0.5), var=1000);

*joint model for observed data;
  pl2 = LOGISTIC(a0+a1*x1);
  MODEL L2 ~ BINARY(pl2);
  py = LOGISTIC(b0 + b1*x1 + b2*x2 + b3*l2);
  MODEL y ~ BINARY(py);

*posterior predictive risk difference(efficient way);
  ARRAY rdi[&N];
  ARRAY spi[&N];
  ARRAY api[&N];
  rdi[IND] = 1/(1+exp(-(b0 + (b1+b2)*1 + b3*1/(1+exp(-(a0 + a1*1))))))-
               1/(1+exp(-(b0 + (b1+b2)*0 + b3*1/(1+exp(-(a0 + a1*0))))));
   spi[IND] = (rdi[ind]-&truerd)**2;
  api[IND] = abs(rdi[ind]-&truerd);
  IF ind=&n THEN DO;
   mspe =  MEAN(of spi:);
   mape =  MEAN(of api:);
   rd =  MEAN(of rdi:);
   bias = rd-&trueRD;
  END;

  LABEL rd="Average difference of potential outcome ";
RUN;
ODS LISTING;
PROC MEANS DATA = dpost NOPRINT;
 BY cohort;
 VAR rd bias mspe mape;
 OUTPUT OUT = bayesres MEAN= STD= / AUTONAME;
RUN;
PROC UNIVARIATE DATA = dpost NOPRINT;
 BY cohort;
 VAR rd;
 OUTPUT OUT = bayesres2 PCTLPRE=prd_ PCTLPTS=2.5 97.5;
RUN;
PROC DATASETS LIBRARY=work; DELETE d dpost;QUIT;


DATA br(KEEP=rd  bias mcse mcvar mse cover cover_p mspe mape); 
 MERGE bayesres bayesres2;
 BY cohort;
 rd = rd_Mean;
 bias = bias_mean;
 mcse = rd_stddev;
 mcvar = rd_stddev*rd_stddev;
 mse = bias*bias + mcse*mcse;
 cover = rd-1.96*mcse < &truerd < rd+1.96*mcse;
 cover_p = prd_2_5 < &truerd < prd_97_5;
 mspe = mspe_mean;
 mape = mape_mean;

RUN;
OPTIONS NOTES;

PROC MEANS DATA = fr0  FW=6 MEAN STD VAR MEDIAN P5 P95 MIN MAX NOLABELS;
  TITLE "Frequentist g-formula risk difference (original data) (truth=&truerd)";
  VAR rd /*py0*/  ;
RUN;
PROC MEANS DATA = fr2  FW=6 MEAN STD VAR MEDIAN P5 P95 MIN MAX  NOLABELS;
  TITLE "Frequentist g-formula risk difference (bootstrap estimates) (truth=&truerd)";
  VAR rd /*py0*/ bias mcse mcvar mse mspe mape cover cover_p;
RUN;
PROC MEANS DATA = br  FW=6 MEAN STD VAR MEDIAN P5 P95 MIN MAX  NOLABELS;
 TITLE "Bayesian g-formula results";
  VAR rd bias mcse mcvar mse mspe mape cover cover_p;
RUN;


