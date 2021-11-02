DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE;
/**********************************************************************************************************************
* Author: Alex Keil
* Program: proc_mcmc_misclass_correction_bayes_2014223.sas
* Language: R
* Date: Sunday, February 23, 2014 at 12:40:35 AM
* Project:  Fit maclehose et al 2009 models 2 and 3
* Tasks:
* Data in: 
* Data out: 
* Description:
* Keywords:
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/
* 3) simulate data according to simplified version of data from MacLehose et al 2009;
*sample size;
%LET n=1000;
*misclassification parameters;
%LET sens0=.82;
%LET fpr0=.08;
%LET sens1=.90;
%LET fpr1=.03;
*true log odds ratio;
%LET lnor=0.301;
*seed value;
%LET seed = 129312;
*seed value for repeatability;

DATA d;
 CALL STREAMINIT(&seed); *set seed;
 DO i = 1 TO &n;
  z = RAND('bernoulli', 0.5);
  xtrue = RAND('bernoulli',  1/(1+exp(-(-1.5 + 3*z))));
  y = RAND('bernoulli', 1/(1+exp(-(-1 + z + &lnor*xtrue))));
  *measured x;
  IF xtrue = 1 AND y = 1 THEN x = RAND('bernoulli', &sens1);
  ELSE IF xtrue = 0 AND y = 1 THEN x = RAND('bernoulli', &fpr1);
  ELSE IF xtrue = 1 AND y = 0 THEN x = RAND('bernoulli', &sens0);
  ELSE IF xtrue = 0 AND y = 0 THEN x = RAND('bernoulli', &fpr0);
  OUTPUT;
 END;
 

 *tabulate data;
PROC FREQ DATA = d;
 TABLES x*xtrue*y;
RUN;

*frequentist models to confirm data generation step;
*exposure model;
PROC LOGISTIC DATA = d DESCENDING;
 TITLE "true exposure vs. z - true propensity score model";
 TITLE2 "These should approximately equal g0 and g2 from Bayes model";
 MODEL xtrue = z;
 ODS SELECT parameterestimates;
PROC LOGISTIC DATA = d DESC;
 TITLE "measured exposure vs. z - observed propensity score model";
 MODEL x = z;
 ODS SELECT parameterestimates;
RUN;
PROC LOGISTIC DATA = d DESC;
 TITLE "y vs. true exposure, z - true outcome model";
 TITLE2 "These should approximately equal b0, b1 and b2 from Bayes model";
 MODEL y = xtrue z;
 ODS SELECT parameterestimates;

PROC LOGISTIC DATA = d DESC;
 TITLE "y vs. measured exposure, z - this is what we observe";
 MODEL y = x z;
 ODS SELECT parameterestimates;
RUN;

*add in missing values for all values that will be imputed;
*also prepare for estimating risk difference from logistic model;
DATA d;
SET d;
 ind = _N_;
 *corrected exposure value;
 x_cor=.;
 *risk under exposed, unexposed;
 y_1 = .;
 y_0= .;
RUN;
*automatic way for getting sample size (needed for risk difference calculation);
PROC MEANS DATA = d NOPRINT;
 VAR ind  ; OUTPUT OUT = ndata N=sampsize;
RUN;
DATA _NULL_; SET ndata;
 CALL SYMPUT('n', PUT(sampsize, best9.));
RUN;

*Bayesian model to correct for exposure misclassification (will work in sas 9.4, not earlier);
ODS GRAPHICS OFF; *leave on to inspect trace plots;
PROC MCMC DATA = d SEED = 12123 NMC = 100000 /*generally, a lot of iterations are needed for models like this */
 OUTPOST=dpost MONITOR=(int b1 b2 g2 or rr rd);
 TITLE "Bayesian model corrected for exposure misclassification";
 TITLE2 "True b1 = &lnor";
 PARMS int 0 b1 b2 g0 g2;
 PRIOR int ~ GENERAL(0); *diffuse uniform prior (requires int to have initial value in parms statement);
 PRIOR g: b: ~ NORMAL(0, var=5);*note that BUGS uses 1/var; *what happens when you vary this?;
 *model 2 setup from MacLehose et al 2009 - misclassification known with certainty (useful for sensitivity analysis);
 /*BEGINCNST;
  a1 = 0.82;a2 = 0.08;a3 = 0.90;a4 = 0.03;
 ENDCNST;*/
 /**/
*model 3 setup from MacLehose et al 2009;
 PARMS a1 a2 a3 a4 ;
 PRIOR a1 ~ BETA(82, 18);*control sensitivity (true positives, false negatives) ;
 PRIOR a2 ~ BETA(8, 92); *control false positive rate (false positives,true negatives) FPR= 1-sp;
 PRIOR a3 ~ BETA(90, 10);*case sensitivity (true positives, false negatives);
 PRIOR a4 ~ BETA(3, 97); *case false positive rate (false positives,true negatives);

 *true exposure propensity model;
  p_score = LOGISTIC(g0 + g2*z);
  MODEL x_cor ~ BINARY(p_score);

 *outcome model for the log odds;
  py = LOGISTIC(int + b1*x_cor + b2*z);
  MODEL y ~ BINARY(prob=py);

  *observed exposure model - based on sensitivity and specificity priors and imputed true exposure;
  px = a1*(x_cor)*(1-y) + a2*(1-x_cor)*(1-y) + a3*(x_cor)*(y) + a4*(1-x_cor)*y;
  MODEL x ~ BINARY(prob=px);

  *odds ratio;
  or = EXP(b1);

  *model standardized risk difference, ratio (if using cohort data);
  *ARRAY r0i[&N];
  *ARRAY r1i[&N];
  *r1i[IND] = 1/(1+exp(-(int + b1*1 + b2*z)));
  *r0i[IND] = 1/(1+exp(-(int + b1*0 + b2*z)));
  *IF ind=&n THEN DO;
   *rd =  MEAN(of r1i:) - MEAN(of r0i:);
   *rr =  MEAN(of r1i:) / MEAN(of r0i:);
  *END;
RUN;
ODS GRAPHICS;


*posterior median (needed for RR and OR, but not for log(OR) or RD);
PROC UNIVARIATE DATA = dpost NOPRINT;
 VAR rr or b1 rd;
 OUTPUT OUT = posterior_dists PCTLPRE=rr_ or_ b1_ rd_ PCTLPTS = 50 2.5 97.5;
RUN;
PROC TRANSPOSE DATA = posterior_dists  OUT=posterior_dists;RUN;
PROC SORT DATA = posterior_dists; BY _NAME_;

PROC PRINT DATA = posterior_dists;
 TITLE 'Posterior median, 95% percentile based credible interval';
RUN;
