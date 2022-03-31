/**********************************************************************************************************************
* Author: Alex Keil
* Program: pop_attributable_fractions.sas
* System notes: SAS 9.4 (SAS/STAT 13.1), Windows 11
* Date: March 31, 2022
* Tasks: valid calculation of population attributable fractions in uncensored data, under confounding
* Description: program simulates data for a small cross-sectional/single time point analysis
   Using a potential outcome defintion of population attributable fractions (PAF), it calculates a PAF 
   as:
     PAF = (r-r0)/r
   where:
     r0 = average risk of disease, had everyone been unexposed
     r  = average risk of disease in the population
    
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
* Next steps: this yields a point estimate only. confidence intervals can be obtained via bootstrapping

**********************************************************************************************************************/

* simulate some data for an example (use your real data);
DATA original_data;
  CALL STREAMINIT(123211);
  DO i = 1 to 1000;
    exposure = RAND("bernoulli", 0.5);
	confounder = RAND("BERNOULLI", 0.5);
	outcome = RAND("BERNOULLI", 0.2 + 0.2*exposure + 0.2*confounder);
   OUTPUT;
  END;
RUN;

DATA predict_data; /* should be 3x the sample size */
   SET original_data;
   par_group = "Original data";
   OUTPUT;
   par_group = "If unexposed";
   exposure = 0;
   outcome = .;
   OUTPUT;
   par_group = "If exposed"; /* mainly included for reference */
   exposure = 1;
   outcome = .;
   OUTPUT;

* model on simulated data (use your real model);
PROC GENMOD DATA = predict_data DESC;
  MODEL outcome = exposure confounder / LINK=LOG D=BINOMIAL; /*POISSON COULD BE USED, but check that predictions are < 1.0 - can also use LOGIT link*/
  OUTPUT OUT=par_data p=predicted;
  ESTIMATE "rr" exposure 1 / EXP;
RUN;


* PAR/PAF calculations;
PROC MEANS DATA = par_data MEAN;
  CLASS par_group;
  VAR outcome predicted exposure;
  OUTPUT OUT = mns1 MEAN=risk_obs risk_potential p_exposure;
RUN;
PROC MEANS DATA = par_data MEAN;
  WHERE outcome=1;
  VAR exposure;
  OUTPUT OUT = mns2 MEAN= p_exposure_cases;
RUN;
DATA mns;
 SET mns1;
 IF _n_ = 1 THEN SET mns2;

DATA par (KEEP = risk_population risk_unexposed risk_exposed rr rd prob_exposure prob_exposure_cases par );
 SET mns  END=EOF;
 RETAIN risk_population risk_unexposed risk_exposed par prob_exposure prob_exposure_cases 0;
 IF par_group = "Original data" THEN DO; risk_population = risk_potential; prob_exposure = p_exposure; prob_exposure_cases=p_exposure_cases; END;
 ELSE IF par_group = "If unexposed" THEN risk_unexposed = risk_potential;
 ELSE IF par_group = "If exposed" THEN risk_exposed = risk_potential;
 par = (risk_population - risk_unexposed)/risk_population;
 rd = risk_exposed - risk_unexposed;
 rr = risk_exposed/risk_unexposed;
 LABEL par = "Population attributable fraction";
 IF EOF then output;

PROC PRINT DATA = par;
 TITLE "Population/potential risk measures";
RUN;
   
