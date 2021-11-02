/**********************************************************************************************************************
* Author: Alex Keil
* Program: err_collapsibility.sas
* Date: 11/12/12
* Project: Dissertation
* Tasks: Explore collapsibility of the excess relative risk
* Data in: simulated
* Data out: NA
* Description: 
**********************************************************************************************************************/
*clear the log window and the output window;
DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120 PAGESIZE=70 SKIP = 2 FORMDLIM = '-' MPRINT NOCENTER;
OPTIONS FORMCHAR = '|----|+|---+=|-/\<>*';
FOOTNOTE;TITLE;

%LET ncohorts = 2000;
%LET n = 3000;
%LET baserate = log(68/100000*8); *baseline rate (from report to nation, divided by number to yield efficiency here);
%LET censtime = 80/8;
%LET truebeta = 0.5;*err, radon;
%LET truezbeta = 10;*err, smoking;
%LET alpha = 0.1;

*create data;
DATA a;
	CALL STREAMINIT(327089823);
	DO cohort = 1 TO &ncohorts;
		DO id = 1 TO &n;
			z = RAND("bernoulli", 0.5); *smoking;
			x = RAND("bernoulli", 1/(1+exp(-(0.5 - log(8)*0.5))));
			rate = exp(&baserate)*(1 + &truebeta*x + &truezbeta*z);*err model;
			*rate = exp(&baserate + log(1+&truebeta)*x + log(1+&truezbeta)*z);*log linear rr model;
			ty = RAND("weibull", 1, 1/rate);
			t = min(ty, &censtime);
			y = (ty<=&censtime);
			outint = FLOOR(t-0.000001);
			OUTPUT;
		END;
	END;
RUN;



*break into person periods;
DATA b last;
	SET a;
	BY cohort id;
	DO in = 0 TO outint;
		out = min(in+1, t);
		d = (in<ty<=out)*y;
		py = out-in;
		OUTPUT b;
		IF in=outint THEN OUTPUT last;
	END;
RUN;

*tabulate;
PROC MEANS DATA = b SUM NOPRINT;
	BY cohort;
	CLASS x z;
	VAR py d;
	OUTPUT OUT = btab(WHERE=(_TYPE_=3) RENAME=(_FREQ_=n))  SUM = pyr events;
RUN;

*generate strata weights(real weight function acts on person level, exposure probability constant over time);

DATA btab2;
	SET btab  ;
	BY cohort x z;
	lpy = LOG(pyr);
	*FOR RISK MODELS;
	*risk = events/n;
RUN;


*run model via nlmixed;
*crude - relative rate;
ODS LISTING CLOSE;
OPTIONS NONOTES;
PROC GENMOD DATA=btab2 DESCENDING;
	BY cohort;
	MODEL events = x / DIST=POISSON LINK=LOG OFFSET=LPY;
	ODS OUTPUT parameterestimates = rr;
DATA rrparms(KEEP = cohort b rr se bias cover power true);
	SET rr(KEEP = cohort parameter estimate stderr RENAME=(stderr=se));
	IF parameter = "x" THEN DO;
		b = estimate;
		rr = exp(b);
		true = LOG(1+&truebeta);
		bias = estimate-true;
		cover = estimate-1.96*se <= 1+&truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";
*crude - relative risk;
PROC GENMOD DATA=last DESCENDING;
	BY cohort;
	MODEL D = x / DIST=BINOMIAL LINK=LOG ;
	ODS OUTPUT parameterestimates = rr2;
DATA rrparms2(KEEP = cohort b rr se bias cover power true);
	SET rr2(KEEP = cohort parameter estimate stderr RENAME=(stderr=se));
	IF parameter = "x" THEN DO;
		b = estimate;
		rr = exp(b);
		true = LOG(1+&truebeta);
		bias = estimate-true;
		cover = estimate-1.96*se <= 1+&truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";


*crude - excess relative rate;

PROC NLMIXED DATA=btab2 ;
	BY cohort;
	eta = beta0 ;
	lambda = exp(eta)*(1 + beta2*x);
	MODEL events ~ POISSON(pyr*lambda);
	ODS OUTPUT parameterestimates = nlparmsuw;
DATA parms(KEEP = cohort b se bias cover power true);
	SET nlparmsuw(KEEP = cohort parameter estimate standarderror RENAME=(standarderror=se));
	IF parameter = "beta2" THEN DO;
		b = estimate;
		true = &truebeta;
		bias = estimate-true;
		cover = estimate-1.96*se <= &truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";

*crude - excess relative risk;
PROC NLP DATA=btab2 SE VARDEF=N;
	BY cohort;
	PARMS beta0=-1, beta2=.1;
	eta = beta0 ;
	lambda = exp(eta)*(1 + beta2*x);
	llambda = log(lambda);
	mlambda = log(1-lambda);
	logL = (events*llambda + (n-events)*mlambda); *binomial likelihood;
	MAX logL;
	ODS OUTPUT parameterestimates = nlpparmswrisk;
RUN;
DATA parms2risk(KEEP = cohort b se bias cover power true);
	SET nlpparmswrisk(KEEP = cohort parameter estimate appstderr gradobj RENAME=(appstderr=se) WHERE=(-0.05<gradobj<0.05));
	IF parameter = "beta2" THEN DO;
		b = estimate;
		true = &truebeta;
		bias = estimate-true;
		cover = estimate-1.96*se <= &truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";
*stratification adjusted - relative rate;
PROC GENMOD DATA=btab2 DESCENDING;
	BY cohort;
	MODEL events = x Z/ DIST=POISSON LINK=LOG OFFSET=LPY;
	ODS OUTPUT parameterestimates = rr0;
DATA rrparms0(KEEP = cohort b rr se bias cover power true);
	SET rr0(KEEP = cohort parameter estimate stderr RENAME=(stderr=se));
	IF parameter = "x" THEN DO;
		b = estimate;
		rr=exp(b);
		true = LOG(1+&truebeta);
		bias = estimate-true;
		cover = estimate-1.96*se <= 1+&truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";
*stratification adjusted - relative risk;
PROC GENMOD DATA=last DESCENDING;
	BY cohort;
	MODEL d = x Z/ DIST=BINOMIAL LINK=LOG;
	ODS OUTPUT parameterestimates = rr02;
DATA rrparms02(KEEP = cohort b rr se bias cover power true);
	SET rr02(KEEP = cohort parameter estimate stderr RENAME=(stderr=se));
	IF parameter = "x" THEN DO;
		b = estimate;
		rr=exp(b);
		true = LOG(1+&truebeta);
		bias = estimate-true;
		cover = estimate-1.96*se <= 1+&truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";
*stratification adjusted - excess relative rate;
PROC NLMIXED DATA=btab2 ;
	BY cohort;
	eta = beta0 ;
	lambda = exp(eta)*(1 + beta1*z + beta2*x);
	MODEL events ~ POISSON(pyr*lambda);
	ODS OUTPUT parameterestimates = nlparmss;
DATA parms0(KEEP = cohort b se bias cover power true);
	SET nlparmss(KEEP = cohort parameter estimate standarderror RENAME=(standarderror=se));
	IF parameter = "beta2" THEN DO;
		b = estimate;
		true = &truebeta;
		bias = estimate-true;
		cover = estimate-1.96*se <= &truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";
*stratification adjusted - excess relative risk using NLP;
PROC NLP DATA=btab2 SE VARDEF=N;
	BY cohort;
	PARMS beta0=-1, beta2=.1, beta3=.1;
	eta = beta0 ;
	lambda = exp(eta)*(1 + beta2*x + beta3*z);
	llambda = log(lambda);
	mlambda = log(1-lambda);
	logL = (events*llambda + (n-events)*mlambda); *binomial log-likelihood;
	MAX logL;
	ODS OUTPUT parameterestimates = nlpparmswrisk3;
RUN;
DATA parms3risk(KEEP = cohort b se bias cover power true);
	SET nlpparmswrisk3(KEEP = cohort parameter estimate appstderr gradobj RENAME=(appstderr=se) WHERE=(-0.05<gradobj<0.05));
	IF parameter = "beta2" THEN DO;
		b = estimate;
		true = &truebeta;
		bias = estimate-true;
		cover = estimate-1.96*se <= &truebeta <= estimate+1.96*se;
		power = (estimate/se)**2>(1.96)**2;
		OUTPUT;
	END;
	LABEL se = " ";
RUN;
OPTIONS NOTES;
ODS LISTING ;
DM ODSRESULTS 'clear;' CONTINUE; *clear ODS generated datasets (clears all results);

PROC MEANS DATA = a FW=5 MAXDEC=2 MEAN MEDIAN P5 P95 STD;
	TITLE "Sample data";
RUN;

PROC MEANS DATA = rrparms MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "GENMOD -  crude RR (rate) estimates for x";
	TITLE2 "Bias term measures non-collapsibility + improper model specification";
	VAR rr b bias se cover power true;
RUN;
PROC MEANS DATA = rrparms0 MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "GENMOD -  stratified RR (rate) estimates for x";
	VAR rr b bias se cover power true;
RUN;
PROC MEANS DATA = rrparms2 MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "GENMOD -  crude RR (risk) estimates for x";
	TITLE2 "Bias term measures non-collapsibility + improper model specification";
	VAR rr b bias se cover power true;
RUN;
PROC MEANS DATA = rrparms02 MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "GENMOD -  stratified RR (risk) estimates for x";
	VAR rr b bias se cover power true;
RUN;
PROC MEANS DATA = parms MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "NLMIXED -  crude ERR (rate) estimates for x";
	TITLE2 "Bias term measures non-collapsibility";
	VAR b bias se cover power true;
RUN;
PROC MEANS DATA = parms0 MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "NLMIXED - stratification adjusted ERR (rate) estimates for x";
	VAR b bias se cover power true;
RUN;

PROC MEANS DATA = parms2risk MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "NLP - crude ERR (risk) estimates for x ";
	VAR b bias se cover power true;
RUN;

PROC MEANS DATA = parms3risk MAXDEC=2 FW=5 MEAN MEDIAN STD P5 P95 N;
	TITLE "NLP - stratification adjusted ERR (risk) estimates for x";
	VAR b bias se cover power true;
RUN;

RUN;QUIT;RUN;
