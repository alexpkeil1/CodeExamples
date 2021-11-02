/**********************************************************************************************************************
* Author: Alex Keil
* Program: weibull_hr_relns.sas
* Date: 8/21/12
* Project: examples
* Tasks: Simple simulation to show relationship between HR and weibull paramaters in AFT model
**********************************************************************************************************************/
*clear the log window and the output window;
DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120 SKIP = 2 FORMDLIM = '-' MPRINT NOCENTER;
OPTIONS FORMCHAR = '|----|+|---+=|-/\<>*';
%LET shape = 2;*1 for exponential, other >0 for weibull;
%LET hr = 5; *hazard ratio;
DATA sim;
	CALL STREAMINIT(1212121);
	DO i = 1 TO 400000;
	shape = &shape; *shape parameter for weibull dist'n (1=exponential);
	scaleint = 10;
	HR = &hr; *hazard ratio for x;

	z = EXP(RAND('normal', 0,1));
	x = EXP(RAND('normal', 0 + .5*log(z) ,1));
	z2 = (RAND('BERNOULLI', 0 + .5 ));
	wscale = exp(scaleint + x*(-log(hr)/(shape)) + 1*z + 0*z2); *non-collapsible if including z2;
	wscale2 = exp(scaleint + x*(-log(hr)/(shape)) + 1*z + log(50)*z2); *non-collapsible if including z2;

	t0 = RAND('weibull', shape, wscale);
	t02 = RAND('weibull', shape, wscale2);
	t = MIN(t0, 5);
	t2 = MIN(t02, 5);
	IF t = 5 THEN d = 0; ELSE d=1;
	IF t2 = 5 THEN d2 = 0; ELSE d2=1;
	OUTPUT;
	END;
RUN;
PROC MEANS DATA = sim FW=5 MAXDEC=2;
RUN;
PROC PHREG DATA = sim;
	TITLE "Cox model, z2-/>T";
	MODEL t*d(0) = x  z/ RL;
	ODS SELECT parameterestimates;
PROC PHREG DATA = sim;
	TITLE "Cox model, z2->T";
	MODEL t2*d2(0) = x  z/ RL;
	ODS SELECT parameterestimates;
PROC LIFEREG DATA = sim;
	*WHERE x2=1;
	TITLE "AFT Weibull model - HR = exp(-betaX)^(Weibull shape), z2-/>T";
	MODEL t*d(0) = x  z;
	ODS SELECT parameterestimates;
	ODS OUTPUT  parameterestimates=wparms;
RUN;
PROC LIFEREG DATA = sim;
	TITLE "AFT Weibull model - HR = exp(-betaX)^(Weibull shape), z2->T";
	MODEL t2*d2(0) = x  z;
	ODS SELECT parameterestimates;
	ODS OUTPUT  parameterestimates=wparms2;
RUN;

DATA wparms (KEEP = label betax wshape hr);
	LENGTH label $30;
	SET wparms END=eof;
	label = "Z2 has no effect on T";
	RETAIN betax wshape HR;
	IF parameter = "x" THEN betax = Estimate;
	IF parameter = "Weibull Shape" THEN wshape = Estimate;
	hr = exp(-1*betax*wshape);
	IF eof THEN OUTPUT;
DATA wparms2 (KEEP = label betax wshape hr);
	LENGTH label $30;
	SET wparms2 END=eof;
	label =  "Z2 affects T";
	RETAIN betax wshape HR;
	IF parameter = "x" THEN betax = Estimate;
	IF parameter = "Weibull Shape" THEN wshape = Estimate;
	hr = exp(-1*betax*wshape);
	IF eof THEN OUTPUT;

DATA parms;
	SET wparms2 wparms;
	betabias = betax-(-log(&hr)/(&shape));
	hrbias = log(hr)-log(&hr);
PROC PRINT DATA = parms;
	TITLE "HR from Weibull model = exp(-beta*x)^weibull shape, assessing collapsibility";
RUN; *HR is not collapsible, beta parameter from AFT is;

RUN;QUIT;RUN;
