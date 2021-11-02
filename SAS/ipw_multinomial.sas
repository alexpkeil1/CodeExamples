*checking param=ref vs. param=glm in proc logistic with multinomial exposure in IPW;


%LET dosecats=4; *number of categories of exposure;


%MACRO dum();
%GLOBAL pxl;
%LET i=2;
%LET pxL=_px_[1];
%DO %WHILE (%EVAL(&i<=&dosecats));
 %LET pxL=&pxL , _px_[&i];
 %LET i = %EVAL(&i+1);
%END;
%MEND;
%dum;
DATA _null_;
 test = "&pxl";
 PUT test;
RUN;

DATA a (KEEP=id j z v x y) ps(KEEP=_:);
 DO j = 1 TO 10;
  DO id  = 1 TO 10000;
   	  z = RAND('bernoulli', 0.5);
   	  v = RAND('table', 0.3, 0.3, 0.4);
      ARRAY _cpx_[&dosecats]; *probability that exposure category is <= i;
      ARRAY _px_[&dosecats]; *probability that exposure category = i;
	  DO p = 1 to &dosecats;
	   IF p>1 THEN _cpx_[p] = _cpx_[p-1] + 1/(1+exp(-(-2.75 + p/&dosecats + 2*z + .5*v - .5*z*v)));
	   ELSE _cpx_[p] = 1/(1+exp(-(-2.75 + p/&dosecats + 2*z + .5*v - .5*z*v)));
	  END;
      DO i = &dosecats TO 2 BY -1; 
        _px_[i] =  (_cpx_[i]-_cpx_[i-1])/MAX(of _cpx_[*]);  
	  END;
      _px_[1] = _cpx_[1]/MAX(of _cpx_[*]);
      *_px_[&dosecats] = 1-SUM(of _px_[*]);
	  x = RAND('table', &pxl);*multinomial distribution;
	  y = RAND('bernoulli', 1/(1+exp(-(-5 + x + 3*z + 2*v + -4*v*z))));
	  OUTPUT a;
	  OUTPUT ps;
  END;
 END;
RUN;

PROC FREQ DATA = a;
 TABLES x*v / MISSING;
RUN;

PROC SORT DATA = a;
 BY j x;
RUN;

ODS LISTING CLOSE;
PROC LOGISTIC DATA = a OUT=nmod;
 BY j;
 TITLE "Working weight model, IP exposure weights - numerator";
 CLASS x  / PARAM=GLM;
 MODEL x = ; 
 OUTPUT OUT=num(KEEP=j ip_:) PREDPROBS=(i);
RUN;


*checking whether PARAM=ref gives same answer: yes, for this application; 
PROC LOGISTIC DATA = a OUT=nmod2;
 BY j;
 CLASS x  / PARAM=REF;
 MODEL x = ; 
 OUTPUT OUT=num2(KEEP=j ip_:) PREDPROBS=(i);
RUN;
 

 
PROC LOGISTIC DATA = a OUT=dmod;
 BY j;
 TITLE "Working weight model, IP exposure weights - denominator";
 CLASS v(REF="1") x  / PARAM=GLM;
 MODEL x = z v / LINK=GLOGIT; 
 OUTPUT OUT=den(KEEP=j ip_:) PREDPROBS=(i);
RUN;

PROC LOGISTIC DATA = a OUT=dmod2;
 BY j;
 CLASS v(REF="1") x  / PARAM=REF;
 MODEL x = z v / LINK=GLOGIT; 
 OUTPUT OUT=den2(KEEP=j ip_:) PREDPROBS=(i);
RUN;
/**/
OPTIONS MERGENOBY=nowarn;
DATA num; SET num; 
 ARRAY ip_[&dosecats];ARRAY iNUMp_[&dosecats];DO i = 1 TO &dosecats; iNUMp_[i] = ip_[i];END; DROP i ip:;
DATA den; SET den; 
 ARRAY ip_[&dosecats];ARRAY iDENp_[&dosecats];DO i = 1 TO &dosecats; iDENp_[i] = ip_[i];END; DROP i ip:;
DATA num2; SET num2; 
 ARRAY ip_[&dosecats];ARRAY iNUMp_[&dosecats];DO i = 1 TO &dosecats; iNUMp_[i] = ip_[i];END; DROP i ip:;
DATA den2; SET den2; 
 ARRAY ip_[&dosecats];ARRAY iDENp_[&dosecats];DO i = 1 TO &dosecats; iDENp_[i] = ip_[i];END; DROP i ip:;
OPTIONS mergenoby=NOWARN;
DATA b;
 MERGE a den num;
 ARRAY iDENp_[&dosecats];
 ARRAY iNUMp_[&dosecats];
 iptw_u = 1/iDENp_[x];
 iptw_s = iNUMp_[x]/iDENp_[x];
 DROP inum: iden:;
RUN;
DATA c;
 MERGE b den2 num2;
 ARRAY iDENp_[&dosecats];
 ARRAY iNUMp_[&dosecats];
 iptw_u2 = 1/iDENp_[x];
 iptw_s2 = iNUMp_[x]/iDENp_[x];
 DROP inum: iden:;
RUN;
OPTIONS mergenoby=WARN;


PROC LOGISTIC DATA = c OUT=cr;
 BY j;
 TITLE "Crude";
 CLASS x / PARAM=glm;
 MODEL y = x;
QUIT;
PROC LOGISTIC DATA = c OUT=p;
 BY j;
 TITLE "MSM, ipw 1";
 CLASS x / PARAM=glm;
 MODEL y = x;
 WEIGHT iptw_s;
QUIT;

PROC LOGISTIC DATA = c OUT=p2;
 BY j;
 TITLE "MSM, ipw 2";
 CLASS x / PARAM=glm;
 MODEL y = x;
 WEIGHT iptw_s2;
QUIT;
ODS LISTING;
