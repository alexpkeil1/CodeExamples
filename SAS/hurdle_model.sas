DATA A;
 CALL STREAMINIT(18817212);
 DO id = 1 to 2000;
  x = rand('bernoulli', 0.5);
  y = rand('normal', 5 + x); *normal outcome;
  y2 = EXP(rand('normal', log(5) + x, 0.5)); *lognormal outcome;
  ly2 = LOG(y2);
  y3 = EXP(rand('normal', -0.15 - 2*x, 0.0005))*RAND('bernoulli',1/(1+exp(-(-6 + x)))); *zero inflated log-normal outcome;
  y3gt0 = (y3>0);
  OUTPUT;
 END;
RUN;;
PROC MEANS DATA = a FW=5 MAXDEC=3;
 TITLE "sample data";
RUN;;


PROC OPTMODEL FD=CENTRAL PRESOLVER=3;
 SET OBS, NINTS={1..1}, Nbeta={1..1}, Nparms={1..2};
 NUMBER x{OBS, NBETA}, y3{OBS}, death{OBS}, cohort{OBS};
 READ DATA a INTO OBS=[_N_] y3;

 READ DATA a INTO OBS=[_N_] {j IN nbeta}  <x[_N_,j]>; 
 READ DATA a INTO OBS=[_N_] cohort;
 NUMBER pi = 3.1415926535897932384626433832795;
 
  *model components;
 VAR beta0  >= -400 <= 400; 
 VAR alpha0  >= -400 <= 400; 
 VAR betaX{Nbeta}  >= -400 <= 400;
 VAR alphaX{Nbeta}  >= -400 <= 400;
 VAR sig INIT 2 ;*>=1 <=200;
 
 IMPVAR eta{i IN OBS} = alpha0 + SUM{j IN Nbeta} x[i,j]*alphaX[j];
 IMPVAR logmu{i IN OBS} = beta0 + SUM{j IN Nbeta} x[i,j]*betaX[j];

 IMPVAR peq0{i IN OBS} = 1/(1+EXP(eta[i]));
 IMPVAR __yhat{i IN OBS} = EXP(logmu[i])*(1-peq0[i]);

 *MLE using log-likelihood;
 IMPVAR lli{i IN OBS}  = IF y3[i]>0
  THEN LOG(1-peq0[i]) - LOG(SQRT(2*pi)*sig) - ((y3[i]-EXP(logmu[i]))**2)/(2*sig**2) 
  ELSE LOG(peq0[i]);
  MAX ll = SUM{i IN OBS} lli[i];
  IMPVAR deviance = -2*LL;

 *storage matrices;
  SOLVE WITH nlp / MULTISTART MSNUMSTARTS=2; *works fine in this setting;
 PRINT ll sig beta0  betax alpha0  alphax  ;
 
QUIT;

PROC NLMIXED DATA = a ;
 TITLE 'ZI log-linear model';
  pi = 3.1415926535897932384626433832795;
  mu = EXP(b0 + b1*x);
  eta = a0 + a1*x;
  peq0 = 1/(1+exp(eta));
  expy=mu*(1-peq0);
  IF Y3>0 THEN 
   ll = LOG(1-peq0) - LOG(SQRT(2*pi)*sig) -((y3-mu)**2)/(2*sig**2); 
  ELSE 
   ll = LOG(peq0);
  MODEL y3 ~ general(LL);
 PREDICT expy OUT=bnl2(RENAME=(PRED=pynl) DROP=df tvalue probt alpha lower upper);
 PREDICT mu OUT=bnl2a(RENAME=(PRED=munl) DROP=df tvalue probt alpha lower upper);
 PREDICT eta OUT=bnl2b(RENAME=(PRED=etanl) DROP=df tvalue probt alpha lower upper);
RUN;

*comparison;
PROC GENMOD data = a ;
 TITLE 'ZI log-linear model (2 pieces), pt 1';
 WHERE Y3 > 0;
 MODEL y3 = x / LINK=log DIST=normal  MAXITER=300 SCORING=50 ;
 OUTPUT OUT = bgm2a p=pya;
RUN;
PROC GENMOD data = a DESC ;
 TITLE 'ZI log-linear model (2 pieces), pt 2';
 MODEL y3gt0 = x / LINK=logit DIST=binomial  MAXITER=300 SCORING=50 ;
OUTPUT OUT = bgm2b p=pyb;
RUN;;


*now with proc fmm (HURDLE MODEL);
PROC FMM DATA=a gconv=0;
 TITLE 'Finite mixture: log-linear hurdle model';
  MODEL y3 = x / DIST=normal LINK=log;
  MODEL y3 = / DIST=constant;
  PROBMODEL x / LINK=logit; *logit/glogit is default;
  ODS OUTPUT ParameterEstimates=hurdle;
  ODS OUTPUT MixingProbs=hurdle2;
RUN;
