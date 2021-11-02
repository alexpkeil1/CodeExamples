/**********************************************************************************************************************
* Author: Alex Keil
* Program: highly_adaptive_lasso.sas
* Date: 20180815
* Project: implement Highly adaptive LASSO (HAL) predictor in SAS
* Description: 1) Implements HAL in a simple curve fitting exercise from Benkeser and van der Laan (2016) IEEE ICDSAA
*              2) implements HAL within super learner for simple multivariate regression
*              3) Compares HAL with approximate HAL, where LASSO parameter is chosen by AICC, rather than cross validation
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
* TODO: speed optimization
**********************************************************************************************************************/
*clear the log window and the output window;
DM LOG 'clear;' CONTINUE; DM OUT 'clear;' CONTINUE; 
OPTIONS MERGENOBY = warn NODATE NONUMBER LINESIZE = 120  PAGESIZE=80 SKIP = 2 FORMDLIM = '-' MPRINT NOCENTER; TITLE1;TITLE2;
OPTIONS FORMCHAR = '|----|+|---+=|-/\<>*';

/*DM ODSRESULTS 'clear;' CONTINUE; *clear ODS generated datasets;*/
ODS LISTING GPATH="C:\Users\akeil\TEMP";


FILENAME slgh URL "https://cirl-unc.github.io/SuperLearnerMacro/super_learner_macro.sas";
%INCLUDE slgh;

* simulation illustration 1 from HAL paper, Benkeser and van der Laan 2016 IEEE ICDSAA;
* simple curve fitting exercise with comparison to random forest;
DATA sim(DROP=i);
 CALL STREAMINIT(1229);
  DO i = 1 to 500;
   x = rand('uniform')*8-4;
   M = rand('uniform')*8-4;
   *my = -2*(x<-3) + 2.55*(x>-2)-2*(x>0)+4*(x>2)-1*(x>3);
   my = 2*SIN(CONSTANT('pi')/2*ABS(x));
   y = my + RAND('gaussian');
   OUTPUT;
  END;
 RUN;

/*
Step 1:
Create the design matrix of indicator functions of all possible interactions, removing duplicate columns
*/
%MACRO __RENAMEDESIGN(VARS);
  /* rename all features
     - by necessity these are very short names "_#", which is needed to allow
       as many terms as will fit.
  */
  %GLOBAL __SLhalvars;
  %LET __slhalvars = ;
  %__SLnote(__RENAMEDESIGN);
  %LOCAL l_l ;
  %LET l_l = 1;
  RENAME 
  %DO %WHILE(%SCAN(&vars, &l_l)^=);
    %LET __slhalvars = &__slhalvars _&l_l;
    %SCAN(&vars, &l_l) =  _&l_l
    %LET l_l = %EVAL(1+&l_l);
  %END;;
%MEND __RENAMEDESIGN;

%MACRO GetHALDesign(y, vars, indata, outdata, drop=TRUE, intdepthlimit=0);
  /*
   steps to create design matrix for HAL (highly adaptive lasso):
    0. preprocessing
    1. create indicator variables for I[i,j,m] = I(X[i,m] <= X[j,m]) for i = 1..N, j=1..N, m=1..P
      - this step is open to some optimization by checking for binary/nominal variables
    2. create all possible interactions of those indicators: e.g. I[i,j,m+1] = I[i,j,1]*I[i,j,2]
    3. remove duplicate columns from the design matrix
   *REQUIRED: training data should have non-missing Y, testing/validation data should have missing Y
    - this is satisfied by super learner macro defaults

  to do:
   1. Optimize selection of unique columns in design matrix (currently involves 2 expensive transpositions)
  */
  %LOCAL haln halnt __hidx  __halmains __nhalmains __halbar __hallimit;
  /*
    0. preprocessing
  */
  DATA __haltm0001_;
    SET &indata(KEEP=&y &vars);
    %__RENAMEDESIGN(&vars);
  PROC SQL NOPRINT; SELECT COUNT(%SCAN(&__slhalvars, 1)) INTO :haln FROM __haltm0001_; QUIT;
  PROC SQL NOPRINT; SELECT COUNT(%SCAN(&__slhalvars, 1)) INTO :halnT FROM __haltm0001_(WHERE=(&Y>.z)); QUIT;
  %LET __hidx = 1;
  %LET __halbar = %SCAN(&__slhalvars, &__hidx) ;
  %DO %WHILE(%SCAN(&__slhalvars, &__hidx+1)^=);
    %LET __hidx = &__hidx+1;
    %LET __halbar = &__halbar | %SCAN(&__slhalvars, &__hidx);
  %END; 
  DATA fake;
    SET __haltm0001_ (OBS=1);
	WHERE &y>.z;
  * main terms only first;
  PROC GLMMOD DATA=fake OUTPARM=__halmains NOPRINT NAMELEN=200;
    MODEL &Y = &__slhalvars / NOINT ;
  PROC SQL NOPRINT;
    SELECT STRIP(EFFNAME) INTO :__halmains SEPARATED BY " " 
      FROM __halmains;
  DATA __halmains;
   SET __halmains END=EOF;
   IF eof THEN CALL SYMPUT("__nhalmains", PUT(_COLNUM_, 9.0));
  RUN;  
  * limiting the depth of interaction terms to speed up processing;
  %IF &intdepthlimit ^=0 %THEN %LET __hallimit =  &intdepthlimit;
  %ELSE %IF &intdepthlimit =0 %THEN %LET __hallimit = &__NHALMAINS;
  %IF %EVAL(&__hallimit < &__NHALMAINS) %THEN %PUT NOTE: HAL - depth of maximum product term limited by user to &__hallimit - design matrix may not be completely saturated;;
  * leverage glmmod for creating all possible interaction terms;
  PROC GLMMOD DATA=fake OUTPARM=__halterms  NOPRINT NAMELEN=200;
    MODEL &Y = &__halbar @&__hallimit / NOINT ;
  DATA __halterms;
   SET __halterms END=EOF;
   IF eof THEN CALL SYMPUT("__nhalterms", PUT(_COLNUM_, 9.0));
  RUN;
  *storing terms in file as code to get around macro length limit;
  FILENAME terms TEMP;
  DATA _NULL_;
    FILE terms;
    IF _N_=1 THEN DO;
	  PUT "*begin creation of full design matrix;";
      PUT "DO n = 1 TO %TRIM(&halnt);";
	END;
    SET __halterms END=eof;
	effname = STRIP(effname);
    t = 1;
	PUT '  __Hi[' _N_ ',n] = 1' @@;
	DO WHILE(SCAN(effname, t, '*')^='');
	  nm = STRIP(SCAN(effname, t, '*'));
      PUT '*(' nm ">= _" nm +(-1)'_[n])' @@;
      t=t+1;
	END;
	PUT ';';
    IF EOF THEN DO;
      PUT 'END;';
	  PUT '*end creation of full design matrix;';
	END;
  RUN;
  /*1. create indicator variables for I[i,j,m] = I(X[i,m] <= X[j,m]) for i = 1..N, j=1..N, m=1..P*/
  DATA __halmat;
	*this step should be done only in the training data!;
    SET __haltm0001_;
    WHERE &y > .z; 
    %LET __hidx = 1;
    %DO %WHILE(%SCAN(&__slhalvars, &__hidx, ' ')^=);
      %LET __x = %SCAN(&__halmains, &__hidx, ' ');
      ARRAY _&__x._[&halnt] ;
      %LET __hidx = %EVAL(&__hidx+1);
	  RETAIN _&__x._:;
	  KEEP _&__x._:;
    %END;
	*outer loop over variables, inner loop over sample;
    %LET __hidx = 1;
    %DO %WHILE(%SCAN(&__slhalvars, &__hidx, ' ')^=);
      DO i = 1 TO &halnt;
        IF i = _N_ THEN _%SCAN(&__slhalvars, &__hidx, ' ')_[i] = %SCAN(&__slhalvars, &__hidx, ' ');
      END;
      %LET __hidx = %EVAL(&__hidx+1);
    %END;
    IF _n_ = &halnt THEN OUTPUT;  
  RUN;
  /*2. create all possible interactions of those indicators: e.g. I[i,j,m+1] = I[i,j,1]*I[i,j,2]*/
  /*
  note: some speedup could be obtained by keeping (as main term) unmodified, all binary variables
  */
  DATA __designdupa(KEEP = __hi:);
  /* main computational slowdown */
    ARRAY __Hi[&__nhalterms, &halnt] ;
    LENGTH __Hi: 3;
    SET __haltm0001_;
    IF _n_=1 THEN SET __halmat;
    %LET __hidx = 1;
    %DO %WHILE(%SCAN(&__slhalvars, &__hidx, ' ')^=);
      %LET __x = %SCAN(&__halmains, &__hidx, ' ');
      ARRAY _&__x._[&halnt] ;
      %LET __hidx = %EVAL(&__hidx+1);
    %END;
    %INCLUDE terms;
  RUN;
  /*
  3. remove duplicate columns from the design matrix;
  */
  PROC TRANSPOSE DATA = __designdupa OUT=__designdup(DROP=_NAME_); 
  PROC SQL NOPRINT;
    DROP TABLE __designdupa;
    CREATE TABLE &outdata AS SELECT DISTINCT * FROM __designdup;
    %IF &DROP=TRUE %THEN DROP TABLE __designdupa, __halmat, __halmains, __halterms, __haltm0001_;;
  QUIT;
  PROC TRANSPOSE DATA = &outdata OUT=&outdata(DROP=_NAME_) PREFIX=__hal_;
  RUN; QUIT; RUN;
%MEND GetHALDesign;

OPTIONS noMPRINT;
%GetHALDesign(y=y, vars=x, indata=sim, outdata=design, drop=FALSE);
OPTIONS NOMPRINT;

OPTIONS MERGENOBY=NOWARN;
DATA an;
 MERGE sim(KEEP=my Y x) design;
RUN;
OPTIONS MERGENOBY=WARN;

/* 
Step 2:
Run LASSO with cross validated selection of shrinkage parameter
*/
TITLE;
* LASSO with cross validation selection of shrinkage parameter (continuous outcome only);
PROC GLMSELECT DATA = an SEED=12039;
  MODEL Y = __hal_: / SELECTION=LASSO(CHOOSE=CVEX STOP=NONE)  CVDETAILS=all CVMETHOD=RANDOM(10);
  OUTPUT OUT =outhal(KEEP=p_hal my x y) PRED=p_hal;
  ODS OUTPUT parameterestimates=lassoparms;
RUN;

* HAL approximation using AICC, method 2;
PROC HPGENSELECT DATA=AN TECH=QUANEW NOSTDERR LASSORHO=0.8 LASSOSTEPS=100 NORMALIZE=YES;
   MODEL Y = __hal_: / DIST=GAUSSIAN LINK=ID;
  SELECTION METHOD=lasso(CHOOSE=AICC STOP=AICC);
  OUTPUT OUT =outhal3(KEEP=p_Ahal) PRED=p_Ahal;
RUN;;

* HAL approximation using reduced form of (non-exhaustive search) CV, method 1;
PROC GLMSELECT DATA = an SEED=12039;
  MODEL Y = __hal_: / SELECTION=LASSO(CHOOSE=CVEX STEPS=100)  ;
  OUTPUT OUT =outhal2(KEEP=p_ahala my x y) PRED=p_ahala;
  ODS OUTPUT parameterestimates=lassoparms;
RUN;


PROC SQL; 
 TITLE "l1-norm bound (truth = 16)";
 SELECT SUM(ABS(estimate)) AS M FROM work.lassoparms;
QUIT;
TITLE;

PROC HPFOREST DATA = OUTHAL SCOREPROLE=OOB;
  ID _all_;
  TARGET y / LEVEL = interval;
  INPUT  x / LEVEL = interval;
  SCORE OUT = OUTPRED(KEEP=p_hal p_y my x y RENAME=(p_y=p_RF));
run;
DATA outpred;
 MERGE outpred outhal2 outhal3;

PROC SORT DATA = outpred;
 BY x;
ODS GRAPHICS / RESET=ALL IMAGENAME = "hal" IMAGEFMT=PNG;
PROC SGPLOT DATA = outpred;
 TITLE 'Simulation 1 Benkeser and van der Laan (2016)';
 LABEL y="Y" my="Truth" p_hal="HAL" p_ahala="Approximate HAL (glmselect)" p_Ahal="Approximate HAL (hpgenselect)" p_rf = "Random Forest";
 SCATTER Y=Y X=X / MARKERATTRS=(COLOR='BLUE' SYMBOL=CircleFilled SIZE=4);
 SERIES y=my X=X / LINEATTRS=(COLOR='BLACK' PATTERN = DASH THICKNESS=2);
 STEP y=p_rf X=X / LINEATTRS=(COLOR='GREEN'  THICKNESS=1);
 STEP y=p_hal X=X / LINEATTRS=(COLOR='RED'  THICKNESS=2);
 STEP y=p_ahala X=X / LINEATTRS=(COLOR='RED'  PATTERN = DOT THICKNESS=2);
 STEP y=p_Ahal X=X / LINEATTRS=(COLOR='RED'  PATTERN=DASH THICKNESS=2);
 YAXIS LABEL="Y";
RUN;

/* 
Step 3:
use in super learner
*/

TITLE "Super learner fit";
DATA train valid ;
  LENGTH id x l 3;
  CALL STREAMINIT(2211887);
  DO id = 1 TO 1100;
    u = RAND("uniform")*0.1 + 0.4;  
    l = RAND("bernoulli", 1/(1+exp(-1 + 2*u)));
    c = RAND("normal", u, 1);
    c2 = RAND("normal", u, .3);
    x = RAND("bernoulli", 1/(1+exp(-1.5 + 2*l + c + c2)));
    y = RAND("NORMAL", 2*u + x, 0.5);
    KEEP x l c c2 y;
    IF id <= 100 THEN OUTPUT train;
    ELSE OUTPUT valid;
  END;
RUN;


%SuperLearner(Y=y, X=x l c c2,
             indata=train, 
             preddata= valid,
             library=bagging rf boost hal ahal ahalb, /* ahal and ahalb are approximations of HAL using a reduced form of the CV lasso to select the lasso parameter */
             folds=10, 
             dist=GAUSSIAN, 
             method=NNLS);


DATA mse(KEEP=__train squarederror:);
  SET sl_out;
  squarederror_sl = (y - p_SL_full)**2;
  squarederror_bagging = (y - p_bagging_full)**2;
  squarederror_rf = (y - p_rf_full)**2;
  squarederror_boost = (y - p_boost_full)**2;
  squarederror_hal = (y - p_hal_full)**2;
  squarederror_ahal = (y - p_ahal_full)**2;
  squarederror_ahalb = (y - p_Ahal_full)**2;
PROC MEANS DATA = mse FW=5 MEAN;
  TITLE 'Mean squared error of predictions in training/validation data';
  VAR squarederror:;
  CLASS __train;
RUN;
