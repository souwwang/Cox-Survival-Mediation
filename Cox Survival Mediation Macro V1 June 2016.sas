*************************************************************************************************************************************************************************** ;
*	ESTIMATING MEDIATION EFFECTS FOR THE COX PROPORTIONAL HAZARDS MODEL                                                                                                     ;
*   Version 1.1, Uploaded 06/15/2016 												                                                                                        ;
*																					               										                                    ;
*	Code written by Wei Wang (2016)													               										                                    ;
* 																					                										                                ;
*	Reference:																		               										                                    ;
* 	Wang W, and Albert JM (2016). Causal Mediation Analysis for the Cox Proportional Hazards Model with a Smooth Baseline Hazard Estimator.                                 ;
* 	Journal of Royal Statistical Society, Series C. Submitted. 														                                                        ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
* 																					                										                                ;
*	PROGRAM DESCRIPTION:															                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	This program implement a mediation formula approach to esimate the total effect, natural indirect effect, and natural direct effect in terms of                         ;
*   survival probability, hazard function as well as restricted  mean survival time for the survival Cox proportional hazards model within the standard                     ;
*   two-stage mediation framework. Simple parametric models (fractional polynomials and restricted cubic splines) are utilized to approximate the                           ;
*   baseline log cumulative hazard function used for the estimation of the related causal quantities. This SAS macro can be applied for one binary exposure,                ;
*   one normally distribued mediator and one survival outcome under the proportional hazards assumption.                                                                    ;
* 																					                                                                                        ;
*	DATASET REQUIREMENTS:															                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	The macro requires a SAS data set with the following variables:                                                                                                         ;
*   (1) Survival outcome variable (SURVtime)                                                                                                                                ;
*   (2) Survival event indicator (SURVevent)                                                                                                                                ;
*   (3) Normally distributed continuous mediator (MED)                                                                                                                      ;
*   (4) List of covariates for the survival outcome and the mediator excluding exposure (Cov)                                                                               ;
*   (5) Binary exposure (EXP).                                                                                                                                              ;
*                                                                                                                                                                           ;
*   Of note, the covariates for the survival outcome and the mediator can be overlapped, or it can be possible that two sets of covariates are identical.                   ;
* 																					                										                                ;
*   Format of data set exmaple:                                                                                                                                             ;
* 																					                										                                ;
*	ID	SURVTIME    SURVevent MED      Cov1     Cov2     Cov3        EXP                                                                                       			    ;
*	1	3.2	        1         22.5     13.2	    1		 21		     0			                                                                                            ;
*	2	10.5	    0         10.2     16.9		0		 23		     1				                                                                                        ;
*	3	7.1	        1         13.5     12.0		0		 23		 	 0					                                                                                    ;
*	4	8.2	        0         21.3     23.5		1		 12		 	 1					                                                                                    ;
*	5	4.5	        0         21.7     21.0		1		 13		 	 1				                                                                                        ;
*	6	12.0	    1         23.3     12.7		1		 45		 	 1						                                                                                ;
*																					                                                                                        ;
*	MODEL SPECIFICATION																                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	X: binary exposure                                                                                                                                                      ;
*   M: normally distributed mediator modelled with linear regression model                                                                                                  ;
*   W: baseline covariate                                                                                        			                                                ;
*   T: time-to-event outcome                                                                                       			                                                ;
*																					                                                                                        ;
*	Mediator, M = ALPHA0 + ALPHA1 * X + ALPHA * W + EPSILON, WHERE EPSILON ~ N(0, SIGMA^2)    	   				                                                            ;
*																					                                                                                        ;
*   Proportional hazards model for time-to-event outcome defined through the hazard function, h(t) = h0(t)exp(beta1 * X + beta2 * M + beta * W)                             ;
*																					                                                                                        ;
*	CAUSAL MEDIATION EFFECTS DEFINITION                                         	                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	POTENTIAL OUTCOME FRAMEWORK WAS USED TO DEFINE THE TOTAL EXPOSURE EFFECTS, NATURAL INDIRECT EFFECT, NATURAL DIRECT EFFECT, AND MEDIATION PROPORTION                   	;
*	IN TERMS OF SURVIVAL PROBABILITY, HAZARD FUNCTION, AND RESTRICTED MEAN SURVIVAL TIME WITHIN TWO-STAGE MEDIATION FRAMEWORK                                               ;
*																					                                                                                        ;
*	Causal Quantities Definition on Restricted Mean Survival Time (RMST)    		                                                                                        ;
*																	        		                                                                                        ;
*	Total Exposure Effect: TE(t*) = mu(t*: 1, M(1)) - mu(t*: 0, M(0))			                                                                                            ;
*																					                                                                                        ;
*	Natural Indirect Effect: IE(t*: x) = mu(t*: x, M(1)) - mu(t*: x, M(0))							                                                                        ;
*																					                                                                                        ;
*	Natural Direct Effect: DE(t*: 1-x) = mu(t*: 1, M(1-x)) - mu(t*: 0, M(1-x))    					                                                                        ;
*																					                                                                                        ;
*	Total Exposure Effect Decomposition: TE(t*) = IE(t*: x) + DE(t*: 1-x)			                                                                                        ;
*																					                                                                                        ;
*	Mediation proportion: P(t*: x) = IE(t*: x) / TE(t*)								                                                                                        ;
*																					                                                                                        ;
*																					                                                                                        ;
*	Causal Quantities Definition on Survival Probability 		     				                                                                                        ;
*																	        		                                                                                        ;
*	Total Exposure Effect: TE(t*) = S(t*: 1, M(1)) - S(t*: 0, M(0))			                                                                                                ;
*																					                                                                                        ;
*	Natural Indirect Effect: IE(t*: x) = S(t*: x, M(1)) - S(t*: x, M(0))							                                                                        ;
*																					                                                                                        ;
*	Natural Direct Effect: DE(t*: 1-x) = S(t*: 1, M(1-x)) - S(t*: 0, M(1-x))    					                                                                        ;
*																					                                                                                        ;
*	Total Exposure Effect Decomposition: TE(t*) = IE(t*: x) + DE(t*: 1-x)			                                                                                        ;
*																					                                                                                        ;
*	Mediation proportion: P(t*: x) = IE(t*: x) / TE(t*)								                                                                                        ;
*																					                                                                                        ;
*																					                                                                                        ;
*	Causal Quantities Definition on hazard function 	     	     				                                                                                        ;
*																	        		                                                                                        ;
*	Total Exposure Effect: TE(t*) = h(t*: 1, M(1)) / h(t*: 0, M(0))			                                                                                                ;
*																					                                                                                        ;
*	Natural Indirect Effect: IE(t*: x) = h(t*: x, M(1)) / h(t*: x, M(0))							                                                                        ;
*																					                                                                                        ;
*	Natural Direct Effect: DE(t*: 1-x) = h(t*: 1, M(1-x)) / h(t*: 0, M(1-x))    					                                                                        ;
*																					                                                                                        ;
*	Total Exposure Effect Decomposition: TE(t*) = IE(t*: x) * DE(t*: 1-x)			                                                                                        ;
*																					                                                                                        ;
*	Mediation proportion: P(t*: x) = log (IE(t*: x)) / log (TE(t*))					                                                                                        ;
*																					                                                                                        ;
*	Of note, for the defined causal mediation effects on each of the RMST, survival function and hazard function scales, we need to choose a time point t*                  ;
*	at which each function is evaluated. The t* (macro variabel t_assess see below) can be chosen as any value that are clincially motivated before the last 		        ;
*	observed event time in order to obtain 'accurate' baseline parametric distribution approximation from the FP and RCS models for the mediaiton effects estimation.       ;
*																					                                                                                        ;
*	MACRO VARIABLES:																                                                                                        ;
*   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- ;
*	DATA: 		    SAS data set to contain all survival mediation analysis variables                                                                              	        ;
*																					                                                                                        ;
*	SURVtime:		time-to-event variable	                                                                                                                                ;
*																					                                                                                        ;
*	SURVevent:		event indicator, event = 1, censor = 0                             		                                                                                ;
*																					                                                                                        ;
*	MED:		    normally distributed continuous mediator            	                                                                                                ;
*																					                                                                                        ;
*	Cov:		    List of covariates for the survival outcome and the mediator excluding exposure                                                                         ;
*																					                                                                                        ;
*	EXP: 		    binary exposure                                      	                                                                                                ;
*																					                                                                                        ;
*	Nboot:          number of bootstrap samples to calculate the 95% CI, default Nboot = 200 (suggested)														            ;
*																					                                                                                        ;
*	T_ACCESS:		t*, time point at which the causal quantities are defined, default (999999) is the last observed event time, the macro users may replace 999999         ;
*					with other arbitrarily designated t* (t* > 0, strongly recommended to be less than the last observed event time, see discussion section of published    ;
*					reference by Wang and Albert JRSS Series C)								                                                                                ;
*																					                                                                                        ;
*	REFERENCE:		can be 0, 1, 2 or 3. The estimated causal quantities are determined by the baseline covariate (predictors for the survival outcome,                     ;
*					heterogeneous causal quantities estimation for different baseline covariates) and the estimated 														;
*					and presented average causal quantities are calculated as the average of causal quantities over the empirical covariates distribution from              ;
*					the designated reference population (see reference for more details). The reference population can be exposed group, unexposed group, total             ;
*					population or any designated population with arbitrarily determined baseline covariate distribution, default vaule for reference is 1.                  ;
*					0: unexposed group as the reference population, and average causal quantities are calculated over the empirical covariates distribution                 ;
*					   from the unexposed group															                                                                    ;
*					1: exposed group as the reference population, and average causal quantities are calculated over the empirical covariates distribution from              ;
*					   from the exposed group															                                                                    ;
*					2: total population as the reference group, and average causal quantities are calculated over the empirical covariates distribution                     ;
*					   from the total population															                                                                ;
*					3: users provided population with pre-specified baseline covariate distribution as the reference group. If reference is 3, the Macro user should        ;
*					   provide a temporary dataset named as REF to show the baseline covariate distribution in which the users want to estimate the causal quantities.      ;
*					   The dataset REF should follow the following format:															                                        ;
*																					                                                                                        ;
*	                          Cov1         Cov2         Cov3         NSUB                                                                                      			    ;
*                             0		       0		    0			 10                                                                                                     ;
*                             0		       0		    1			 10                                                                                                     ;
*                             0		       1		    0			 20                                                                                                     ;
*                             0		       1		    1			 20                                                                                                     ;
*                             1		       0		    0			 10                                                                                                     ;
*                             1		       0		    1			 10                                                                                                     ;
*                             1		       1		    0			 10                                                                                                     ;
*                             1		       1		    1			 10                                                                                                     ;
*																					                                                                                        ;
*						If the dataset REF is defined as above, the reference population is defined as following, in this reference population                              ;
*						number of subjects with Cov1 = 0, Cov2 = 0 and Cov3 = 0 is 10, number of subjects with Cov1 = 0, Cov2 = 0 and                                       ;
*						Cov3 = 1 is 10,	number of subjects with Cov1 = 0, Cov2 = 1 and Cov3 = 0 is 20 etc.                                                                  ;
*																					                                                                                        ;
*					Of note, the survival model can include continuous covariates, however, since the required avarage causal quantities are                                ;
*					calculated over the empirical baseline covariate distribution from the reference population, which means causal quantities will be calculated           ;
*					for each distinct baseline covariate combination, so when continuous covariates are included, the computation time will be pretty long.                 ;
*					With above REF data set, 8 distinct binary covariate combinations are used. If we want to calculate average causal quantities                           ;
*					using reference population with 100 subjects including one continuous covariate, 100 distinct covariate combinations will be used. The running time     ;
*					will be around 12.5 folds longer (100/8). So for dataset with continuous covariates, the author suggests that the users provide arbitrarily defined     ;
*					reference group and do not use original population as reference to save the computation time. For example, if we have two survival model covariates,    ;
*					one is continuous age ~ N(55, 10), the other covariate is binary, the authors suggest using REFERENCE = 3 option, and provide dataset REF as following, ;
*																					                                                                                        ;
*	                          Cov1         Cov2           NSUB                                                                                      			            ;
*                             55		   0		      1                                                                                                                 ;
*                             55		   1		      1                                                                                                                 ;
*																					                                                                                        ;
*					The reported causal quantities should be clearly specified over a reference population with age (Cov1) 55 years old and half subjects have Cov2         ;
*					as 0 and half subjects have Cov2 as 1.																                                                    ;
*																					                                                                                        ;
*	NLMOPTS:		User specified PROC NLMIXED procedure option, if missing then default option for NLMIXED                                                                ;
*                   procedure is used for model fitting														                                                                ;
*																					                                                                                        ;
*	PRINTALL:		Outpout data set print out or not	                                                                                                                    ;
*					F: output data set is not printed out								                                                                                    ;
*					T: output data set is printed out, default												                                                                ;
*																					                                                                                        ;
*	OUT:	        Output dataset                                                                                                                                          ;
*																					                                                                                        ;
*																					                                                                                        ;
*	EXAMPLE CODE:																	                                                                                        ;
*																					                                                                                        ;
*   %include 'SURVPHMED.sas'														                                                                                        ;
*																					                                                                                        ;
*   %SURVPHMED (DATA             = JHSdata,                                                                                                                                 ;
*               SURVtime         = months,                                                                                                                                  ;
*               SURVevent        = stroke,                                                                                                                                  ;
*               MED              = sbp,                                                                                                                                     ;
*               Cov              = bmic agec male,                                                                                                                          ;
*               EXP              = idealhealthsmk,                                                                                                                          ;
*               Nboot            = 200,                                                                                                                                     ;
*				T_ACCESS         = 102.28,                                                                                                                                  ;
*				REFERENCE        = 1,                                                                                                                                       ;
*               Nlmopts          = ,                                                                                                                                        ;
*               Printall         = T,                                                                                                                                       ;
*               Out              = out01                                                                                                                                    ;
*              )                                                                                                                                                            ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
* 																					                                                                                        ;
*************************************************************************************************************************************************************************** ;
***********************************************************************************Final Macro******************************************************************************;
*************************************************************************************************************************************************************************** ;

dm  'log;clear;out;clear;';

OPTIONS ls=130  ps=57 NOCENTER DATE PAGENO=1; 

** Define SubMacro;

** RCS is sub Macro to fit the survival outcome with the restricted cubic splines with 1 or 2 or 3 interior knots to approximate the baseline log cumulative hazard function. 
Parameter estimates and fit criteria (AIC and BIC) are provided in the output data set;

%macro rcs();

  proc nlmixed data=_dat11 cov;

    parms / bydata data=SurvInit1;
    linpredSURV  = &linpredSURV; 
    linpredcrs = gam10 + gam11 * lnt + gam12 * nu11;
    linpred = linpredsurv + linpredcrs;
		  
    G_t    = -exp(linpred);
    g      = G_t - lnt  + linpred + log(gam11 + gam12 * dnu11);
    loglikeSURV = (_event=1)*g + (_event=0)*G_t;
    model lnt ~ general(loglikeSURV);

    ods output ParameterEstimates=temp15;
    ods output FitStatistics = temp16; 
    ods output ConvergenceStatus = temp17 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  %if %sysfunc(exist(temp16)) %then %do;

    data temp18;
      set temp16 temp17 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp18;
      by replicate;
    run;

    proc transpose data = temp18 out = temp19 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data gg45;
      set temp19;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

    data ff45;
      set temp15;
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

    proc datasets lib = work;
      delete temp15 temp16 temp17 temp18 temp19;
    run;

  %end;

  %else %do;

    data gg45;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

    data ff45;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 1;
      ID = 45;
    run;

  %end;

  proc nlmixed data=_dat11 cov;

    parms / bydata data=SurvInit2;
    linpredSURV  = &linpredSURV; 
    linpredcrs = gam20 + gam21 * lnt + gam22 * nu21 + gam23 * nu22;
    linpred = linpredsurv + linpredcrs;
		  
    G_t    = -exp(linpred);
    g      = G_t - lnt  + linpred + log(gam21 + gam22 * dnu21 + gam23 * dnu22);
    loglikeSURV = (_event=1)*g + (_event=0)*G_t;
    model lnt ~ general(loglikeSURV);

    ods output ParameterEstimates=temp15;
    ods output FitStatistics = temp16; 
    ods output ConvergenceStatus = temp17 (keep = replicate status rename = (status = value));
    by replicate; 
  run;

  %if %sysfunc(exist(temp16)) %then %do;

    data temp18;
      set temp16 temp17 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp18;
      by replicate;
    run;

    proc transpose data = temp18 out = temp19 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data gg46;
      set temp19;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

    data ff46;
      set temp15;
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

    proc datasets lib = work;
      delete temp15 temp16 temp17 temp18 temp19;
    run;

  %end;

  %else %do;

    data gg46;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

    data ff46;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 2;
      ID = 46;
    run;

  %end;

  proc nlmixed data=_dat11 cov;

    parms / bydata data=SurvInit3;
    linpredSURV  = &linpredSURV; 
    linpredcrs = gam30 + gam31 * lnt + gam32 * nu31 + gam33 * nu32 + gam34 * nu33;
    linpred = linpredsurv + linpredcrs;
		  
    G_t    = -exp(linpred);
    g      = G_t - lnt  + linpred + log(gam31 + gam32 * dnu31 + gam33 * dnu32 + gam34 * dnu33);
    loglikeSURV = (_event=1)*g + (_event=0)*G_t;
    model lnt ~ general(loglikeSURV);

    ods output ParameterEstimates=temp15;
    ods output FitStatistics = temp16; 
    ods output ConvergenceStatus = temp17 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  %if %sysfunc(exist(temp16)) %then %do;

    data temp18;
      set temp16 temp17 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp18;
      by replicate;
    run;

    proc transpose data = temp18 out = temp19 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data gg47;
      set temp19;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

    data ff47;
      set temp15;
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

    proc datasets lib = work;
      delete temp15 temp16 temp17 temp18 temp19;
    run;

  %end;

  %else %do;

    data gg47;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

    data ff47;
      format converge $3.;
      converge = 'No';
      method = 'RCS';
      nknot = 3;
      ID = 47;
    run;

  %end;

%mend;

** SINGLE is sub Macro to fit the survival outcome with the fractional polynomials with degree 1 to approximate the baseline log cumulative hazard function. 
Parameter estimates and fit criteria (AIC and BIC) are provided in the output data set;

%macro single(in1, in2, out1, out3, out2, n, id);

/*  dm  'log;clear;out;clear;';*/

  data temp21;
    set &in2;
  run;

  proc glm data = temp21;
    model logbasecumhaz = t&n;
    by replicate;
    ods output ParameterEstimates  = temp211;
  run;

  data temp212;
    set temp211;
    if parameter = 'Intercept' then parameter1 = 'lamd0';
    if parameter = "t&n" then parameter1 = 'lamd1';
	drop parameter;
	rename parameter1 = parameter;
    keep replicate parameter1 estimate;
  run;

  data temp213;
    set SURVests2 (keep = replicate model parameter estimate effect) temp212 (in = a);
    if a then model = 'FP';
  run;

  proc sort data = temp213;
    by replicate;
  run;

  proc nlmixed data=&in1 cov;
    parms /bydata data=temp213;

    linpredSURV  = &linpredSURV; 
    linpredcrs = lamd0 + lamd1 * t&n;
    linpred = linpredsurv + linpredcrs;
		  
    G_t    = -exp(linpred);
    g      = G_t + linpred + log(lamd1 * dt&n);
    loglikeSURV = (_event=1)*g + (_event=0)*G_t;
    model lnt ~ general(loglikeSURV);

    ods output ParameterEstimates=temp25;
    ods output FitStatistics = temp26; 
    ods output ConvergenceStatus = temp27 (keep = replicate status rename = (status = value));
	by replicate;
  run;

  %if %sysfunc(exist(temp26)) %then %do;

    data temp28;
      set temp26 temp27 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp28;
      by replicate;
    run;

    proc transpose data = temp28 out = temp29 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data &out2;
      set temp29;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

    data &out1;
      set temp25;
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

    proc datasets lib = work;
      delete temp25 temp26 temp27 temp28 temp29;
    run;

  %end;

  %else %do;

    data &out2;
      format converge $3.;
      converge = 'No';
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

    data &out1;
      method = 'FP';
      P1 = &n;
      ID = &ID;
    run;

  %end;

  proc datasets lib = work;
    delete temp21 temp211 temp212 temp213;
  run;

%mend;

** DOUBLE is sub Macro to fit the survival outcome with the fractional polynomials with degree 2 to approximate the baseline log cumulative hazard function. 
Parameter estimates and fit criteria (AIC and BIC) are provided in the output data set;

%macro double(in1, in2, out1, out3, out2, m, n, id);

/*  dm  'log;clear;out;clear;';*/

  data temp31;
    set &in2;
  run;

  proc glm data = temp31;
    model logbasecumhaz = t&m t&n;
    by replicate;
    ods output ParameterEstimates  = temp311;
  run;

  data temp312;
    set temp311;
    if parameter = 'Intercept' then parameter1 = 'lamd0';
    if parameter = "t&m" then parameter1 = 'lamd1';
    if parameter = "t&n" then parameter1 = 'lamd2';
	drop parameter;
	rename parameter1 = parameter;
    keep replicate parameter1 estimate;
  run;

  data temp313;
    set SURVests2 (keep = replicate model parameter estimate effect) temp312 (in = a);
    if a then model = 'FP';
  run;

  proc sort data = temp313;
    by replicate;
  run;

  proc nlmixed data=&in1 cov;

    parms /bydata data=temp313;

    linpredSURV  = &linpredSURV; 
    linpredcrs = lamd0 + lamd1 * t&m + lamd2 * t&n;
    linpred = linpredsurv + linpredcrs;
		  
    G_t    = -exp(linpred);
    g      = G_t + linpred + log(lamd1 * dt&m + lamd2 * dt&n);
    loglikeSURV = (_event=1)*g + (_event=0)*G_t;
    model lnt ~ general(loglikeSURV);

    ods output ParameterEstimates=temp35;
    ods output FitStatistics = temp36; 
    ods output ConvergenceStatus = temp37 (keep = replicate status rename = (status = value));
    by replicate;
  run;

  %if %sysfunc(exist(temp36)) %then %do;

    data temp38;
      set temp36 temp37 (in = a);
      if a then Descr = 'Con';
    run;

    proc sort data = temp38;
      by replicate;
    run;

    proc transpose data = temp38 out = temp39 (rename = (AIC__smaller_is_better_ = AIC AICC__smaller_is_better_ = AICC BIC__smaller_is_better_ = BIC));
      id descr;
      var value;
      by replicate;
    run;

    data &out2;
      set temp39;
      drop _name_;
      if con = 0.0 then converge = 'Yes';
      else converge = 'No';
      drop con;
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

    data &out1;
      set temp35;
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

    proc datasets lib = work;
      delete temp35 temp36 temp37 temp38 temp39;
    run;

  %end;

  %else %do;

    data &out2;
      format converge $3.;
      converge = 'No';
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

    data &out1;
      method = 'FP';
      P1 = &m;
      P2 = &n;
      ID = &ID;
    run;

  %end;

  proc datasets lib = work;
    delete temp31 temp311 temp312 temp313;
  run;

%mend;

** MED is sub Macro to calculate the mean potential outcomes used for the causal quantities definition under the potential outcome framework. In the MED sub Macro, it also called 16 
sub sub Macros (INTE11, INTE211, INTE212, INTE213, INTE214, INTE311, INTE312, INTE313, INTE314, INTE315, INTE411, INTE412, INTE413, INTE414, INTE415 and INTE416) to calculate the 
mean potential outcomes for FP baseline hazard function (1 submcro), RCS with degree of freedom 2 (4 submacros), RCS with degree of freedom 3 (5 submacros) and RCS with degree of 
freedom 4 (6 submacros) baseline hazard function respectively.;

/* FP Potential Outcome Mean Calculation */

%macro inte11(indata, outdata, x1, x2);

  data temp62;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp62;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp62;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    if &id >= 1 & &id <= 8 then do;
      lamd0 = baseline[1, 1];
      lamd1 = baseline[1, 2];
    end;

    if &id > 8 & &id <= 44 then do;
      lamd0 = baseline[1, 1];
      lamd1 = baseline[1, 2];
      lamd2 = baseline[1, 3];
    end;

/*print sigma alpha beta alpha2 lamd0 lamd1 lamd2 cov survx medx;*/

    out00=j(obs,1,.);
    out10=j(obs,1,.);
    out30=j(obs,1,.);

/* S(tmax) potential outcome */

    start fun00(m) global(w, beta, alpha, alpha2, sigma, sur, me, lamd0, lamd1, lamd2, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32); 
      s = &rmst;
      n = min(((&sx) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

/* f(tmax) potential outcome */

    start fun10(m) global(w, beta, alpha, alpha2, sigma, sur, me, lamd0, lamd1, lamd2, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32); 
      s = &rmst;
      n = min(((&sx) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) * exp(n) * (&dsx)/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

/* Restricted Mean Survival Time(tmax) potential outcome */

    start fun30(m) global(w, beta, alpha, alpha2, sigma, sur, me, lamd0, lamd1, lamd2, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, s); 
      n = min(((&sx) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out30(x) global (w, beta, alpha, alpha2, sigma, sur, me, lamd0, lamd1, lamd2, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun30", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a = { 0  &rmst };
      b = { .M  .P }; 
      call quad(z00, "fun00", b);
      out00[sub] = z00;
      call quad(z10, "fun10", b);
      out10[sub] = z10;
      call quad(z30, "out30", a);
      out30[sub] = z30;
    end;

    postprobs=subjid||out00||out10||out30;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 

    quit;

%mend;

/* RCS Mean Potential Outcome Calculation */

/* RCS DF = 2*/

/* Assessed t* <= smin */

%macro inte211(indata, outdata, x1, x2);

  data temp72;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp72;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp72;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outcome13=j(obs,1,.);
    outcome51=j(obs,1,.);
    outcome52=j(obs,1,.);

    start fun51(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun52(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx1) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun13", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &rmst };
      b = { .M  .P }; 

      call quad(z51,"fun51",b); 
      outcome51[sub] = z51;

      call quad(z52,"fun52",b); 
      outcome52[sub] = z52;

      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
    end;

    outcome1 = outcome51;
    outcome2 = outcome52;
    outcome3 = outcome13;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed smin < t* <= s11 */

%macro inte212(indata, outdata, x1, x2);

  data temp72;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp72;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp72;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome51=j(obs,1,.);
    outcome52=j(obs,1,.);

    start fun51(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun52(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx2) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun23", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin };
      a2 = { &smin &rmst};
      b = { .M  .P }; 

      call quad(z51,"fun51",b); 
      outcome51[sub] = z51;

      call quad(z52,"fun52",b); 
      outcome52[sub] = z52;

      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
    end;

    outcome1 = outcome51;
    outcome2 = outcome52;
    outcome3 = outcome13 + outcome23;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s11 < t* <= smax */

%macro inte213(indata, outdata, x1, x2);

  data temp72;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp72;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp72;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome51=j(obs,1,.);
    outcome52=j(obs,1,.);

    start fun51(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun52(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx3) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun33", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin };
      a2 = { &smin &s11};
      a3 = { &s11 &rmst };
      b = { .M  .P }; 

      call quad(z51,"fun51",b); 
      outcome51[sub] = z51;

      call quad(z52,"fun52",b); 
      outcome52[sub] = z52;

      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
    end;

    outcome1 = outcome51;
    outcome2 = outcome52;
    outcome3 = outcome13 + outcome23 + outcome33;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed t* > smax */

%macro inte214(indata, outdata, x1, x2);

  data temp72;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp72;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp72;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam10 = baseline[1, 1];
    gam11 = baseline[1, 2];
    gam12 = baseline[1, 3];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda11 = rcs[1, 3];
    k11 = rcs[1, 9];

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome43=j(obs,1,.);
    outcome51=j(obs,1,.);
    outcome52=j(obs,1,.);

    start fun51(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun52(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11); 
      s = &rmst;
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx4) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun33", b);
      return (u);
    finish;

    start fun43(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);;
      return(v); 
    finish; 

    start out43(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam30, gam31, gam32, k_min, k_max, lamda11, k11, s); 
      s = x;
      b = { .M  .P };  
      call quad(u, "fun43", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin };
      a2 = { &smin &s11};
      a3 = { &s11 &smax };
      a4 = { &smax &rmst };
      b = { .M  .P }; 

      call quad(z51,"fun51",b); 
      outcome51[sub] = z51;

      call quad(z52,"fun52",b); 
      outcome52[sub] = z52;

      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
      call quad(z43,"out43",a4); 
      outcome43[sub] = z43;
    end;

    outcome1 = outcome51;
    outcome2 = outcome52;
    outcome3 = outcome13 + outcome23 + outcome33 + outcome43;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* RCS DF = 3*/

/* Assessed t* <= smin */

%macro inte311(indata, outdata, x1, x2);

  data temp82;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	/* With baseline covariate specification */
    %if &stat = 1 %then %do;
	  %put 1;
      use temp82;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

	/* Without baseline covariate specification */
    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp82;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];

    k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outcome71=j(obs,1,.);
    outcome72=j(obs,1,.);

    outcome13=j(obs,1,.);

    start fun71(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun72(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) * exp(n) * (&dsx1) / s /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &rmst }; 
      b = { .M  .P }; 

      call quad(z71,"fun71",b); 
      outcome71[sub] = z71;
      call quad(z72,"fun72",b); 
      outcome72[sub] = z72;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
    end;

    outcome1 = outcome71;
    outcome2 = outcome72;
    outcome3 = outcome13;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed smin < t* <= s21 */

%macro inte312(indata, outdata, x1, x2);

  data temp82;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	/* With baseline covariate specification */
    %if &stat = 1 %then %do;
	  %put 1;
      use temp82;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

	/* Without baseline covariate specification */
    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp82;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];

    k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outcome71=j(obs,1,.);
    outcome72=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);

    start fun71(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun72(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) * exp(n) * (&dsx2) / s /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &rmst }; 
      b = { .M  .P }; 

      call quad(z71,"fun71",b); 
      outcome71[sub] = z71;
      call quad(z72,"fun72",b); 
      outcome72[sub] = z72;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
    end;

    outcome1 = outcome71;
    outcome2 = outcome72;
    outcome3 = outcome13 + outcome23;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s21 < t* <= s22 */

%macro inte313(indata, outdata, x1, x2);

  data temp82;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	/* With baseline covariate specification */
    %if &stat = 1 %then %do;
	  %put 1;
      use temp82;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

	/* Without baseline covariate specification */
    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp82;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];

    k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outcome71=j(obs,1,.);
    outcome72=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);

    start fun71(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun72(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) * exp(n) * (&dsx3) / s /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &s21 }; 
      a3 = { &s21 &rmst }; 
      b = { .M  .P }; 

      call quad(z71,"fun71",b); 
      outcome71[sub] = z71;
      call quad(z72,"fun72",b); 
      outcome72[sub] = z72;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
    end;

    outcome1 = outcome71;
    outcome2 = outcome72;
    outcome3 = outcome13 + outcome23 + outcome33;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed s22 < t* <= smax */

%macro inte314(indata, outdata, x1, x2);

  data temp82;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	/* With baseline covariate specification */
    %if &stat = 1 %then %do;
	  %put 1;
      use temp82;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

	/* Without baseline covariate specification */
    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp82;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];

    k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outcome71=j(obs,1,.);
    outcome72=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome43=j(obs,1,.);

    start fun71(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun72(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) * exp(n) * (&dsx4) / s /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    start fun43(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out43(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun43", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &s21 }; 
      a3 = { &s21 &s22 }; 
      a4 = { &s22  &rmst }; 
      b = { .M  .P }; 

      call quad(z71,"fun71",b); 
      outcome71[sub] = z71;
      call quad(z72,"fun72",b); 
      outcome72[sub] = z72;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
      call quad(z43,"out43",a4); 
      outcome43[sub] = z43;
    end;

    outcome1 = outcome71;
    outcome2 = outcome72;
    outcome3 = outcome13 + outcome23 + outcome33 + outcome43;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* Assessed t* > smax */

%macro inte315(indata, outdata, x1, x2);

  data temp82;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	/* With baseline covariate specification */
    %if &stat = 1 %then %do;
	  %put 1;
      use temp82;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

	/* Without baseline covariate specification */
    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp82;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam20 = baseline[1, 1];
    gam21 = baseline[1, 2];
    gam22 = baseline[1, 3];
    gam23 = baseline[1, 4];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda21 = rcs[1, 4];
    lamda22 = rcs[1, 5];

    k21 = rcs[1, 10];
    k22 = rcs[1, 11];

    outcome71=j(obs,1,.);
    outcome72=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome43=j(obs,1,.);
    outcome53=j(obs,1,.);

    start fun71(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun72(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22); 
      s = &rmst;
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) * exp(n) * (&dsx5) / s /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    start fun43(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out43(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun43", b);
      return (u);
    finish;

	start fun53(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v = exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out53(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda21, lamda22, k21, k22,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun53", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      a1 = { 0  &smin }; 
      a2 = { &smin &s21 }; 
      a3 = { &s21 &s22 }; 
      a4 = { &s22  &smax }; 
      a5 = { &smax &rmst };
      b = { .M  .P }; 

      call quad(z71,"fun71",b); 
      outcome71[sub] = z71;
      call quad(z72,"fun72",b); 
      outcome72[sub] = z72;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
      call quad(z43,"out43",a4); 
      outcome43[sub] = z43;
      call quad(z53,"out53",a5); 
      outcome53[sub] = z53;
    end;

    outcome1 = outcome71;
    outcome2 = outcome72;
    outcome3 = outcome13 + outcome23 + outcome33 + outcome43 + outcome53;

    postprobs=subjid||outcome1||outcome2||outcome3;
    cname = {"id" "S_T" "F_T" "RMST"};

    create &outdata from postprobs  [ colname=cname ];
    append from postprobs; 
  quit;

%mend;

/* RCS DF = 4*/

/* Assessed t* <= smin */

%macro inte411(indata, outdata, x1, x2);

  data temp92;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp92;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp92;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outcome91=j(obs,1,.);
    outcome92=j(obs,1,.);

    outcome13=j(obs,1,.);

    start fun91(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun92(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx1) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      b   = { .M  .P };
      a1   = { 0  &rmst }; 
      call quad(z91,"fun91",b); 
      outcome91[sub] = z91;
      call quad(z92,"fun92",b); 
      outcome92[sub] = z92;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
    end;

  outcome1 = outcome91;
  outcome2 = outcome92;
  outcome3 = outcome13;

  postprobs=subjid||outcome1||outcome2||outcome3;
  cname = {"id" "S_T" "F_T" "RMST"};

  create &outdata from postprobs  [ colname=cname ];
  append from postprobs; 
  quit;

%mend;

/* Assessed smin < t* <= s31 */

%macro inte412(indata, outdata, x1, x2);

  data temp92;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp92;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp92;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outcome91=j(obs,1,.);
    outcome92=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);

    start fun91(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun92(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx2) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      b   = { .M  .P };
      a1   = { 0  &smin }; 
      a2   = { &smin &rmst }; 
      call quad(z91,"fun91",b); 
      outcome91[sub] = z91;
      call quad(z92,"fun92",b); 
      outcome92[sub] = z92;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
    end;

  outcome1 = outcome91;
  outcome2 = outcome92;
  outcome3 = outcome13 + outcome23;

  postprobs=subjid||outcome1||outcome2||outcome3;
  cname = {"id" "S_T" "F_T" "RMST"};

  create &outdata from postprobs  [ colname=cname ];
  append from postprobs; 
  quit;

%mend;

/* Assessed s31 < t* <= s32 */

%macro inte413(indata, outdata, x1, x2);

  data temp92;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp92;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp92;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outcome91=j(obs,1,.);
    outcome92=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);

    start fun91(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun92(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx3) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      b   = { .M  .P };
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &rmst }; 
      call quad(z91,"fun91",b); 
      outcome91[sub] = z91;
      call quad(z92,"fun92",b); 
      outcome92[sub] = z92;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
    end;

  outcome1 = outcome91;
  outcome2 = outcome92;
  outcome3 = outcome13 + outcome23 + outcome33;

  postprobs=subjid||outcome1||outcome2||outcome3;
  cname = {"id" "S_T" "F_T" "RMST"};

  create &outdata from postprobs  [ colname=cname ];
  append from postprobs; 
  quit;

%mend;

/* Assessed s32 < t* <= s33 */

%macro inte414(indata, outdata, x1, x2);

  data temp92;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp92;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp92;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outcome91=j(obs,1,.);
    outcome92=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome43=j(obs,1,.);

    start fun91(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun92(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx4) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    start fun43(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out43(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun43", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      b   = { .M  .P };
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &s32 }; 
      a4   = { &s32 &rmst }; 
      call quad(z91,"fun91",b); 
      outcome91[sub] = z91;
      call quad(z92,"fun92",b); 
      outcome92[sub] = z92;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
      call quad(z43,"out43",a4); 
      outcome43[sub] = z43;
    end;

  outcome1 = outcome91;
  outcome2 = outcome92;
  outcome3 = outcome13 + outcome23 + outcome33 + outcome43;

  postprobs=subjid||outcome1||outcome2||outcome3;
  cname = {"id" "S_T" "F_T" "RMST"};

  create &outdata from postprobs  [ colname=cname ];
  append from postprobs; 
  quit;

%mend;

/* Assessed s33 < t* <= smax */

%macro inte415(indata, outdata, x1, x2);

  data temp92;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp92;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp92;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outcome91=j(obs,1,.);
    outcome92=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome43=j(obs,1,.);
    outcome53=j(obs,1,.);

    start fun91(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun92(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx5) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    start fun43(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out43(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun43", b);
      return (u);
    finish;

    start fun53(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out53(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun53", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      b   = { .M  .P };
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &s32 }; 
      a4   = { &s32 &s33 }; 
      a5   = { &s33  &rmst }; 
      call quad(z91,"fun91",b); 
      outcome91[sub] = z91;
      call quad(z92,"fun92",b); 
      outcome92[sub] = z92;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
      call quad(z43,"out43",a4); 
      outcome43[sub] = z43;
      call quad(z53,"out53",a5); 
      outcome53[sub] = z53;
    end;

  outcome1 = outcome91;
  outcome2 = outcome92;
  outcome3 = outcome13 + outcome23 + outcome33 + outcome43 + outcome53;

  postprobs=subjid||outcome1||outcome2||outcome3;
  cname = {"id" "S_T" "F_T" "RMST"};

  create &outdata from postprobs  [ colname=cname ];
  append from postprobs; 
  quit;

%mend;

/* Assessed t* > smax */

%macro inte416(indata, outdata, x1, x2);

  data temp92;
    set &indata;
  run;

  proc iml;
    use medpar2r;
    read all into med;
    use temp60r;
    read all into surv;
    use temp61r;
    read all into baseline;
    use temp13r;
    read all into rcs;

	%if &stat = 1 %then %do;
	  %put 1;
      use temp92;
      read all var{id} into subjid;
      read all var{&Cov} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0) || cov;
      medx = j(obs,1,1) || j(obs,1,&x2) || cov;
    %end;

    %else %if &stat = 2 %then %do;
	  %put 2;
      use temp92;
      read all var{id} into subjid;
      read all var{replicate} into cov;
      replicate = cov[,1];
      obs = nrow(replicate);
	  survx = j(obs,1,&x1) || j(obs,1,0);
      medx = j(obs,1,1) || j(obs,1,&x2);
    %end;

    sigma = med[1, ncol(med)];
    beta = med[1, 1:(ncol(med)-1)];
    alpha = surv[1, 1:ncol(surv)];
    alpha2 = alpha[1, 2];

    gam30 = baseline[1, 1];
    gam31 = baseline[1, 2];
    gam32 = baseline[1, 3];
    gam33 = baseline[1, 4];
    gam34 = baseline[1, 5];

    k_min = rcs[1, 1];
    k_max = rcs[1, 2];
    lamda31 = rcs[1, 6];
    lamda32 = rcs[1, 7];
    lamda33 = rcs[1, 8];

    k31 = rcs[1, 12];
    k32 = rcs[1, 13];
    k33 = rcs[1, 14];

    outcome91=j(obs,1,.);
    outcome92=j(obs,1,.);

    outcome13=j(obs,1,.);
    outcome23=j(obs,1,.);
    outcome33=j(obs,1,.);
    outcome43=j(obs,1,.);
    outcome53=j(obs,1,.);
    outcome63=j(obs,1,.);

    start fun91(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx6) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun92(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33); 
      s = &rmst;
      n = min(((&sx6) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n)) * exp(n) * (&dsx6) / s  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start fun13(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx1) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out13(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun13", b);
      return (u);
    finish;

    start fun23(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx2) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out23(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun23", b);
      return (u);
    finish;

    start fun33(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx3) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out33(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun33", b);
      return (u);
    finish;

    start fun43(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx4) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))/sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out43(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun43", b);
      return (u);
    finish;

    start fun53(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx5) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out53(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun53", b);
      return (u);
    finish;

	start fun63(m) global(beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      n = min(((&sx6) + sur*t(alpha) + alpha2 * m), 450);
      v =  exp((-1)*exp(n))  /sqrt(2*constant('pi')*sigma**2) * exp((-1)/(2 * sigma**2)*(m - (me*t(beta)))**2);
      return(v); 
    finish; 

    start out63(x) global (beta, alpha, alpha2, sigma, sur, me, gam10, gam11, gam12, gam20, gam21, gam22, gam23, gam30, gam31, gam32, gam33, gam34, k_min, k_max, lamda31, lamda32,lamda33,  k31, k32, k33,  s); 
      s = x;
      b = { .M  .P }; 
      call quad(u, "fun63", b);
      return (u);
    finish;

    do sub = 1 to obs;
      sur = survx[sub,];
      me = medx[sub,];
      b   = { .M  .P };
      a1   = { 0  &smin }; 
      a2   = { &smin &s31 }; 
      a3   = { &s31 &s32 }; 
      a4   = { &s32 &s33 }; 
      a5   = { &s33  &smax}; 
      a6 = { &smax &rmst };

      call quad(z91,"fun91",b); 
      outcome91[sub] = z91;
      call quad(z92,"fun92",b); 
      outcome92[sub] = z92;
      call quad(z13,"out13",a1); 
      outcome13[sub] = z13;
      call quad(z23,"out23",a2); 
      outcome23[sub] = z23;
      call quad(z33,"out33",a3); 
      outcome33[sub] = z33;
      call quad(z43,"out43",a4); 
      outcome43[sub] = z43;
      call quad(z53,"out53",a5); 
      outcome53[sub] = z53;
      call quad(z63,"out63",a6); 
      outcome63[sub] = z63;
    end;

  outcome1 = outcome91;
  outcome2 = outcome92;
  outcome3 = outcome13 + outcome23 + outcome33 + outcome43 + outcome53 + outcome63;

  postprobs=subjid||outcome1||outcome2||outcome3;
  cname = {"id" "S_T" "F_T" "RMST"};

  create &outdata from postprobs  [ colname=cname ];
  append from postprobs; 
  quit;

%mend;

%macro med();

  %do i = 1 %to (&nboot + 1);

/*    dm  'log;clear;out;clear;';*/

    data temp13r;
      set temp13;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data _dat32;
      set _dat05;
      if replicate = (&i - 1);
    run;

    data temp46r;
      set temp46;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data medpar2r;
      set medpar2;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data temp60r;
      set temp60;
      if replicate = (&i - 1);
      drop replicate;
    run;

    data _null_;
      set temp13r;
      call symput('smin', exp(k_min));
      call symput('smax', exp(k_max));
      call symput('s11', exp(k11));
      call symput('s21', exp(k21));
      call symput('s22', exp(k22));
      call symput('s31', exp(k31));
      call symput('s32', exp(k32));
      call symput('s33', exp(k33));
    run;

    data _null_;
      set temp46r;
      call symput('ID', ID);
    run;

	data temp14r;
	  set temp13r;
	  t_access = &t_access;
	  if t_access = 999999 then max_time = exp(k_max);
	  else max_time = t_access;
	  id = &id;
	  smin = &smin;
	  smax = &smax;
	  s11 = &s11;
	  s21 = &s21;
	  s22 = &s22;
	  s31 = &s31;
	  s32 = &s32;
	  s33 = &s33;
	  keep t_access max_time id smin smax s11 s21 s22 s31 s32 s33; 
	run;

	data temp15r;
	  set temp14r;
	  if id ge 1 and id le 44 then teind = 0;
	  if id = 45 then do;
	    if max_time le smin then teind = 11;
		else if max_time le s11 then teind = 12;
		else if max_time le smax then teind = 13;
		else teind = 14;
		if t_access = 999999 then teind = 13;
	  end;
	  if id = 46 then do;
	    if max_time le smin then teind = 21;
		else if max_time le s21 then teind = 22;
		else if max_time le s22 then teind = 23;
		else if max_time le smax then teind = 24;
		else teind = 25;
		if t_access = 999999 then teind = 24;
	  end;
	  if id = 47 then do;
	    if max_time le smin then teind = 31;
		else if max_time le s31 then teind = 32;
		else if max_time le s32 then teind = 33;
		else if max_time le s33 then teind = 34;
		else if max_time le smax then teind = 35;
		else teind = 36;
		if t_access = 999999 then teind = 35;
	  end;
    run;

	data _null_;
      set temp15r;
      call symput('rmst', max_time);
      call symput('teind', teind);
    run;

	%put &teind &rmst &smin &smax &s11 &s21 &s22 &s31 &s32 &s33;

    %if &id ge 1 and &id le 44 %then %do;

      %if &id ge 1 and &id le 8 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep lamd0 lamd1;
        run;

      %end;

      %if &id ge 9 and &id le 44 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep lamd0 lamd1 lamd2;
        run;

      %end;

      data _null_;
        format parameter parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 $250.;
        if &id = 1 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3)';
          call symput('dsx',parameter1);
        end;
        if &id = 2 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2)';
          call symput('dsx',parameter1);
        end;
        if &id = 3 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5)';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 4 then do;
          parameter = 'lamd0 + lamd1 * log(s)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s';
          call symput('dsx',parameter1);
        end;
        if &id = 5 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5)';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 6 then do;
          parameter = 'lamd0 + lamd1 * s ';
          call symput('sx',parameter);
          parameter1 = 'lamd1';
          call symput('dsx',parameter1);
        end;
        if &id = 7 then do;
          parameter = 'lamd0 + lamd1 * s **2';
          call symput('sx',parameter);
          parameter1 = '2 * lamd1 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 8 then do;
          parameter = 'lamd0 + lamd1 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '3 * lamd1 * s ** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 9 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** (-2) * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + lamd2 * s** (-3) + (-2) * lamd2 * s ** (-3) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 10 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** (-1)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + (-1) * lamd2 * s ** (-2)';
          call symput('dsx',parameter1);
        end;
        if &id = 11 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** (-0.5)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + (-0.5) * lamd2 * s** (-1.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 12 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + lamd2 /s';
          call symput('dsx',parameter1);
        end;
        if &id = 13 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** 0.5';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + (0.5) * lamd2 * s** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 14 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 15 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + 2*lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 16 then do;
          parameter = 'lamd0 + lamd1 * s ** (-2) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '(-2) * lamd1 * s ** (-3) + 3 * lamd2 * s** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 17 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** (-1) * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + lamd2 * s** (-2) + (-1) * lamd2 * s ** (-2) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 18 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** (-0.5)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + (-0.5) * lamd2 * s ** (-1.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 19 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + lamd2 /s';
          call symput('dsx',parameter1);
        end;
        if &id = 20 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** 0.5';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + 0.5 * lamd2 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 21 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 22 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 23 then do;
          parameter = 'lamd0 + lamd1 * s ** (-1) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '(-1) * lamd1 * s ** (-2) + 3 * lamd2 * s** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 24 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** (-0.5) * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + lamd2 * s** (-1.5) + (-0.5) * lamd2 * s ** (-1.5) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 25 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + lamd2 / s';
          call symput('dsx',parameter1);
        end;
        if &id = 26 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** 0.5';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + 0.5 * lamd2 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 27 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 28 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 29 then do;
          parameter = 'lamd0 + lamd1 * s ** (-0.5) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '(-0.5) * lamd1 * s ** (-1.5) + 3 * lamd2 * s ** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 30 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * log(s) * log(s)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + 2 * lamd2 * log(s)/s';
          call symput('dsx',parameter1);
        end;
        if &id = 31 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s ** (0.5)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + 0.5 * lamd2 * s ** (-0.5)';
          call symput('dsx',parameter1);
        end;
        if &id = 32 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s ';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + lamd2 ';
          call symput('dsx',parameter1);
        end;
        if &id = 33 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s  + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 34 then do;
          parameter = 'lamd0 + lamd1 * log(s) + lamd2 * s**3';
          call symput('sx',parameter);
          parameter1 = 'lamd1 / s + lamd2 * s**2 * 3';
          call symput('dsx',parameter1);
        end;

        if &id = 35 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s ** (0.5) * log(s)';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + lamd2 * s** (-0.5) + 0.5 * lamd2 * s ** (-0.5) * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 36 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + lamd2';
          call symput('dsx',parameter1);
        end;
        if &id = 37 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 38 then do;
          parameter = 'lamd0 + lamd1 * s ** (0.5) + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = '0.5 * lamd1 * s ** (-0.5) + 3 * lamd2 * s ** 2';
          call symput('dsx',parameter1);
        end;

        if &id = 39 then do;
          parameter = 'lamd0 + lamd1 * s + lamd2 * s * log(s)';
          call symput('sx',parameter);
          parameter1 = 'lamd1 + lamd2 + lamd2 * log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 40 then do;
          parameter = 'lamd0 + lamd1 * s + lamd2 * s ** 2';
          call symput('sx',parameter);
          parameter1 = 'lamd1 + 2 * lamd2 * s';
          call symput('dsx',parameter1);
        end;
        if &id = 41 then do;
          parameter = 'lamd0 + lamd1 * s + lamd2 * s ** 3';
          call symput('sx',parameter);
          parameter1 = 'lamd1 + 3 * lamd2 * s**2';
          call symput('dsx',parameter1);
        end;

        if &id = 42 then do;
          parameter = 'lamd0 + lamd1 * s**2 + lamd2 * s**2 * log(s)';
          call symput('sx',parameter);
          parameter1 = '2*lamd1*s + lamd2*s + 2*lamd2 * s* log(s)';
          call symput('dsx',parameter1);
        end;
        if &id = 43 then do;
          parameter = 'lamd0 + lamd1 * s**2 + lamd2 * s**3';
          call symput('sx',parameter);
          parameter1 = '2*lamd1*s + 3*lamd2*s**2';
          call symput('dsx',parameter1);
        end;

        if &id = 44 then do;
          parameter = 'lamd0 + lamd1 * s**3 + lamd2 * s**3 * log(s)';
          call symput('sx',parameter);
          parameter1 = '3*lamd1*s**2 + lamd2*s**2 + 3*lamd2 * s**2* log(s)';
          call symput('dsx',parameter1);
        end;
      run;

      %inte11(_dat32, final100, 0, 0);
      %inte11(_dat32, final101, 0, 1);
      %inte11(_dat32, final110, 1, 0);
      %inte11(_dat32, final111, 1, 1);

    %end;

    %if &id = 45 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep gam10 gam11 gam12;
        run;

        data _null_;
          format parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 $250.;
          parameter1 = 'gam10 + gam11 * log(s)';
          call symput('sx1',parameter1);
          parameter2 = 'gam11';
          call symput('dsx1',parameter2);
          parameter3 = 'gam10 + gam11 * log(s) + gam12 * (- lamda11 * (log(s) - k_min)**3)';
          call symput('sx2',parameter3);
          parameter4 = 'gam11 + gam12 * (- lamda11 * (log(s) - k_min)**2*3)';
          call symput('dsx2',parameter4);
          parameter5 = 'gam10 + gam11 * log(s) + gam12 * ((log(s) - k11)**3 - lamda11 * (log(s) - k_min)**3)';
          call symput('sx3',parameter5);
          parameter6 = 'gam11 + gam12 * ((log(s) - k11)**2*3 - lamda11 * (log(s) - k_min)**2*3)';
          call symput('dsx3',parameter6);
          parameter7 = 'gam10 + gam11 * log(s) + gam12 * ((log(s) - k11)**3 - lamda11 * (log(s) - k_min)**3 - (1 - lamda11) * (log(s) - k_max)**3)';
          call symput('sx4',parameter7);
          parameter8 = 'gam11 + gam12 * ((log(s) - k11)**2*3 - lamda11 * (log(s) - k_min)**2*3 - (1 - lamda11) * (log(s) - k_max)**2*3)';
          call symput('dsx4',parameter8);
        run;

		%put &rmst &smin &smax &s11 &s21 &s22 &s31 &s32 &s33;

        %if &teind = 11 %then %do;

          %inte211(_dat32, final100, 0, 0);
          %inte211(_dat32, final101, 0, 1);
          %inte211(_dat32, final110, 1, 0);
          %inte211(_dat32, final111, 1, 1);

        %end;

        %else %if &teind = 12 %then %do;

          %inte212(_dat32, final100, 0, 0);
          %inte212(_dat32, final101, 0, 1);
          %inte212(_dat32, final110, 1, 0);
          %inte212(_dat32, final111, 1, 1);

        %end;

		%else %if &teind = 13 %then %do;

          %inte213(_dat32, final100, 0, 0);
          %inte213(_dat32, final101, 0, 1);
          %inte213(_dat32, final110, 1, 0);
          %inte213(_dat32, final111, 1, 1);

        %end;

    	%else %if &teind = 14 %then %do;

          %inte214(_dat32, final100, 0, 0);
          %inte214(_dat32, final101, 0, 1);
          %inte214(_dat32, final110, 1, 0);
          %inte214(_dat32, final111, 1, 1);

        %end;

      %end;

      %if &id = 46 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep gam20 gam21 gam22 gam23;
        run;

        data _null_;
          format parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 parameter9 parameter10 $500.;
          parameter1 = 'gam20 + gam21 * log(s)';
          call symput('sx1',parameter1);
          parameter2 = 'gam21';
          call symput('dsx1',parameter2);

          parameter3 = 'gam20 + gam21 * log(s) - gam22 * (lamda21 * (log(s) - k_min)**3) - gam23 * (lamda22 * (log(s) - k_min)**3)';
          call symput('sx2',parameter3);
          parameter4 = 'gam21 + gam22 * (- lamda21 * (log(s) - k_min)**2*3) + gam23 * (- lamda22 * (log(s) - k_min)**2*3)';
          call symput('dsx2',parameter4);

          parameter5 = 'gam20 + gam21 * log(s) + gam22 * ((log(s) - k21)**3 - lamda21 * (log(s) - k_min)**3) - gam23 * (lamda22 * (log(s) - k_min)**3)';
          call symput('sx3',parameter5);
          parameter6 = 'gam21 + gam22 * ((log(s) - k21)**2*3 - lamda21 * (log(s) - k_min)**2*3) - gam23 * (3 * lamda22 * (log(s) - k_min)**2)';
          call symput('dsx3',parameter6);

          parameter7 = 'gam20 + gam21 * log(s) + gam22 * ((log(s) - k21)**3 - lamda21 * (log(s) - k_min)**3) + gam23 * ((log(s) - k22)**3 - lamda22 * (log(s) - k_min)**3)';
          call symput('sx4',parameter7);
          parameter8 = 'gam21 + gam22 * ((log(s) - k21)**2*3 - lamda21 * (log(s) - k_min)**2*3) + gam23 * (3 * (log(s) - k22)**2 - 3 * lamda22 * (log(s) - k_min)**2)';
          call symput('dsx4',parameter8);

          parameter9 = 'gam20 + gam21 * log(s) + gam22 * ((log(s) - k21)**3 - lamda21 * (log(s) - k_min)**3 - (1 - lamda21) * (log(s) - k_max)**3) 
          + gam23 * ((log(s) - k22)**3 - lamda22 * (log(s) - k_min)**3 - (1 - lamda22) * (log(s) - k_max)**3)';
          call symput('sx5',parameter9);
          parameter10 = 'gam21 + gam22 * ((log(s) - k21)**2*3 - lamda21 * (log(s) - k_min)**2*3 - (1 - lamda21) * (log(s) - k_max)**2*3) 
          + gam23 * ((log(s) - k22)**2*3 - lamda22 * (log(s) - k_min)**2*3 - (1 - lamda22) * (log(s) - k_max)**2*3)';
          call symput('dsx5',parameter10);

        run;

		%put &rmst &smin &smax &s11 &s21 &s22 &s31 &s32 &s33;

		%if &teind = 21 %then %do;

          %inte311(_dat32, final100, 0, 0);
          %inte311(_dat32, final101, 0, 1);
          %inte311(_dat32, final110, 1, 0);
          %inte311(_dat32, final111, 1, 1);

        %end;

        %else %if &teind = 22 %then %do;

          %inte312(_dat32, final100, 0, 0);
          %inte312(_dat32, final101, 0, 1);
          %inte312(_dat32, final110, 1, 0);
          %inte312(_dat32, final111, 1, 1);

        %end;

		%else %if &teind = 23 %then %do;

          %inte313(_dat32, final100, 0, 0);
          %inte313(_dat32, final101, 0, 1);
          %inte313(_dat32, final110, 1, 0);
          %inte313(_dat32, final111, 1, 1);

        %end;

    	%else %if &teind = 24 %then %do;

          %inte314(_dat32, final100, 0, 0);
          %inte314(_dat32, final101, 0, 1);
          %inte314(_dat32, final110, 1, 0);
          %inte314(_dat32, final111, 1, 1);

        %end;

		%else %if &teind = 25 %then %do;

          %inte315(_dat32, final100, 0, 0);
          %inte315(_dat32, final101, 0, 1);
          %inte315(_dat32, final110, 1, 0);
          %inte315(_dat32, final111, 1, 1);

        %end;

      %end;

      %if &id = 47 %then %do;

        data temp61r;
          set temp61;
          if replicate = (&i - 1);
          keep gam30 gam31 gam32 gam33 gam34;
        run;

        data _null_;
          format parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 parameter7 parameter8 parameter9 parameter10 parameter11 parameter12 $500.;
          parameter1 = 'gam30 + gam31 * log(s)';
          call symput('sx1',parameter1);
          parameter2 = 'gam31';
          call symput('dsx1',parameter2);

          parameter3 = 'gam30 + gam31 * log(s) + gam32 * (- lamda31 * (log(s) - k_min)**3) + gam33 * (-lamda32 * (log(s) - k_min)**3) + gam34 * (-lamda33 * (log(s) - k_min)**3)';
          call symput('sx2',parameter3);
          parameter4 = 'gam31 + gam32 * (- lamda31 * (log(s) - k_min)**2*3) + gam33 * (- lamda32 * (log(s) - k_min)**2*3)+ gam34 * (- lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx2',parameter4);

          parameter5 = 'gam30 + gam31 * log(s) + gam32 * ((log(s) - k31)**3 - lamda31 * (log(s) - k_min)**3) + gam33 * (-lamda32 * (log(s) - k_min)**3)+ gam34 * (-lamda33 * (log(s) - k_min)**3)';
          call symput('sx3',parameter5);
          parameter6 = 'gam31 + gam32 * ((log(s) - k31)**2*3 - lamda31 * (log(s) - k_min)**2*3) + gam33 * (-lamda32 * (log(s) - k_min)**2*3) + gam34 * (-lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx3',parameter6);

          parameter7 = 'gam30 + gam31 * log(s) + gam32 * ((log(s) - k31)**3 - lamda31 * (log(s) - k_min)**3) + gam33 * ((log(s) - k32)**3 - lamda32 * (log(s) - k_min)**3) + gam34 * (-lamda33 * (log(s) - k_min)**3)';
          call symput('sx4',parameter7);
          parameter8 = 'gam31 + gam32 * ((log(s) - k31)**2*3 - lamda31 * (log(s) - k_min)**2*3) + gam33 * (3 * (log(s) - k32)**2 - 3 * lamda32 * (log(s) - k_min)**2) + gam34 * (-lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx4',parameter8);

          parameter9 = 'gam30 + gam31 * log(s) + gam32 * ((log(s) - k31)**3-lamda31*(log(s)-k_min)**3)+gam33 *((log(s)-k32)**3-lamda32 *(log(s)-k_min)**3)+gam34 *((log(s) - k33)**3-lamda33*(log(s)-k_min)**3)';
          call symput('sx5',parameter9);
          parameter10 = 'gam31+gam32*((log(s)-k31)**2*3-lamda31*(log(s)-k_min)**2*3)+gam33*(3 * (log(s) - k32)**2 - 3 * lamda32 * (log(s) - k_min)**2) + gam34 * (3 * (log(s) - k33)**2 -lamda33 * (log(s) - k_min)**2*3)';
          call symput('dsx5',parameter10);

		  parameter11 = 'gam30+gam31*log(s)+gam32*((log(s)-k31)**3-lamda31*(log(s)-k_min)**3-(1-lamda31)*(log(s)-k_max)**3)+gam33*((log(s)-k32)**3-lamda32*(log(s)-k_min)**3-(1-lamda32)*(log(s)-k_max)**3)+gam34*((log(s)-k33)**3-lamda33*(log(s)-k_min)**3-(1-lamda33)*(log(s)-k_max)**3)';
          call symput('sx6',parameter11);
          parameter12 = 'gam31+gam32*((log(s)-k31)**2*3-lamda31*(log(s)-k_min)**2*3-(1-lamda31)*(log(s)-k_max)**2*3)+gam33*((log(s)-k32)**2*3-lamda32*(log(s)-k_min)**2*3-(1-lamda32)*(log(s)-k_max)**2*3)+gam34*((log(s)-k33)**2*3-lamda33*(log(s)-k_min)**2*3-(1-lamda33)*(log(s)-k_max)**2*3)';
          call symput('dsx6',parameter12);
        run;

		%put &sx6;
		%put &rmst &smin &smax &s11 &s21 &s22 &s31 &s32 &s33;

		%if &teind = 31 %then %do;

          %inte411(_dat32, final100, 0, 0);
          %inte411(_dat32, final101, 0, 1);
          %inte411(_dat32, final110, 1, 0);
          %inte411(_dat32, final111, 1, 1);

        %end;

        %else %if &teind = 32 %then %do;

          %inte412(_dat32, final100, 0, 0);
          %inte412(_dat32, final101, 0, 1);
          %inte412(_dat32, final110, 1, 0);
          %inte412(_dat32, final111, 1, 1);

        %end;

		%else %if &teind = 33 %then %do;

          %inte413(_dat32, final100, 0, 0);
          %inte413(_dat32, final101, 0, 1);
          %inte413(_dat32, final110, 1, 0);
          %inte413(_dat32, final111, 1, 1);

        %end;

    	%else %if &teind = 34 %then %do;

          %inte414(_dat32, final100, 0, 0);
          %inte414(_dat32, final101, 0, 1);
          %inte414(_dat32, final110, 1, 0);
          %inte414(_dat32, final111, 1, 1);

        %end;

		%else %if &teind = 35 %then %do;

          %inte415(_dat32, final100, 0, 0);
          %inte415(_dat32, final101, 0, 1);
          %inte415(_dat32, final110, 1, 0);
          %inte415(_dat32, final111, 1, 1);

        %end;

		%else %if &teind = 36 %then %do;

          %inte416(_dat32, final100, 0, 0);
          %inte416(_dat32, final101, 0, 1);
          %inte416(_dat32, final110, 1, 0);
          %inte416(_dat32, final111, 1, 1);

        %end;

      %end;

      data final;
        set final final100 (in = a) final101 (in = b) final110 (in = c) final111 (in = d);
        if a then ind = 1;
        if b then ind = 2;
        if c then ind = 3;
        if d then ind = 4;
        if a then replicate = &i - 1;
        if b then replicate = &i - 1;
        if c then replicate = &i - 1;
        if d then replicate = &i - 1;
        hazard = F_T/S_T;
      run;

      proc sort data = final;
        by replicate id;
      run;

      proc datasets lib = work;
        delete final100 final101 final110 final111 temp13r temp14r temp15r temp60r temp61r temp46r medpar2r _dat32;
      run;

  %end;

%mend;

*******************************************************************************************************
                                             Main Macro;
*******************************************************************************************************;

%macro SURVPHMED (data             = '', 
                  SURVtime         = '',
                  SURVevent        = '',
                  MED              = '', 
                  Cov              = '', 
                  EXP              = '',
                  Nboot            = 200,
				  T_ACCESS         = 999999,
				  REFERENCE        = 1,
                  Nlmopts          = ,
                  Printall         = T,
                  Out              = out01
                  );

data zzz0;
  indata = "&data";
  survtime = "&survtime";
  survevent = "&survevent";
  med = "&med";
  cov = "&cov";
  exp = "&exp";
  reference = &reference;
  if indata = '' or SURVtime=' ' or SURVevent=' ' or MED=' ' or EXP = ' ' 
  or compress(indata, '') = "''" or compress(SURVtime, '') = "''" or compress(SURVevent, '') = "''" or compress(MED, '') = "''" or compress(EXP, '') = "''"  
  then sta = 0; 
  else if (cov = '' or compress(cov, '') = "''") and reference ne 3 then sta = 2; 
  else if (cov = '' or compress(cov, '') = "''") and reference = 3 then sta = 3; 
  else sta = 1;
run;

data _null_;
  set zzz0;
  call symput('stat', sta);
run;

/* Make sure input data set name, outcome, mediator and expsoure variables are not missing 
stat = 0: miss data set name, exposure, mediator or outcome variables
stat = 1: cov exist
stat = 2: cov does not exist and reference is not 3
stat = 3: cov does not exist and reference is 3*/

%if &stat = 0 %then %do;
  proc iml;
    print,"WARNING: NEED TO SPECIFY DATASET NAME, EXPOSURE, MEDIATOR AND OUTCOME VARIABLE","PROGRAM WILL TERMINATE";
  quit;
%end;

%else %if &stat = 3 %then %do;
  proc iml;
    print,"WARNING: NO COVARIATES FOR THE SURVIVAL OUTCOME AND THE MEDIATOR, HOWEVER REFERENCE IS SPECIFIED AS 3, CONFLICT EXISTS FOR VARIABLE SPECIFICATION","PROGRAM WILL TERMINATE";;
  quit;
%end;

%else %do;

ods listing close;

proc surveyselect data= &data out=boots
     method = urs
	 samprate = 1 outhits rep = &Nboot seed = 7345678;
run;

data zzz2;
set &data;
replicate = 0;
run;

data bootsample1;
retain replicate;
set zzz2 boots;
drop NumberHits;
run;

*******************************************************************************************************
Determine Parametric Survival Model;
*******************************************************************************************************;

/* RCS */

data _dat; 
  set bootsample1; 
  _time = &SURVtime;
  _event = &SURVevent;
  lnt = log(_time);
  i = 1;
run;

data _dat1;
  set _dat;
  if _event = 1;
run;

proc univariate data = _dat1 noprint;
  var lnt;
  output out=temp11 pctlpts=50 33 67 25 75 pctlpre = k pctlname = _50 _33 _67 _25 _75;
  by replicate;
run;

proc means data = _dat1 noprint;
  var lnt;
  output out = temp12 min(lnt) = k_min max(lnt) = k_max;
  by replicate;
run;

data temp13;
  merge temp11 temp12;
  by replicate;
  drop _type_ _freq_;
  lamda11 = (k_max - k_50)/(k_max - k_min);
  lamda21 = (k_max - k_33)/(k_max - k_min);
  lamda22 = (k_max - k_67)/(k_max - k_min);
  lamda31 = (k_max - k_25)/(k_max - k_min);
  lamda32 = (k_max - k_50)/(k_max - k_min);
  lamda33 = (k_max - k_75)/(k_max - k_min);
  k11 = k_50;
  k21 = k_33;
  k22 = k_67;
  k31 = k_25;
  k32 = k_50;
  k33 = k_75;
  i = 1;
  drop k_50 k_33 k_67 k_25 k_75;
run;

data _dat11;
  merge _dat temp13;
  by replicate i;

  if lnt le k_min then do;
    nu11 = 0;
    dnu11 = 0;
  end;
  else if lnt le k11 then do;
    nu11 = - lamda11 * (lnt - k_min)**3;
    dnu11 = - lamda11 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu11 = (lnt - k11)**3 - lamda11 * (lnt - k_min)**3; 
    dnu11 = (lnt - k11)**2*3 - lamda11 * (lnt - k_min)**2*3; 
  end;
  else do;
    nu11 = (lnt - k11)**3 - lamda11 * (lnt - k_min)**3 - (1 - lamda11) * (lnt - k_max)**3; 
    dnu11 = (lnt - k11)**2*3 - lamda11 * (lnt - k_min)**2*3 - (1 - lamda11) * (lnt - k_max)**2*3; 
  end;

  if lnt le k_min then do;
    nu21 = 0;
    dnu21 = 0;
    nu22 = 0;
    dnu22 = 0;
  end;
  else if lnt le k21 then do;
    nu21 = - lamda21 * (lnt - k_min)**3;
    nu22 = - lamda22 * (lnt - k_min)**3;
    dnu21 = - lamda21 * (lnt - k_min)**2*3;
    dnu22 = - lamda22 * (lnt - k_min)**2*3;
  end;
  else if lnt le k22 then do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3;
    nu22 =  - lamda22 * (lnt - k_min)**3;
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3;
    dnu22 =  - lamda22 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3;
    nu22 = (lnt - k22)**3 - lamda22 * (lnt - k_min)**3;
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3;
    dnu22 = (lnt - k22)**2*3 - lamda22 * (lnt - k_min)**2*3;
  end;
  else do;
    nu21 = (lnt - k21)**3 - lamda21 * (lnt - k_min)**3 - (1 - lamda21) * (lnt - k_max)**3; 
    dnu21 = (lnt - k21)**2*3 - lamda21 * (lnt - k_min)**2*3 - (1 - lamda21) * (lnt - k_max)**2*3; 
    nu22 = (lnt - k22)**3 - lamda22 * (lnt - k_min)**3 - (1 - lamda22) * (lnt - k_max)**3; 
    dnu22 = (lnt - k22)**2*3 - lamda22 * (lnt - k_min)**2*3 - (1 - lamda22) * (lnt - k_max)**2*3; 
  end;

  if lnt le k_min then do;
    nu31 = 0;
    dnu31 = 0;
    nu32 = 0;
    dnu32 = 0;
    nu33 = 0;
    dnu33 = 0;
  end;
  else if lnt le k31 then do;
    nu31 = - lamda31 * (lnt - k_min)**3;
    nu32 = - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = - lamda31 * (lnt - k_min)**2*3;
    dnu32 = - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k32 then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k33 then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3;
    nu33 = - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3;
    dnu33 = - lamda33 * (lnt - k_min)**2*3;
  end;
  else if lnt le k_max then do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3;
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3;
    nu33 = (lnt - k33)**3 - lamda33 * (lnt - k_min)**3;
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3;
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3;
    dnu33 = (lnt - k33)**2*3 - lamda33 * (lnt - k_min)**2*3;
  end;
  else do;
    nu31 = (lnt - k31)**3 - lamda31 * (lnt - k_min)**3 - (1 - lamda31) * (lnt - k_max)**3; 
    dnu31 = (lnt - k31)**2*3 - lamda31 * (lnt - k_min)**2*3 - (1 - lamda31) * (lnt - k_max)**2*3; 
    nu32 = (lnt - k32)**3 - lamda32 * (lnt - k_min)**3 - (1 - lamda32) * (lnt - k_max)**3; 
    dnu32 = (lnt - k32)**2*3 - lamda32 * (lnt - k_min)**2*3 - (1 - lamda32) * (lnt - k_max)**2*3; 
    nu33 = (lnt - k33)**3 - lamda33 * (lnt - k_min)**3 - (1 - lamda33) * (lnt - k_max)**3; 
    dnu33 = (lnt - k33)**2*3 - lamda33 * (lnt - k_min)**2*3 - (1 - lamda33) * (lnt - k_max)**2*3; 
  end;

  drop i;
run;

data _dat21;
  set _dat11;
  if _event = 1;
 * drop age -- months _event;
  drop &Cov &exp _event;
run;

proc phreg data=_dat;
*  class male (ref="0" param=ref) idealhealthsmk (ref="0" param=ref) agec (ref="0" param=ref) bmic (ref="0" param=ref);
  class &exp (ref="0" param=ref);
  model _time*_event(0)= &exp &med &Cov;
  ods output ParameterEstimates=SURVests1(rename=(Parameter=Effect));
  baseline out = bb5 survival = _all_ CUMHAZ=_ALL_;
  by replicate;
run;

data tt51;
  set survests1;
  if replicate = 0;
  i = 1;
  np1 = _n_;
run;

data tt52;
  set tt51;
  by i;
  if last.i;
  keep np1;
run; 

data _null_;
  set tt52;
  call symput('NP1', np1);
run;

*Survival model parameters;
data SURVests2; 
  length model $4 parameter $8; 
  set SURVests1;
  format row 3.0;
  row = (_n_ - replicate * &np1)-1;
  parameter = "a" || left(row);
  if DF=0 then delete; *remove unestimated parameters;
  model="Surv";
  if _n_ ne 1 then linpred =  "+" || cats(parameter) || "*" || effect;
  else linpred =  cats(parameter) || "*" || effect;
  i = 1;
run; 

data SURVests3; 
  set SURVests2;
  if replicate = 0;
  if _N_=1 then call symput('linpredSURV',linpred);
  else call symput('linpredSURV',trim(resolve('&linpredSURV'))||linpred);
run;
%put &linpredSURV;
title;

proc sort data = _dat21;
  by replicate _time;
run;

data bb51;
  set bb5;
  if _time ne 0;
  keep replicate &exp &med &Cov _time CumHaz i;
  i = 1;
run;

proc sort data = bb51;
  by replicate _time;
run;

proc transpose data = SURVests2 out = survests4;
  var estimate;
  id parameter;
  by replicate i;
run;

proc phreg data=_dat;
*  class male (ref="0" param=ref) idealhealthsmk (ref="0" param=ref) agec (ref="0" param=ref) bmic (ref="0" param=ref);
  class &exp (ref="0" param=ref);
  model _time*_event(0)= &exp &Cov;
  ods output ParameterEstimates=ests1(rename=(Parameter=Effect));
  by replicate;
run;

data ests2; 
  length model $4 parameter $8; 
  set ests1;
  format row 3.0;
  row = (_n_ - replicate * (&np1 -1))-1;
  parameter = "c" || left(row);
  if DF=0 then delete; *remove unestimated parameters;
  model="Surv";
  if _n_ ne 1 then linpred =  "+" || cats(parameter) || "*" || effect;
  else linpred =  cats(parameter) || "*" || effect;
  i = 1;
run; 

proc transpose data = ests2 out = ests4;
  var estimate;
  id parameter;
  by replicate i;
run;

data bb53;
  merge survests4 bb51;
  by replicate i;
  logbaseCumHaz = log(CumHaz) - (&linpredSURV);
  keep replicate _time logbaseCumHaz;
run;

data _dat22;
  merge _dat21 bb53;
  by replicate _time;
run;

proc glm data = _dat22;
  model logbasecumhaz = lnt nu11;
  by replicate;
  ods output ParameterEstimates  = tt11;
run;

proc transpose data = tt11 out = tt12 (drop = _name_ rename = (col1 = gam10 col2 = gam11 col3 = gam12));
  by replicate;
  var estimate;
run;

proc glm data = _dat22;
  model logbasecumhaz = lnt nu21 nu22;
  by replicate;
  ods output ParameterEstimates  = tt21;
run;

proc transpose data = tt21 out = tt22 (drop = _name_ rename = (col1 = gam20 col2 = gam21 col3 = gam22 col4 = gam23));
  by replicate;
  var estimate;
run;

proc glm data = _dat22;
  model logbasecumhaz = lnt nu31 nu32 nu33;
  by replicate;
  ods output ParameterEstimates  = tt31;
run;

proc transpose data = tt31 out = tt32 (drop = _name_ rename = (col1 = gam30 col2 = gam31 col3 = gam32 col4 = gam33 col5 = gam34));
  by replicate;
  var estimate;
run;

data bb54;
  merge tt12 tt22 tt32;
  by replicate;
run;

proc transpose data = bb54 out = bb55 (rename = (_name_ = parameter col1 = estimate));
  var gam10 -- gam34;
  by replicate;
run;

data SurvInit;
  set SURVests2 (keep = replicate model parameter estimate effect) bb55 (in = a);
  if a then model = 'RCS';
run;

data SurvInit;
  set SurvInit;
  if parameter in ('gam10','gam20', 'gam30') then effect = 'Intercept';
  if parameter in ('gam11','gam21', 'gam31') then effect = 'lnt';
  if parameter in ('gam12') then effect = 'nu11';
  if parameter in ('gam22') then effect = 'nu21';
  if parameter in ('gam23') then effect = 'nu22';
  if parameter in ('gam32') then effect = 'nu31';
  if parameter in ('gam33') then effect = 'nu32';
  if parameter in ('gam34') then effect = 'nu33';
run;

proc sort data = survinit;
  by replicate;
run;

data SurvInit1;
  set SurvInit;
  if substr(parameter, 1, 4) in ('gam2', 'gam3') then delete;
run; 

data SurvInit2;
  set SurvInit;
  if substr(parameter, 1, 4) in ('gam1', 'gam3') then delete;
run; 

data SurvInit3;
  set SurvInit;
  if substr(parameter, 1, 4) in ('gam2', 'gam1') then delete;
run; 

%rcs;

/* Fractional Polynomials */

data _dat31;
  set _dat;
  t1 = _time ** (-2);
  t2 = _time ** (-1);
  t3 = _time ** (-0.5);
  t4 = log(_time);
  t5 = _time ** (0.5);
  t6 = _time;
  t7 = _time ** 2;
  t8 = _time ** 3;

  dt1 = _time ** (-3) * (-2);
  dt2 = _time ** (-2) * (-1);
  dt3 = _time ** (-1.5) * (-0.5);
  dt4 = 1/_time;
  dt5 = _time ** (-0.5) * (0.5);
  dt6 = 1;
  dt7 = _time * 2;
  dt8 = _time ** 2 * 3;

  t9 = _time ** (-2) * log(_time);
  t10 = _time ** (-1) * log(_time);
  t11 = _time ** (-0.5) * log(_time);
  t12 = log(_time) * log(_time);
  t13 = _time ** (0.5) * log(_time);
  t14 = _time * log(_time);
  t15 = _time ** 2 * log(_time);
  t16 = _time ** 3 * log(_time);

  dt9 = _time ** (-3) + _time ** (-3) * (-2) * log(_time);
  dt10 = _time ** (-2) + _time ** (-2) * (-1) * log(_time);
  dt11 = _time ** (-1.5) + _time ** (-1.5) * (-0.5) * log(_time);
  dt12 = 2 * log(_time) /_time;
  dt13 = _time ** (-0.5) + _time ** (-0.5) * (0.5) * log(_time);
  dt14 = 1 + log(_time);
  dt15 = _time + _time * 2 * log(_time);
  dt16 = _time ** 2 + _time ** 2 * 3 * log(_time);

  drop i;
run;

data bb56;
  set bb53;
  t1 = _time ** (-2);
  t2 = _time ** (-1);
  t3 = _time ** (-0.5);
  t4 = log(_time);
  t5 = _time ** (0.5);
  t6 = _time;
  t7 = _time ** 2;
  t8 = _time ** 3;

  dt1 = _time ** (-3) * (-2);
  dt2 = _time ** (-2) * (-1);
  dt3 = _time ** (-1.5) * (-0.5);
  dt4 = 1/_time;
  dt5 = _time ** (-0.5) * (0.5);
  dt6 = 1;
  dt7 = _time * 2;
  dt8 = _time ** 2 * 3;

  t9 = _time ** (-2) * log(_time);
  t10 = _time ** (-1) * log(_time);
  t11 = _time ** (-0.5) * log(_time);
  t12 = log(_time) * log(_time);
  t13 = _time ** (0.5) * log(_time);
  t14 = _time * log(_time);
  t15 = _time ** 2 * log(_time);
  t16 = _time ** 3 * log(_time);

  dt9 = _time ** (-3) + _time ** (-3) * (-2) * log(_time);
  dt10 = _time ** (-2) + _time ** (-2) * (-1) * log(_time);
  dt11 = _time ** (-1.5) + _time ** (-1.5) * (-0.5) * log(_time);
  dt12 = 2 * log(_time) /_time;
  dt13 = _time ** (-0.5) + _time ** (-0.5) * (0.5) * log(_time);
  dt14 = 1 + log(_time);
  dt15 = _time + _time * 2 * log(_time);
  dt16 = _time ** 2 + _time ** 2 * 3 * log(_time);

run;

%single(_dat31, bb56, ff01, hh01, gg01, 1, 1);
%single(_dat31, bb56, ff02, hh02, gg02, 2, 2);
%single(_dat31, bb56, ff03, hh03, gg03, 3, 3);
%single(_dat31, bb56, ff04, hh04, gg04, 4, 4);
%single(_dat31, bb56, ff05, hh05, gg05, 5, 5);
%single(_dat31, bb56, ff06, hh06, gg06, 6, 6);
%single(_dat31, bb56, ff07, hh07, gg07, 7, 7);
%single(_dat31, bb56, ff08, hh08, gg08, 8, 8);


%double(_dat31, bb56, ff09, hh09, gg09, 1, 9, 9);
%double(_dat31, bb56, ff10, hh10, gg10, 1, 2, 10);
%double(_dat31, bb56, ff11, hh11, gg11, 1, 3, 11);
%double(_dat31, bb56, ff12, hh12, gg12, 1, 4, 12);
%double(_dat31, bb56, ff13, hh13, gg13, 1, 5, 13);
%double(_dat31, bb56, ff14, hh14, gg14, 1, 6, 14);
%double(_dat31, bb56, ff15, hh15, gg15, 1, 7, 15);
%double(_dat31, bb56, ff16, hh16, gg16, 1, 8, 16);

%double(_dat31, bb56, ff17, hh17, gg17, 2, 10, 17);
%double(_dat31, bb56, ff18, hh18, gg18, 2, 3, 18);
%double(_dat31, bb56, ff19, hh19, gg19, 2, 4, 19);
%double(_dat31, bb56, ff20, hh20, gg20, 2, 5, 20);
%double(_dat31, bb56, ff21, hh21, gg21, 2, 6, 21);
%double(_dat31, bb56, ff22, hh22, gg22, 2, 7, 22);
%double(_dat31, bb56, ff23, hh23, gg23, 2, 8, 23);

%double(_dat31, bb56, ff24, hh24, gg24, 3, 11, 24);
%double(_dat31, bb56, ff25, hh25, gg25, 3, 4, 25);
%double(_dat31, bb56, ff26, hh26, gg26, 3, 5, 26);
%double(_dat31, bb56, ff27, hh27, gg27, 3, 6, 27);
%double(_dat31, bb56, ff28, hh28, gg28, 3, 7, 28);
%double(_dat31, bb56, ff29, hh29, gg29, 3, 8, 29);

%double(_dat31, bb56, ff30, hh30, gg30, 4, 12, 30);
%double(_dat31, bb56, ff31, hh31, gg31, 4, 5, 31);
%double(_dat31, bb56, ff32, hh32, gg32, 4, 6, 32);
%double(_dat31, bb56, ff33, hh33, gg33, 4, 7, 33);
%double(_dat31, bb56, ff34, hh34, gg34, 4, 8, 34);

%double(_dat31, bb56, ff35, hh35, gg35, 5, 13, 35);
%double(_dat31, bb56, ff36, hh36, gg36, 5, 6, 36);
%double(_dat31, bb56, ff37, hh37, gg37, 5, 7, 37);
%double(_dat31, bb56, ff38, hh38, gg38, 5, 8, 38);

%double(_dat31, bb56, ff39, hh39, gg39, 6, 14, 39);
%double(_dat31, bb56, ff40, hh40, gg40, 6, 7, 40);
%double(_dat31, bb56, ff41, hh41, gg41, 6, 8, 41);

%double(_dat31, bb56, ff42, hh42, gg42, 7, 15, 42);
%double(_dat31, bb56, ff43, hh43, gg43, 7, 8, 43);

%double(_dat31, bb56, ff44, hh44, gg44, 8, 16, 44);

data temp40;
  format method $3.0;
  set gg01 gg02 gg03 gg04 gg05 gg06 gg07 gg08 gg09 gg10
  gg11 gg12 gg13 gg14 gg15 gg16 gg17 gg18 gg19 gg20
  gg21 gg22 gg23 gg24 gg25 gg26 gg27 gg28 gg29 gg30
  gg31 gg32 gg33 gg34 gg35 gg36 gg37 gg38 gg39 gg40
  gg41 gg42 gg43 gg44 gg45 gg46 gg47;
*  set gg46;
run;

data temp50;
  format method $3.0;
  set ff01 ff02 ff03 ff04 ff05 ff06 ff07 ff08 ff09 ff10
  ff11 ff12 ff13 ff14 ff15 ff16 ff17 ff18 ff19 ff20
  ff21 ff22 ff23 ff24 ff25 ff26 ff27 ff28 ff29 ff30
  ff31 ff32 ff33 ff34 ff35 ff36 ff37 ff38 ff39 ff40
  ff41 ff42 ff43 ff44 ff45 ff46 ff47;
*  set ff46;
run;

data temp40;
  set temp40;
  where converge = 'Yes' and AIC ne .;
  i = 1;
run;

proc sort data = temp40;
  by replicate method aic;
run;

proc sort data = temp50;
  by replicate id;
run;

data temp41 (drop = i) temp42 (drop = i);
  set temp40;
  by replicate method aic;
  if first.method and method = 'FP' then output temp41;
  if first.method and method = 'RCS' then output temp42;
run;

proc sort data = temp40;
  by replicate i aic;
run;

data temp46;
  set temp40;
  by replicate i aic;
  if first.i;
  drop i; 
run;

data temp43;
  merge temp46 (in = a) temp50;
  by replicate id;
  if a;
  keep replicate id p1 p2 nknot parameter estimate StandardError;
run;

proc mixed data=_dat noclprint method=ml cl covtest;
  model &med = &exp &Cov / s cl covb; 
  ods output solutionf=beta;
  ods output covparms=cov;
  ods output covb = covb;
  by replicate;
run;

data tt53;
  set beta;
  if replicate = 0;
  i = 1;
  np2 = _n_;
run;

data tt54;
  set tt53;
  by i;
  if last.i;
  keep np2;
run; 

data _null_;
  set tt54;
  call symput('NP2', np2);
run;

data beta2; 
  length parameter $8; 
  set beta END=lastobs;
  if replicate = 0;
  format row 3.0 linpred $255.;
  row = _n_-1;
  parameter = "b" || left(row);
  linpred =  "+" || cats(parameter) || "*" || effect;
  if _N_=1 then call symput('linpredMED',parameter);
  else call symput('linpredMED',trim(resolve('&linpredMED'))||linpred);
run; 

data beta1; 
  length parameter $8; 
  set beta END=lastobs;
  format row 3.0;
  row = (_n_ - replicate * &np2)-1;
  parameter = "b" || left(row);
run; 

%put &linpredMED;

data cov1; 
  length parameter $8; 
  set cov;
  format row 3.0;
  row = _n_-1;
  if covparm="Residual" then do;
    parameter="Sigma";
    estimate = sqrt(estimate);
  end;
run; 

data medpar1;
  set beta1 cov1;
run;

proc sort data = medpar1;
  by replicate;
run;

proc transpose data = medpar1 out = medpar2 (drop = _name_);
  id parameter;
  var estimate;
  by replicate;
run;

proc transpose data = temp43 out = temp60 (drop = _name_ id);
  by replicate id;
  id parameter;
  var estimate;
  where substr(parameter, 1, 1) = 'a';
run;

proc transpose data = temp43 out = temp61 (drop = _name_ id);
  by replicate id;
  id parameter;
  var estimate;
  where substr(parameter, 1, 1) in ('g', 'l');
run;

data final;
  replicate = -1;
run;

%if &reference = 0 %then %do;

  proc sort data = _dat out = _dat00;
    by replicate &Cov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 0;
    output out = _dat01(drop = _type_ _freq_) n(&exp) = nsub;
    by replicate &Cov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 0;
    output out = _dat03(drop = _type_ _freq_) n(&exp) = total;
    by replicate;
  run;

%end;

%if &reference = 1 %then %do;

  proc sort data = _dat out = _dat00;
    by replicate &Cov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 1;
    output out = _dat01(drop = _type_ _freq_) n(&exp) = nsub;
    by replicate &Cov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    where &exp = 1;
    output out = _dat03(drop = _type_ _freq_) n(&exp) = total;
    by replicate;
  run;

%end;

%if &reference = 2 %then %do;

  proc sort data = _dat out = _dat00;
    by replicate &Cov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    output out = _dat01(drop = _type_ _freq_) n(&exp) = nsub;
    by replicate &Cov;
  run;

  proc means data = _dat00 noprint;
    var &exp;
    output out = _dat03(drop = _type_ _freq_) n(&exp) = total;
    by replicate;
  run;

%end;

%if &reference = 3 %then %do;

  proc sort data = ref out = _dat00;
    by &Cov;
  run;

  data _dat01;
    set _dat00;
    do i = 1 to (&nboot + 1);
    replicate = i - 1;
    output;
    end;
    drop i;
  run;

  proc sort data = _dat01;
    by replicate &Cov;
  run;

  proc means data = _dat01 noprint;
    var nsub;
    output out = _dat03(drop = _type_ _freq_) sum(nsub) = total;
    by replicate;
  run;

%end;

data _dat05;
  merge _dat01 _dat03;
  by replicate;
  percent = nsub/total * 100;
  drop nsub total;
  id = _n_;
run;

%med;

data final0;
  set final;
  if replicate ne (-1);
run;

data final1;
  merge final0 _dat05;
  by replicate id;
  S_T = S_T * percent/100;
  rmst = rmst * percent/100;
  hazard = hazard * percent/100;
run;

proc sort data = final1;
  by replicate ind;
run;

proc means data = final1 noprint;
  var s_t;
  by replicate ind;
  output out = final2 sum(s_t) = S_t sum(rmst) = rmst sum(hazard) = hazard;
run;

proc transpose data = final2 out = final3 (rename = (_1 = t_1 _2 = t_2 _3 = t_3 _4 = t_4));
  var s_t;
  id ind;
  by replicate;
run;

proc transpose data = final2 out = final4 (rename = (_1 = rm_1 _2 = rm_2 _3 = rm_3 _4 = rm_4));
  var rmst;
  id ind;
  by replicate;
run;

proc transpose data = final2 out = final9 (rename = (_1 = hazard_1 _2 = hazard_2 _3 = hazard_3 _4 = hazard_4));
  var hazard;
  id ind;
  by replicate;
run;

data final5;
  merge final3 final4 final9;
  by replicate;
  tes = (t_4 - t_1) * 100;
  ie0s = (t_2 - t_1) * 100;
  de1s = (t_4 - t_2) * 100;
  p0s = ie0s/tes * 100;
  ie1s = (t_4 - t_3) * 100;
  de0s = (t_3 - t_1) * 100;
  p1s = ie1s/tes * 100;

  term = rm_4 - rm_1;
  ie0rm = rm_2 - rm_1;
  de1rm = rm_4 - rm_2;
  p0rm = ie0rm/term * 100;
  ie1rm = rm_4 - rm_3;
  de0rm = rm_3 - rm_1;
  p1rm = ie1rm/term * 100;

  tehh = hazard_4 / hazard_1;
  ie0hh = hazard_2 / hazard_1;
  de1hh = hazard_4 / hazard_2;
  p0hh = log(ie0hh)/log(tehh) * 100;
  ie1hh = hazard_4 / hazard_3;
  de0hh = hazard_3 / hazard_1;
  p1hh = log(ie1hh)/log(tehh) * 100;
run;

data ff0;
  set final5;
  if replicate = 0;
/*  format ies des tes ps ierm derm term prm iehh dehh tehh phh ieh deh teh ph iedh dedh tedh pdh 8.2;*/
run;

data ff1;
  set final5;
  if replicate ne 0;
run;

proc univariate data = ff1 noprint;
  var tes ie0s de1s p0s ie1s de0s p1s term ie0rm de1rm p0rm ie1rm de0rm p1rm tehh ie0hh de1hh p0hh ie1hh de0hh p1hh;
  output out=ff2 pctlpts=2.5 97.5  pctlpre = tes ie0s de1s p0s ie1s de0s p1s term ie0rm de1rm p0rm ie1rm de0rm p1rm tehh ie0hh de1hh p0hh ie1hh de0hh p1hh pctlname = _025 _975;
run;

data ff30 (keep = tes ie0s de1s p0s ie1s de0s p1s tes_025 ie0s_025 de1s_025 p0s_025 ie1s_025 de0s_025 p1s_025 tes_975 ie0s_975 de1s_975 p0s_975 ie1s_975 de0s_975 p1s_975
rename = (tes = te ie0s = ie0 de1s = de1 p0s = p0 ie1s = ie1 de0s = de0 p1s = p1
tes_025 = te_025 ie0s_025 = ie0_025 de1s_025 = de1_025 p0s_025 = p0_025 ie1s_025 = ie1_025 de0s_025 = de0_025 p1s_025 = p1_025 
tes_975 = te_975 ie0s_975 = ie0_975 de1s_975 = de1_975 p0s_975 = p0_975 ie1s_975 = ie1_975 de0s_975 = de0_975 p1s_975 = p1_975)) 
ff32 (keep = term ie0rm de1rm p0rm ie1rm de0rm p1rm term_025 ie0rm_025 de1rm_025 p0rm_025 ie1rm_025 de0rm_025 p1rm_025 term_975 ie0rm_975 de1rm_975 p0rm_975 ie1rm_975 de0rm_975 p1rm_975
rename = (term = te ie0rm = ie0 de1rm = de1 p0rm = p0 ie1rm = ie1 de0rm = de0 p1rm = p1
term_025 = te_025 ie0rm_025 = ie0_025 de1rm_025 = de1_025 p0rm_025 = p0_025 ie1rm_025 = ie1_025 de0rm_025 = de0_025 p1rm_025 = p1_025 
term_975 = te_975 ie0rm_975 = ie0_975 de1rm_975 = de1_975 p0rm_975 = p0_975 ie1rm_975 = ie1_975 de0rm_975 = de0_975 p1rm_975 = p1_975))
ff33(keep = tehh ie0hh de1hh p0hh ie1hh de0hh p1hh tehh_025 ie0hh_025 de1hh_025 p0hh_025 ie1hh_025 de0hh_025 p1hh_025 tehh_975 ie0hh_975 de1hh_975 p0hh_975 ie1hh_975 de0hh_975 p1hh_975
rename = (tehh = te ie0hh = ie0 de1hh = de1 p0hh = p0 ie1hh = ie1 de0hh = de0 p1hh = p1
tehh_025 = te_025 ie0hh_025 = ie0_025 de1hh_025 = de1_025 p0hh_025 = p0_025 ie1hh_025 = ie1_025 de0hh_025 = de0_025 p1hh_025 = p1_025 
tehh_975 = te_975 ie0hh_975 = ie0_975 de1hh_975 = de1_975 p0hh_975 = p0_975 ie1hh_975 = ie1_975 de0hh_975 = de0_975 p1hh_975 = p1_975));
  merge ff0 ff2;
run;

data ff31;
  retain te te_025 te_975 ie0 ie0_025 ie0_975 de1 de1_025 de1_975 p0 p0_025 p0_975 ie1 ie1_025 ie1_975 de0 de0_025 de0_975 p1 p1_025 p1_975;
  set ff30;
run;

data temp14;
  set temp13;
  if replicate = 0;
  t_access = &t_access;
  max_time = exp(k_max);
  if t_access ne 999999 and t_access gt max_time then sta = 1;
  else sta = 0;
run;

data _null_;
  set temp14;
  call symput('stata', sta);
run;

%if &stata = 1 %then %do;
  proc iml;
    print,"WARNING: THE USER REQUIRED MEDIATION EFFECTS ASSESSMENT AFTER THE LAST OBSERVED EVENT TIME (NOT SUGGESTED, AND THE REASON IS THAT BASELINE PARAMETERIC DISTRIUBTION APPROXIMATION WHICH IS REQUIRED FOR THE MEDIATION EFFECTS ESTIMATION MAY NOT BE ACCURATE, SEE REFERENCE FOR DETAILS)";
  quit;
%end;

proc datasets lib = work;
  save ff31 ff32 ff33;
run; 

data &out;
  retain scale;
  set ff31 (in = a) ff32 (in = b) ff33 (in = c);
  format scale $20.;
  if a then scale = 'Survival Scale';
  if b then scale = 'RMST Scale';
  if c then scale = 'Hazard Scale';
run;

ods listing;

%if &printall=T %then %do;
  /*dm  'log;clear;out;clear;';*/
  proc print data = &out noobs;
  run;

  proc datasets lib = work;
    delete ff31 ff32 ff33;
  run; 

%end;  *-end printall=T option;

%end;

%mend;

*******************************************************************************************************
Output Data Set Variable Explanation:

Survival Scale: 
TE(t*) = mu(t*: 1, M(1)) - mu(t*: 0, M(0))
IE0(t*) = mu(t*: 0, M(1)) - mu(t*: 0, M(0))
DE1(t*) = mu(t*: 1, M(1)) - mu(t*: 0, M(1))
P0(t*) = IE0(t*) / TE(t*)
IE1(t*) = mu(t*: 1, M(1)) - mu(t*: 1, M(0))
DE0(t*) = mu(t*: 1, M(0)) - mu(t*: 0, M(0))
P1(t*) = IE1(t*) / TE(t*)

RMST Scale: 
TE(t*) = rm(t*: 1, M(1)) - rm(t*: 0, M(0))
IE0(t*) = rm(t*: 0, M(1)) - rm(t*: 0, M(0))
DE1(t*) = rm(t*: 1, M(1)) - rm(t*: 0, M(1))
P0(t*) = IE0(t*) / TE(t*)
IE1(t*) = rm(t*: 1, M(1)) - rm(t*: 1, M(0))
DE0(t*) = rm(t*: 1, M(0)) - rm(t*: 0, M(0))
P1(t*) = IE1(t*) / TE(t*)

Hazard Scale: 
TE(t*) = h(t*: 1, M(1)) / h(t*: 0, M(0))
IE0(t*) = h(t*: 0, M(1)) / h(t*: 0, M(0))
DE1(t*) = h(t*: 1, M(1)) / h(t*: 0, M(1))
P0(t*) = log (IE0(t*)) / log (TE(t*))
IE1(t*) = h(t*: 1, M(1)) / h(t*: 1, M(0))
DE0(t*) = h(t*: 1, M(0)) / h(t*: 0, M(0))
P1(t*) = log (IE1(t*)) / log (TE(t*))
*******************************************************************************************************;

proc import datafile = 'C:\Users\Wei Wang\Box Sync\Survival Mediation Analysis\2-programs\SAS Macro Upload\example.xls' out = ttemp1 replace;
run; 

libname te 'C:\Users\Wei Wang\Desktop\SM_Results\Temp4';

%SURVPHMED (data             = TTEMP1, 
            SURVtime         = SURVTIME,
            SURVevent        = SURVevent,
            MED              = MED, 
            Cov              = COV1 COV2 COV3, 
            EXP              = EXP,
            Nboot            = 5,
			T_ACCESS         = 999999,
			REFERENCE        = 1,
            Nlmopts          = ,
            Printall         = T,
            Out              = te.out01
           );

proc import datafile = 'C:\Users\Wei Wang\Box Sync\Survival Mediation Analysis\2-programs\SAS Macro Upload\example.xls' out = ttemp1 replace;
run; 

data ref;
input COV1 COV2 COV3 nsub;
datalines;
1 1 0 10
1 0 0 10
;
run;


		   %SURVPHMED (data             = TTEMP1, 
            SURVtime         = SURVTIME,
            SURVevent        = SURVevent,
            MED              = MED, 
            Cov              = COV1 COV2 COV3, 
            EXP              = EXP,
            Nboot            = 5,
			T_ACCESS         = 999999,
			REFERENCE        = 3,
            Nlmopts          = ,
            Printall         = T,
            Out              = te.out02
           );
