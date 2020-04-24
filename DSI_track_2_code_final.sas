*DSI competition Track 2;
*import dataset;
PROC IMPORT out=covid
	datafile='/home/u45096072/DSI competition database - Final.csv'
	DBMS=CSV replace;
	GETNAMES=YES;
RUN;

proc contents data=covid;

*Part1: model for incidence rate;

*Step 1: look at distribution of all variables and potential transformation needed.
*histograms;

proc sgscatter data=covid;
	matrix incidence_rate black--unemployment / diagonal=(histogram);
run;

proc sgscatter data=covid;
	matrix incidence_rate older_adults--hospital_land / diagonal=(histogram);;
run;

proc sgscatter data=covid;
	matrix incidence_rate hospital_10k--food_insecurity / diagonal=(histogram);;
run;

proc sgscatter data=covid;
	matrix incidence_rate google_trend--life_exp / diagonal=(histogram);;
run;

*Step 2: log transformation;
data covid2;
	set covid;
	log_incidence=log10(incidence_rate);
	log_black=log10(black);
	log_education=log10(education);
	log_income=log10(income);
	log_older_adults=log10(older_adults);
	log_uninsured=log10(uninsured);
	log_hospital1=log10(total_hospitals+1);
	log_ranking=log10(national_rankings+1);
	log_hospital2=log10(hospital_land+1);
	log_hospital3=log10(hospital_10k+1);
	log_food1=log10(food_store);
	log_food2=log10(food_store_land);
	log_food3=log10(food_retail_1k);
	log_google=log10(google_trend+1);
run;

*Check scatterplots after log transformation;
proc sgscatter data=covid2;
	matrix log_incidence poverty unemployment food_insecurity / diagonal=(histogram);
run;

proc sgscatter data=covid2;
	matrix log_incidence diabetes_prev hypertension_prev land_area life_exp / diagonal=(histogram);;
run;

proc sgscatter data=covid2;
	matrix log_incidence log_black log_education log_income log_older_adults / diagonal=(histogram);;
run;

proc sgscatter data=covid2;
	matrix log_incidence log_uninsured log_hospital1 log_ranking log_hospital2 log_hospital3 / diagonal=(histogram);;
run;

proc sgscatter data=covid2;
	matrix log_incidence log_food1 log_food2 log_food3 log_google / diagonal=(histogram);;
run;

*step 3: fit simple regression on all remaining potential predictors;
*simple regression;
proc reg data=covid2;
	model log_incidence = log_black;
	model log_incidence = log_education;
	model log_incidence = poverty;
	model log_incidence= log_income;
	model log_incidence = unemployment;
	model log_incidence = log_older_adults;
	model log_incidence= log_uninsured;
	model log_incidence = log_hospital1;
	model log_incidence = log_hospital3;
run;

proc reg data=covid2;
	model log_incidence= food_insecurity;
	model log_incidence = log_food1;
	model log_incidence = log_food2;
	model log_incidence = log_food3;
	model log_incidence = log_google;
	model log_incidence = Diabetes_prev;
	model log_incidence = Hypertension_prev;
	model log_incidence = land_area;
	model log_incidence = life_exp;
run;

*step 4: For all potential variables, check if there is any multicollinearity problem;
*grouped scatterplots and vif check;
proc corr data=covid2 plots(maxpoints=100000000)=matrix(nvar=all);
	var log_incidence log_black log_education log_income unemployment;
run;

proc corr data=covid2 plots(maxpoints=100000000)=matrix(nvar=all);
	var log_incidence log_hospital3 land_area life_exp;
run;

proc corr data=covid2 plots(maxpoints=100000000)=matrix(nvar=all);
	var log_incidence log_food1 log_food2 food_insecurity land_area;
run;

proc corr data=covid2 plots(maxpoints=100000000)=matrix(nvar=all);
	var log_incidence diabetes_prev hypertension_prev log_older_adults life_exp;
run;

proc reg data=covid2;
	model log_incidence = log_black log_education log_income unemployment log_older_adults
						  log_hospital3 food_insecurity log_food1 log_food2 
						  diabetes_prev hypertension_prev land_area life_exp / vif r;
run; quit;

*step 5: Fit lasso regression model with the remaining variables;
*standardize data;
PROC STANDARD DATA=covid2 MEAN=0 STD=1 OUT=covid_standardized;
	VAR log_black log_education log_income unemployment log_hospital3 
		food_insecurity diabetes_prev hypertension_prev life_exp;
RUN;

*cross validation;
PROC SURVEYSELECT DATA=covid_standardized OUT=TRAINTEST SEED=123
	SAMPRATE=0.7 METHOD=SRS OUTALL;
RUN;

PROC GLMSELECT DATA=TRAINTEST PLOTS=ALL SEED=123;
	PARTITION ROLE=SELECTED(TRAIN='1' TEST='0');
	MODEL log_incidence = log_black log_education log_income unemployment log_hospital3 
						  food_insecurity diabetes_prev hypertension_prev life_exp 
						  / SELECTION=LAR (CHOOSE=CV STOP=NONE) CVMETHOD=RANDOM(10);
RUN;

*lasso regression;
proc glmselect data=covid_standardized plots=all;
	model log_incidence = log_black log_education log_income unemployment log_hospital3 
						  food_insecurity diabetes_prev hypertension_prev life_exp /selection=lasso(stop=none choose=bic);
run; 

*The remaining variables after lasso regression;
proc reg data=covid2;
	model log_incidence = log_black log_income log_hospital3 food_insecurity life_exp / vif r;
run; quit;

*Examine non-siginificant variables' effects on the model;
*Remove log_income;
proc reg data=covid2;
	model log_incidence = log_black log_hospital3 food_insecurity life_exp / vif r;
run; quit;
*Remove log_hospital3;
proc reg data=covid2;
	model log_incidence = log_black log_income food_insecurity life_exp / vif r;
run; quit;
*Remove both;
proc reg data=covid2;
	model log_incidence = log_black food_insecurity life_exp / vif r;
run; quit;

*Potential confounder based on knowledge;
proc reg data=covid2;
	model log_incidence = food_insecurity;
	output out=regout p=yhat r=resid;
run; quit;

proc reg data=covid2;
	model log_incidence = unemployment food_insecurity / vif r;
	output out=regout p=yhat r=resid;
run; quit;

*Unemployment is a confounder for food_insecurity and log_incidence;
*We would include it in the model;

*Final model;
proc reg data=covid2;
	model log_incidence = log_black unemployment life_exp food_insecurity;
	output out=regout p=yhat r=resid;
run; quit;

* Residual: normality check;
proc univariate data=regout normal;
 var resid;
 histogram resid;
 qqplot resid;
run;

*Part2: model for fatality rate;

*Step 1: look at distribution of all variables and potential transformation needed.
*histograms;

proc sgscatter data=covid;
	matrix fatality_rate black--unemployment / diagonal=(histogram);
run;

proc sgscatter data=covid;
	matrix fatality_rate older_adults--hospital_land / diagonal=(histogram);
run;

proc sgscatter data=covid;
	matrix fatality_rate hospital_10k--food_insecurity / diagonal=(histogram);
run;

proc sgscatter data=covid;
	matrix fatality_rate google_trend--life_exp / diagonal=(histogram);
run;

*Step 2: Transformation of outcome;
data covid3;
	set covid;
	sqrt_fatality=sqrt(fatality_rate);
	c_fatality=(fatality_rate+10);
	log_fatality=log10(fatality_rate+1);
run;

proc sgscatter data=covid3;
	matrix fatality_rate sqrt_fatality c_fatality log_fatality / diagonal=(histogram);
run;

*Step 2: log transformation of potential predictors;
data covid3;
	set covid3;
	log_black=log10(black);
	log_education=log10(education);
	log_income=log10(income);
	log_older_adults=log10(older_adults);
	log_uninsured=log10(uninsured);
	log_hospital1=log10(total_hospitals+1);
	log_ranking=log10(national_rankings+1);
	log_hospital2=log10(hospital_land+1);
	log_hospital3=log10(hospital_10k+1);
	log_food1=log10(food_store);
	log_food2=log10(food_store_land);
	log_food3=log10(food_retail_1k);
	log_google=log10(google_trend+1);
run;

*step 3: fit simple regression on all potential predictors;
*simple regression;
proc reg data=covid3;
	model sqrt_fatality = log_black;
	model sqrt_fatality = log_education;
	model sqrt_fatality = poverty;
	model sqrt_fatality = log_income;
	model sqrt_fatality = unemployment;
	model sqrt_fatality = log_older_adults;
	model sqrt_fatality = log_uninsured;
	model sqrt_fatality = log_ranking;
	model sqrt_fatality = log_hospital1;
	model sqrt_fatality = log_hospital2;
	model sqrt_fatality = log_hospital3;
run;

proc reg data=covid3;
	model sqrt_fatality = food_insecurity;
	model sqrt_fatality = log_food1;
	model sqrt_fatality = log_food2;
	model sqrt_fatality = log_food3;
	model sqrt_fatality = log_google;
	model sqrt_fatality = Diabetes_prev;
	model sqrt_fatality = Hypertension_prev;
	model sqrt_fatality = land_area;
	model sqrt_fatality = life_exp;
run;

*step 4: For all potential variables, check if there is any multicollinearity problem;
*scatterplots and vif check;
proc corr data=covid3 plots(maxpoints=100000000)=matrix(nvar=all);
	var sqrt_fatality log_black log_education log_income unemployment log_older_adults;
run;

proc corr data=covid3 plots(maxpoints=100000000)=matrix(nvar=all);
	var sqrt_fatality log_uninsured log_ranking log_hospital1 log_food1 log_food2;
run;

proc corr data=covid3 plots(maxpoints=100000000)=matrix(nvar=all);
	var sqrt_fatality diabetes_prev land_area;
run;

proc reg data=covid3;
	model sqrt_fatality = log_black log_education log_income unemployment log_older_adults
						  log_uninsured log_ranking log_hospital1 log_food1 log_food2
						  diabetes_prev land_area / vif r;
run; quit;

*step 5: Fit lasso multiple regression model with the remaining variables;
*standardize data;
PROC STANDARD DATA=covid3 MEAN=0 STD=1 OUT=covid_std2;
	VAR sqrt_fatality log_black log_education log_income unemployment log_older_adults
		log_uninsured log_ranking log_hospital1 log_food2 diabetes_prev land_area;
RUN;

*cross validation;
PROC SURVEYSELECT DATA=covid_std2 OUT=TRAINTEST2 SEED=123
	SAMPRATE=0.7 METHOD=SRS OUTALL;
RUN;

PROC GLMSELECT DATA=TRAINTEST2 PLOTS=ALL SEED=123;
	PARTITION ROLE=SELECTED(TRAIN='1' TEST='0');
	MODEL sqrt_fatality=log_black log_education log_income unemployment log_older_adults
					    log_uninsured log_ranking log_hospital1 log_food2 
					    diabetes_prev land_area / SELECTION=LAR (CHOOSE=CV STOP=NONE) CVMETHOD=RANDOM(10);
RUN;

*lasso regression;
proc glmselect data=covid_std2 plots=all;
	model sqrt_fatality = log_black log_education log_income unemployment log_older_adults
					      log_uninsured log_ranking log_hospital1 log_food2 
					      diabetes_prev land_area /selection=lasso(stop=none choose=bic);
run; 

*Final model;
proc reg data=covid3;
	model sqrt_fatality = log_black unemployment log_older_adults log_food2;
	output out=regout2 p=yhat r=resid;
run; quit;

* Residual: normality check;
proc univariate data=regout2 normal;
 var resid;
 histogram resid;
 qqplot resid;
run;

*None of the predictor is significant;

*Gamma for fatality was tried;
PROC GENMOD DATA = covid3;
	MODEL fatality_rate = log_black log_education poverty unemployment log_older_adults
					  	  log_ranking log_hospital1 log_hospital2 log_food1 log_food2
					  	  diabetes_prev hypertension_prev land_area life_exp / DIST = GAMMA LINK = LOG TYPE1;
*WEIGHT SURVEY_WEIGHT;
RUN;

*Nothing particularly interesting were found;


