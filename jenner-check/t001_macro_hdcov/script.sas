/************************************************************************************************************************/
/* hdiCF compatibility bundle: the %hdcov recurrence-assessment core                                                    */
/*                                                                                                                      */
/* This exercises the recurrence-assessment logic from %hdcov in SAS/generate_HD_variables.sas                          */
/* (Wang T, Pate V, et al., "High-dimensional Iterative Causal Forest", Am J Epidemiol 2024). It implements             */
/* Steps 2 & 3 of the hdPS covariate build: rank a claims dimension by prevalence, derive the median / Q3 / P90         */
/* count cutpoints, then assign the ordinal covariate (0/1/2/3) and the cov_once / cov_sporadic / cov_frequent          */
/* recurrence flags exactly as the manuscript describes.                                                                */
/*                                                                                                                      */
/* The original driver reads Medicare claims from external SAS libraries; here the same macro logic runs against a      */
/* small synthetic dx3 inpatient dimension (all_dx3_inpt) and a matching cohort, built inline so the recurrence         */
/* assessment runs end to end with no external dependencies.                                                            */
/************************************************************************************************************************/

/* ---- synthetic inputs matching the macro's contract ----                                                            */
/* cohort(bene_id indexdate), and all_<dimension>(bene_id indexdate <dimension> label count) as produced by %getdx     */
data cohort;
  input bene_id $ indexdate :date9.;
  format indexdate date9.;
  datalines;
B001 01JAN2016
B002 01JAN2016
B003 15FEB2016
B004 15FEB2016
B005 03MAR2016
B006 03MAR2016
B007 20APR2016
B008 20APR2016
B009 08MAY2016
B010 08MAY2016
;
run;

data all_dx3_inpt;
  length label $40;
  input bene_id $ indexdate :date9. dx3_inpt $ label $ count;
  format indexdate date9.;
  datalines;
B001 01JAN2016 I50 Heart_failure 5
B002 01JAN2016 I50 Heart_failure 1
B003 15FEB2016 I50 Heart_failure 3
B004 15FEB2016 I50 Heart_failure 2
B005 03MAR2016 I50 Heart_failure 4
B001 01JAN2016 I25 Chronic_ischemic_heart_dis 2
B002 01JAN2016 I25 Chronic_ischemic_heart_dis 1
B006 03MAR2016 I25 Chronic_ischemic_heart_dis 6
B007 20APR2016 I25 Chronic_ischemic_heart_dis 1
B003 15FEB2016 E11 Type2_diabetes 8
B004 15FEB2016 E11 Type2_diabetes 2
B005 03MAR2016 E11 Type2_diabetes 1
B008 20APR2016 E11 Type2_diabetes 3
B009 08MAY2016 E11 Type2_diabetes 1
B010 08MAY2016 N18 Chronic_kidney_disease 4
B009 08MAY2016 N18 Chronic_kidney_disease 2
B006 03MAR2016 N18 Chronic_kidney_disease 1
;
run;

/* N_cohort is the prevalence denominator the manuscript uses; the driver sets it before calling the macro */
proc sql noprint; select count(*) into :N_cohort from cohort; quit;

/*======================================================================================================================*/
/* %hdcov recurrence assessment -- the prevalence-ranking / percentile-cutpoint / ordinal + cov-flag logic taken        */
/* verbatim from %hdcov in SAS/generate_HD_variables.sas                                                                */
/*======================================================================================================================*/
%macro hdcov_recur(dimension=, max_num=200, min_prev=0.01);

	*prep for specific cohort;
	proc sql;
		create table &dimension as select distinct a.*
			from all_&dimension as a inner join cohort as b on a.bene_id=b.bene_id and a.indexdate=b.indexdate;
	quit;

	*get number of beneficiaries with at least one date with a claim for given dimension;
	proc freq data=&dimension order=freq;
		tables &dimension / out=&dimension._freq;
	run;

	*keep &MAX_NUM most prevalent;
	data &dimension._freq2;
		set &dimension._freq;
		pop_percent = count / &N_cohort;
		if pop_percent>=0.5 then do; pop_percent=1-pop_percent; count=&N_cohort-count; end;
		keep &dimension pop_percent count;
	run;

	proc sort data=&dimension._freq2; by descending pop_percent; run;

	data &dimension._freq3;
		set &dimension._freq2;
		prev_order=_N_;
		if _N_<=&max_num then output;
	run;

	*get median, 75th and 90th percentiles of counts among those with at least one claim;
	proc means data = &dimension         median        q3    p90; class &dimension; var count;
		 output out = &dimension._median median=median q3=q3 p90=p90;
	run;

	proc sql;
		create table &dimension._median2 as select a.* from &dimension._median as a
			  inner join &dimension._freq3   as b on a.&dimension = b.&dimension;
	quit;

/*--------------this is the same as the hdiCF manuscript described--------------*/
	proc sql;
		create table &dimension._ordinalVar as select distinct a.*, b.median, b.q3, b.p90,

			case when        a.count = 0 then 0
				 when        a.count = 1 then 1
		 		 when    1 < a.count <= b.q3 OR (b.q3=1 and a.count<=b.p90) then 2
				 when b.q3 < a.count then 3
				 else . end as ordinal_val,

			/*create variables for HDPS*/
			case when        a.count >= 1        then 1 else 0 end as cov_once,
			case when        a.count >= b.median then 1 else 0 end as cov_sporadic,
			case when        a.count >= b.q3     then 1 else 0 end as cov_frequent

		from &dimension as a inner join &dimension._median2 as b on a.&dimension = b.&dimension;
	quit;
/*-----------------------------------------------------------------------------*/

	*Pull reference file for the dimension (prevalence order + cutpoints);
	proc sql;
		create table ref_&dimension  as
			select distinct "&dimension" as dimension,
							a.&dimension  as code,
							a.prev_order, a.pop_percent as prevalence, a.count, b.median, b.q3, b.p90
			from &dimension._freq3 as a left join &dimension._median2 as b on a.&dimension=b.&dimension
			order by prev_order;
	quit;
%mend;

/* Run the published recurrence assessment on the synthetic dx3 inpatient dimension */
%hdcov_recur(dimension=dx3_inpt);

/* The prevalence-ordered reference: one row per code with its median/Q3/P90 count cutpoints */
proc print data=ref_dx3_inpt;
  title "hdiCF hdcov: recurrence reference (prevalence order + median/Q3/P90 cutpoints) by ICD-10 code";
run;

/* The per-claim ordinal covariate + cov_once/sporadic/frequent recurrence flags */
proc print data=dx3_inpt_ordinalVar;
  title "hdiCF hdcov: ordinal_val + cov_once/cov_sporadic/cov_frequent recurrence flags";
run;
title;
