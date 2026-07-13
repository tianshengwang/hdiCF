/************************************************************************************************************************/
/* hdiCF compatibility bundle: the %select_vars variable-selection logic                                                */
/*                                                                                                                      */
/* This exercises the variable-selection logic from SAS/select_variables_for_hdiCF.sas                                  */
/* (Wang T, Pate V, et al., "High-dimensional Iterative Causal Forest", Am J Epidemiol 2024). For each claims           */
/* dimension it reads the recurrence reference table, keeps the codes whose prevalence clears the c cutoff and whose    */
/* prevalence rank is within the top max_num, builds the "<dimension>_<code>" variable-name list with a PROC SQL        */
/* INTO ... SEPARATED BY, then keeps exactly those columns on the wide covariate dataset -- the feature set that gets   */
/* handed to the causal forest.                                                                                         */
/*                                                                                                                      */
/* The original reads the permanent ref_hdCov / hdCov datasets produced by generate_HD_variables.sas from external SAS  */
/* libraries; here both are built inline (ref_hdCov matches the ref_<dimension> tables that %hdcov emits), so the       */
/* selection and keep steps run end to end with no external dependencies.                                               */
/************************************************************************************************************************/

/* ---- recurrence reference: one row per (dimension, code) with prevalence rank + prevalence, as %hdcov emits ---- */
data ref_hdCov;
  length dimension $10 code $7;
  input dimension $ code $ prev_order prevalence;
  datalines;
dx3_inpt I50 1 0.42
dx3_inpt I25 2 0.31
dx3_inpt E11 3 0.08
dx3_inpt N18 4 0.005
dx3_outpt I50 1 0.55
dx3_outpt E11 2 0.22
dx3_outpt I48 3 0.011
cpt5_outpt 85610 1 0.18
cpt5_outpt 80053 2 0.09
atc4_outpt C03C 1 0.29
;
run;

/* ---- wide covariate dataset: bene_id/indexdate + one ordinal column per candidate variable ---- */
data hdCov;
  input bene_id $ indexdate :date9.
        dx3_inpt_I50 dx3_inpt_I25 dx3_inpt_E11 dx3_inpt_N18
        dx3_outpt_I50 dx3_outpt_E11 dx3_outpt_I48
        cpt5_outpt_85610 cpt5_outpt_80053
        atc4_outpt_C03C;
  format indexdate date9.;
  datalines;
B001 01JAN2016 3 2 1 0 3 1 0 2 0 3
B002 01JAN2016 1 0 2 1 2 3 1 0 1 2
B003 15FEB2016 2 3 0 0 1 2 0 3 2 1
;
run;

/*======================================================================================================================*/
/* %select_vars_sub -- verbatim from SAS/select_variables_for_hdiCF.sas: build the top-max_num, prevalence>=c            */
/* variable-name list for one dimension                                                                                 */
/*======================================================================================================================*/
%macro select_vars_sub(dimension=, max_num=100, min_prev=1);
	%GLOBAL list_&dimension N_&dimension;

	proc sql noprint;
		select distinct "&dimension" || "_" || strip(code) into :list_&dimension separated by ' '
			from ref_hdCov
			where lowcase(dimension)="&dimension" and prev_order<=&max_num and prevalence>=&min_prev/100;
		%LET N_&dimension = &SqlObs;
	quit;
%mend;

/* Build the selected-variable list per dimension (c=1% cutoff, top 100 by prevalence rank) */
%select_vars_sub(dimension=dx3_inpt)
%select_vars_sub(dimension=dx3_outpt)
%select_vars_sub(dimension=cpt5_outpt)
%select_vars_sub(dimension=atc4_outpt)

%put NOTE: selected dx3_inpt   = &list_dx3_inpt;
%put NOTE: selected dx3_outpt  = &list_dx3_outpt;
%put NOTE: selected cpt5_outpt = &list_cpt5_outpt;
%put NOTE: selected atc4_outpt = &list_atc4_outpt;

/* Keep only the selected variables on the wide covariate dataset -- the feature matrix for the causal forest.       */
/* Codes below the c=1% cutoff (dx3_inpt_N18 at 0.5%) are dropped exactly as the manuscript's selection prescribes.  */
data hdCov_selected;
	set hdCov;
	keep bene_id indexdate
	     &list_dx3_inpt
	     &list_dx3_outpt
	     &list_cpt5_outpt
	     &list_atc4_outpt
	;
run;

proc print data=hdCov_selected;
  title "hdiCF select_vars: covariate matrix after prevalence-cutoff variable selection";
run;

proc contents data=hdCov_selected varnum;
  title "hdiCF select_vars: selected feature columns kept for the causal forest";
run;
title;
