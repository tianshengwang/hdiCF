/************************************************************************************************************************/
/* Program: /local/projects/medicare/sglt_cvs/programs/4c_hc_covariates.sas                                             */
/* Purpose:  select variables for hdiCF 	                                                                                */
/* Author: Virginia Pate                                                                                                */
/************************************************************************************************************************/

options sasautos=(SASAUTOS "/local/projects/medicare/sglt_cvs/programs/macros");
%setup(full, 4_covariates, saveLog=N);

%LET drug1=SGLT;
%LET drug2=GLP;

%setup(full, select_variables_for_causal_forest, saveLog=N);

/*%setup(full, select_variables_for_causal_forest, saveLog=N)*/

%macro select_vars(dxGroup=3, atcGroup=4, max_num=100, min_prev=1);

	*Select top &max_num prevalent variables from each dimension;
	%macro select_vars_sub(dimension=);
		%GLOBAL list_&dimension N_&dimension;
	
		proc sql noprint;
			select distinct "&dimension" || "_" || strip(code) into :list_&dimension separated by ' '
				from out.ref_hdCov_&drug1.v&drug2
				where lowcase(dimension)="&dimension" and prev_order<=&max_num and prevalence>=&min_prev/100;
			%LET N_&dimension = &SqlObs;
		quit;
	%mend;
	
	%select_vars_sub(dimension=dx&dxGroup._inpt)
	%select_vars_sub(dimension=dx&dxGroup._outpt)
	%select_vars_sub(dimension=cpt5_inpt)
	%select_vars_sub(dimension=cpt5_outpt)
	%select_vars_sub(dimension=atc&atcGroup._outpt)


	*Keep only those variables on the cohort dataset;
	data hdCov2_&drug1.v&drug2._prev&min_prev._num&max_num;
		set out.hdCov_&drug1.v&drug2;
		keep bene_id indexdate
			  &&&list_dx&dxGroup._inpt
			  &&&list_dx&dxGroup._outpt
			  &list_cpt5_inpt
			  &list_cpt5_outpt
			  &&&list_atc&atcGroup._outpt
		;
	run;
	
	*Add demographic variables to cohort datasets;
 	proc sql;  
 		create table hdcov1_&drug1.v&drug2._prev&min_prev._num&max_num as select a.*, b.*,		 
 				case when mdy(month(b.dob),day(b.dob),2000)<mdy(month(a.indexdate),day(a.indexdate),2000) 
 	                   then year(a.indexdate) - year(b.dob) 
 	              else year(a.indexdate) - year(b.dob) + 1 end as age 
 			from hdCov2_&drug1.v&drug2._prev&min_prev._num&max_num as a 
				left join temp.newusers_&drug1.v&drug2 as b 
 					on a.bene_id=b.bene_id and a.indexdate=b.indexdate; 
 	quit; 
 	
	*Add outcomes and output final permanent dataset with all variables needed for the HDiCF analysis;
 	proc sql;
		create table out.hdicf_&drug1.v&drug2._p&min_prev._n&max_num._i&dxGroup._A&atcgroup. as
			select a.*, b.* 
		from temp.outcomes_sgltvglp as a 
			left join hdCov1_&drug1.v&drug2._prev&min_prev._num&max_num as b
		on a.bene_ID=b.bene_ID and a.indexdate=b.indexdate;
   quit;
 	
%mend;
options mprint;

/*primary analysis & sensitivity anlaysis 1 & 2*/
%select_vars(dxGroup=3, atcgroup=3,  min_prev=1, max_num=200);
/*Sensitivity analysis 3: First, we utilized the 100 most prevalent codes (n=100) in each dimension.*/
%select_vars(dxGroup=3, atcgroup=3,  min_prev=1, max_num=100) 
/*Sensitivity analysis 4: Second, we applied a cutoff c=0.02 for these codes.*/ 
%select_vars(dxGroup=3, atcgroup=3,  min_prev=2, max_num=200);
/*Sensitivity analysis 5: Second, we applied a cutoff c=0.05 for these codes.*/ 
%select_vars(dxGroup=3, atcgroup=3,  min_prev=5, max_num=200);
/*Sensitivity analysis 6: Third, we employed a 4-digit convention for the ICD-10 codes.*/
%select_vars(dxGroup=4, atcgroup=3,  min_prev=1, max_num=200); 
/*Sensitivity analysis 7: Fourth, we used the 4th level for the ATC codes. */
%select_vars(dxGroup=3, atcgroup=4, min_prev=1, max_num=200);
/*Sensitivity analysis 8: Fourth, we used the 4th level for the ATC codes. */
%select_vars(dxGroup=3, atcgroup=4, min_prev=2, max_num=200);
/*Sensitivity analysis 9: Fourth, we used the 4th level for the ATC codes. */
%select_vars(dxGroup=3, atcgroup=4, min_prev=5, max_num=200);
/*Sensitivity analysis 10: Fourth, we used the 4th level for the ATC codes. */
%select_vars(dxGroup=3, atcgroup=4, min_prev=5, max_num=100);
/*Sensitivity analysis 11: Fourth, we used the 4th level for the ATC codes. */
%select_vars(dxGroup=3, atcgroup=4, min_prev=5, max_num=50);


proc print data=out.hdicf_sgltvglp_p1_n200_i3_a4 (obs=1);run;
proc print data=out.hdicf_sgltvglp_p1_n200_i3_a4 (obs=1);run;
proc print data=out.hdicf_sgltvglp_p1_n200_i3_a5 (obs=1);run;
