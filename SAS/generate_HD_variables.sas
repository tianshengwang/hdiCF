
/************************************************************************************************************************/
/* Program: /local/projects/medicare/sglt_cvs/programs/4c_hc_covariates.sas                                             */
/* Purpose:  Get baseline characteristics                                                                               */
/* Author: Virginia Pate                                                                                                */
/* Updates:                                                                                                             */
/*    12/21/2022 - copied from adrd program                                                                             */
/*    3/6/2023   - corrected prevalence >50%, kept 200 variables, no max prevalence, added code to increase length of dx*/
/*    5/17/2023  - changed ordinal variable - if count = Q3, assign value of 2 rather than 3                          */
/*    5/31/2023  - corrected inpatient vs outpatient identification of procedures and diagnoses                       */
/*    4/6/2024  - Tian fixed the semantics issue (atc3 & atc4 previously named atc4 & atc5) and missing real atc4 data*/
/*    7/17/2024 - added 3 variables for use in determining HDPS covariates (CovX_once, CovX_sporadic, CovX_frequent)  */
/*    9/18/2024 - Tian revised "A code that appeared above the 75th percentile number of times would have a “true” value for all 3 recurrence variables. 
If any of the values were equal, the variable representing the higher cutpoint was dropped."*/
/**********************************************************************************************************************/
options source source2 msglevel=I mprint mcompilenote=all mautosource 
	sasautos=(SASAUTOS "/local/projects/medicare/sglt_cvs/programs/macros");


%LET drug1=SGLT;
%LET drug2=GLP;
%setup(full, 4c_hd_covariates_&drug1.v&drug2, saveLog=Y);


%LET startYear=2015;
%LET endYear=2019;
%LET bldays=365;

data cohort(keep=bene_id indexdate); 
	merge temp.inter_sgltvglp(in=a where=(indexdate>="15OCT2015"d+&bldays 
												AND filldate2 ne . and filldate2<="31DEC2017"d
												and excludeFlag_prevalentUser=0 
												and excludeFlag_sameDayInitiator=0 
                                    and excludeFlag_preFill2Initiator=0
												and excludeFlag_age=0 
												and excludeFlag_blenroll=0))
			temp.covariates_sgltvglp(in=b where=( bl_ckdstagef=0) keep=bene_id indexdate bl_ckdstagef);
	by bene_id indexdate;
	if a and b;
run;

/******************************************************************************************************/
/* STEP 1 & 2 (hdPS). Specify Data Sources AND Identify Candidate Empirical Covariates (granularity). */
/******************************************************************************************************/
/*DIAGNOSES*/
%macro getdx(startYr=&startYear, endYr=&endYear, dxGroup=3);
  proc sql;
      *COUNT = number of days beneficiary has at least one code within dimension;
      create table dx&dxGroup as
      select distinct bene_id, indexdate, inpt, dx&dxGroup, put(dx&dxGroup,$dxten.) as label, count(distinct from_dt) as count
		from ( 
			select distinct bene_id, indexdate, from_dt, dx&dxGroup, max(inpt=1) as inpt
			from (
				%DO yr=&startYr %TO &endYr;
		         select distinct a.bene_id, a.indexdate, b.from_dt, substr(b.dx&yr,1,&dxGroup) as dx&dxGroup, 
						case when upcase(b.source)='MEDPAR' then 1 else 0 end as inpt
		         from cohort as a 
		               inner join der.alldx10&yr as b 
								on a.bene_id=b.bene_id and a.indexdate-&bldays<=b.from_dt<=a.indexdate
	         %IF &yr<&endYr %THEN union corresponding; %END;  
				) 
			group by bene_id, indexdate, from_dt, dx&dxGroup
      	) 
		group by bene_id, indexdate, inpt, dx&dxGroup
		order by bene_id, indexdate, inpt, dx&dxGroup;
   quit;
%mend;
%getdx(dxGroup=3);
%getdx(dxGroup=4);


/* PROCEDURES */
%macro allproc(startYr=&startYear, endYr=&endYear);
   proc sql;
		create table cpt as
	      	select distinct bene_id, indexdate, inpt, cpt, put(cpt,$cpt.) as label, count(distinct proc_dt) as count
			from 
				(select distinct bene_id, indexdate, proc_dt, cpt, max(inpt=1) as inpt 
				from 
					(
						%DO yr=&startYr %TO &endYr;
				         select distinct a.bene_id, a.indexdate, b.thru_dt as proc_dt, b.hcpcs_cd as cpt,
								case when b.plcsrvc='21' then 1 else 0 end as inpt
				         from cohort as a 
				            inner join raw.bcarrier_line&yr as b
				               on a.bene_id = b.bene_id and a.indexdate-&bldays<=b.thru_dt<=a.indexdate

							union corresponding
						%END;

				      select distinct a.bene_id, a.indexdate, b.proc_dt, b.proc as cpt, 0 as inpt
				      from cohort as a 
				         inner join der.allcpt as b
				            on a.bene_id = b.bene_id and a.indexdate-&bldays<=b.proc_dt<=a.indexdate
		         ) 
				group by bene_id, indexdate, proc_dt, cpt
			)
			group by bene_id, indexdate, inpt, cpt
			order by bene_id, indexdate, inpt, cpt;
   quit;
%mend;
%allproc();


/* PRESCRIPTIONS */
%macro get_rx(startYr=&startYear, endYr=&endYear, atcGroup=3);
   proc sql;
      create table atc&atcGroup._outpt as
         select distinct a.bene_id, a.indexdate, a.atc&atcGroup as atc&atcGroup._outpt, b.atc_label as label, count(distinct rx_date) as count
         from (%DO yr=&startYr %TO &endYr;
            select distinct a.bene_id, a.indexdate, b.srvc_dt as rx_date format=date9., 
            substr(c.atc,1, (&atcGroup+1)) /*4/6/2024 Tian added +1 to fix the semantics issue*/ as atc&atcGroup
            from cohort as a
               inner join raw.pde_saf_file&yr as b on a.bene_id=b.bene_id and a.indexdate-&bldays<=b.srvc_dt<=a.indexdate
               inner join atc.atc_ndc as c on substr(b.prdsrvid,1,9) = c.ndc9
            %IF &yr<&endYr %THEN union all corresponding; %END;
         ) as a left join atc.atc_ndc as b on a.atc&atcGroup. = b.atc
		group by bene_id, indexdate, atc&atcGroup._outpt
		order by bene_id, indexdate, atc&atcGroup._outpt;
   quit;
%mend;
%get_rx(atcGroup=3);
%get_rx(atcGroup=4);

data all_dx3_inpt(rename=(dx3=dx3_inpt)) all_dx3_outpt(rename=(dx3=dx3_outpt)); 
	set dx3; 
	if inpt=1 then output all_dx3_inpt; else output all_dx3_outpt;
run;

data all_dx4_inpt(rename=(dx4=dx4_inpt)) all_dx4_outpt(rename=(dx4=dx4_outpt)); 
	set dx4; 
	if inpt=1 then output all_dx4_inpt; else output all_dx4_outpt;
run;

data all_cpt5_inpt(rename=(cpt=cpt5_inpt)) all_cpt5_outpt(rename=(cpt=cpt5_outpt)); 
	set cpt;
	if inpt=1 then output all_cpt5_inpt; else output all_cpt5_outpt;
run;

proc datasets lib=work nolist nodetails; 
	delete dx3 dx4 cpt; run;
	change atc3_outpt = all_atc3_outpt atc4_outpt = all_atc4_outpt;
run; quit;


/****************************************************************************************/
/* STEP 2 & 3 (hdPS)  Identify Candidate Empirical Covariates & Assess Recurrence       */    
/****************************************************************************************/
/*%let dimension=dx3_inpt;%let max_num=200; %let min_prev=0.01;
proc print data=&dimension(obs=10);run;
proc print data=&dimension._freq;run;
proc print data=&dimension._freq2;run;
proc print data=&dimension._freq3;run;*/

%macro hdcov(dimension=, max_num=200, min_prev=0.01);

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
		*if count >=100; /*9/28/24, Tian added this according to the hdPS method paper Page 4, Step 2:
		If fewer than 100 patients were identified with a covariate, the covariate was dropped.*/
		pop_percent = count / &N_cohort;
		if pop_percent>=0.5 then do; pop_percent=1-pop_percent; count=&N_cohort-count; end;
		keep &dimension pop_percent count;
	run;

	proc sort data=&dimension._freq2; by descending pop_percent; run;

	data &dimension._freq3;
		set &dimension._freq2;
		prev_order=_N_;
		if _N_<=&max_num /*and pop_percent>=&min_prev*/ then output;
	run;

	*create ordinal variable on cohort dataset;
	*get median, 75th and 90th percentiles -- CURRENTLY THESE VALUES ARE THE MEDIAN/Q3 OF COUNTS AMONG THOSE WITH AT LEAST ONE DX/CPT/RX!;
	proc means data = &dimension         median        q3    p90; class &dimension; var count; 	
		 output out = &dimension._median median=median q3=q3 p90=p90; 
	run;
	
	proc sql;
		create table &dimension._median2 as select a.* from &dimension._median as a 
			  inner join &dimension._freq3   as b on a.&dimension = b.&dimension; 
	quit;

/*--------------6/6/2024 this is the same as the hdiCF manuscript described--------------*/
	proc sql;
		create table &dimension._ordinalVar as select distinct a.*, b.median, b.q3, b.p90,

			case when        a.count = 0 then 0	  
				 when        a.count = 1 then 1
		 		 when    1 < a.count <= b.q3 OR (b.q3=1 and a.count<=b.p90) then 2
				 when b.q3 < a.count then 3
				 else . end as ordinal_val,

			/*create variables for HDPS, 9/18/2024: revised > to >=*/
			case when        a.count >= 1        then 1 else 0 end as cov_once,
			case when        a.count >= b.median then 1 else 0 end as cov_sporadic,
			case when        a.count >= b.q3     then 1 else 0 end as cov_frequent
    
		from &dimension as a inner join &dimension._median2 as b on a.&dimension = b.&dimension;
	quit;
/*---------------------------------------------------------------------------------------*/

	*get label from original file;
	proc sql noprint; 
		select distinct a.&dimension, b.label into :&dimension.1-:&dimension.&max_num, :label1-:label&max_num
			  from &dimension._freq3 as a left join &dimension as b on a.&dimension=b.&dimension; 
			  %LET Ndimension=&SqlObs; 
	quit;

	proc sql;
		create table cohort_&dimension as select distinct a.bene_id, a.indexdate,
			%DO dim=1 %TO &Ndimension; 
				max(case when &dimension = "&&&dimension.&dim" then ordinal_val  else 0 end) as &dimension._&&&dimension.&dim label="&&label&dim, &dimension",			
				max(case when &dimension = "&&&dimension.&dim" then cov_once     else 0 end) as &dimension._once_&&&dimension.&dim,
				max(case when &dimension = "&&&dimension.&dim" then cov_sporadic else 0 end) as &dimension._sporadic_&&&dimension.&dim,
				max(case when &dimension = "&&&dimension.&dim" then cov_frequent else 0 end) as &dimension._frequent_&&&dimension.&dim
				%IF &dim<&Ndimension %THEN ,; %END;
		from cohort as a left join &dimension._ordinalVar as b on a.bene_id=b.bene_id and a.indexdate=b.indexdate
		group by a.bene_id, a.indexdate
		order by bene_id, indexdate;
	quit;

	*Pull reference file for dimensions;
	proc sql;
		create table ref_&dimension  as 
			select distinct "&dimension" as dimension, 
							a.&dimension  as code, 
							a.prev_order, a.pop_percent as prevalence, a.count, b.median, b.q3, b.p90
			from &dimension._freq3 as a left join &dimension._median2 as b on a.&dimension=b.&dimension
			order by prev_order;
	quit;

	*delete intermediate files;
	proc datasets lib=work nolist nodetails; delete &dimension.:; run; quit;
%mend;

/*a procedure for Step 3 Assess Recurrence of hdPS (ie, "If any of the values were equal, the
variable representing the higher cutpoint was dropped." is implemented by the %drop_vars 
in the sister program "select_variables_for_hdPS.sas")*/

*Get N from cohort;

proc sql noprint; select count(*) into :N_cohort from cohort; quit;


*Run macro for each dimension;
%hdcov(dimension=dx3_inpt);
%hdcov(dimension=dx3_outpt);
%hdcov(dimension=dx4_inpt);
%hdcov(dimension=dx4_outpt);
%hdcov(dimension=cpt5_inpt);
%hdcov(dimension=cpt5_outpt);
%hdcov(dimension=atc3_outpt);
/*4/6/2024, added atc5 by Tian because the preivous program doesn't contain "atc5" (actually 4th level) */
%hdcov(dimension=atc4_outpt);

/*proc print data=ana.dx3_inpt_median2; where q3=median;run;
proc print data=dx3_inpt_ordinalvar; where q3=median;run;
proc contents data=dx3_inpt_ordinalvar;run;*/


data out.hdCov_&drug1.v&drug2.;
  merge cohort_dx3_inpt   cohort_dx3_outpt 
		  cohort_dx4_inpt   cohort_dx4_outpt
		  cohort_cpt5_inpt  cohort_cpt5_outpt 
		  cohort_atc3_outpt
		  cohort_atc4_outpt;
	by bene_id indexdate;
run;

data out.ref_hdCov_&drug1.v&drug2;
	length dimension $10 code $7;
	set ref_dx3_inpt ref_dx3_outpt ref_dx4_inpt ref_dx4_outpt 
		ref_cpt5_inpt ref_cpt5_outpt 
	    ref_atc3_outpt ref_atc4_outpt/*4/6/2024 added by Tian*/;
run;


*add descriptive labels to the variables;
proc contents data=out.hdCOV_&drug1.v&drug2 noprint out=labels; run;

data labels2; set labels; length label2 $250;
	label2 = strip(scan(label,2,':')) || ', ' || strip(scan(label,1,'_')) || ' ' || strip(scan(label,2,'_ '));
	keep name label label2;
run;
%macro relabel();
	proc sql noprint; 
		select distinct name, label2 into :var1-:var1000, :label1-:label1000 from labels2;
		%LET N=&SqlObs;
	quit;

	proc datasets lib=out nolist nodetails; modify hdCov_&drug1.v&drug2;
		label %DO i=1 %TO &N; &&var&i = "&&label&i" %END;;
	run;quit;
%mend;
%relabel()


