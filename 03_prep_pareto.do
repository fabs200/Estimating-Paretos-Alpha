***	3. Prepare Pareto Distribution ***

/*
Here we generate the CDF (F), the CCDF (P=1-F) and take logs to generate ln(P) 
for each implicate of both datasets, Pretest and SOEP.
*/


use "${outpath}soep_pretest_1.dta", clear

* loop over all imputed net wealth variables

forval imp=1(1)5 {	

	* generate log net wealth
	gen ln_nw_`imp'_ = ln(_`imp'_nw)
	
	****************************************************************************	
	*		
	*	3.1 Pretest (_pt)
	*		
	****************************************************************************
	
	sort D_pt _`imp'_nw
	qui sum W_pt if D_pt==1
	scalar sc_N_`imp'_pt_ = r(sum)
	gen F_pt_`imp'_ = sum(W_pt)/sc_N_`imp'_pt_ if D_pt==1
	replace F_pt_`imp'_ = round(F_pt_`imp'_, 1e-12)
	gen P_pt_`imp'_ = 1 - F_pt_`imp'_ if D_pt==1
	gen lnP_pt_`imp'_ = ln(P_pt_`imp'_) if D_pt==1
	
	/*
	NOTE: round()
	Without rounding, the highest value of F is displayed in Stata as: 1, but 
	internally, it is saved as a value 0.99999999... If we then calculate ln(1-F), 
	we get an extreme value of around -30. All other values of ln(1-F) range from 
	0.0 to around -12.0. This would result in a very distorted estimation.
	To refrain from data editing (manually assigning a certain value of F 
	(e.g. 0.99999999.... -> 0.9999) to the richest obs.), we decided to round 
	F to 12th decimal point.
	*/

	****************************************************************************
	*		
	*	3.2 SOEP (_sp)
	*		
	****************************************************************************

	sort D_pt _`imp'_nw
	qui sum W_sp if D_pt==0
	scalar sc_N_`imp'_sp_ = r(sum)
	gen F_sp_`imp'_ = sum(W_sp)/sc_N_`imp'_sp_ if D_pt==0
	replace F_sp_`imp'_ = round(F_sp_`imp'_, 1e-12)
	gen P_sp_`imp'_ = 1 - F_sp_`imp'_ if D_pt==0
	gen lnP_sp_`imp'_ = ln(P_sp_`imp'_) if D_pt==0
}

* check all generated variables
foreach data in sp pt {
	di in red _newline "+++++ `data' +++++"
	sum ln_nw_*_
	sum F_`data'_*_
	sum lnP_`data'_*_
	sum P_`data'_*_
}

* save dataset
save "${outpath}soep_pretest_2.dta", replace

***
