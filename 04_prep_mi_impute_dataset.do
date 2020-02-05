***	4. Generate MI dataset ***

/*
In this file another version of the dataset is generated in order to use Stata's 
internal capabilities of dealing with imputed datasets. We only use this for our 
internal robustness check and do not refer to it in our paper.
*/

use "${outpath}soep_pretest_2.dta", clear

********************************************************************************	
*
* 4.1 Preparation for MI dataset
*
********************************************************************************

foreach type in pt sp {
	forval imp=1(1)5 {

		/* 
		NOTE: the variables with _#_ in the name (at end or start) seems to be 
		protected by Stata when converting into mi dataset. Therefore, the last 
		underscore "_" is removed on the generated variables, so that it can be 
		later converted into a mi passive variable.
		*/
		
		rename F_`type'_`imp'_ F_`type'_`imp'
		drop P_`type'_`imp'_ lnP_`type'_`imp'_

	}
}

***

* ln_nw_#_ also must be dropped (or renamed) in order to convert into mi dataset
drop ln_nw_*_

* _mi_miss seems to be necessary to unset MI. It indicates if obs is missing
* if imputed variable are all equal --> not missing
gen _mi_miss = 1
replace _mi_miss = 0 if _1_nw == _2_nw & _1_nw == _3_nw & _1_nw == _4_nw & _1_nw == _5_nw 

* generate "raw" nw, that is if not missing
gen nw = .
replace nw = _1_nw if _mi_miss == 0

mi unset 		// dataset was previously defined as mi dataset
drop mi_miss
rename (_1_nw _2_nw _3_nw _4_nw _5_nw) (nw1 nw2 nw3 nw4 nw5)	
* renaming seems to be necesssary because "mi import" tries to redefine _1_nw, 
* _2_nw... variables

* recreate dataset indicating which var are imputed
mi import wide, imputed(nw=nw1 nw2 nw3 nw4 nw5) drop clear

********************************************************************************	
*
*	4.2 Generate mi impute variables
*
********************************************************************************

* with "mi passive" imputed variables are automatically generated
mi passive: gen ln_nw = ln(nw) 

********************************************************************************	
* 4.2 a) Pretest (_pt)
********************************************************************************

gen F_pt = .
replace F_pt = F_pt_1 if _mi_miss == 0
mi register passive F_pt

forval imp=1(1)5 {
	replace _`imp'_F_pt = F_pt_`imp'
}
mi passive: gen P_pt = 1 - F_pt
mi passive: gen lnP_pt = ln(P_pt)


********************************************************************************	
*
* 4.2 b) SOEP (_sp)
*
********************************************************************************

gen F_sp = .
replace F_sp = F_sp_1 if _mi_miss == 0
mi register passive F_sp

forval imp=1(1)5 {
	replace _`imp'_F_sp = F_sp_`imp'
}
mi passive: gen P_sp = 1 - F_sp
mi passive: gen lnP_sp = ln(P_sp)


local varlist W lnP
foreach var in `varlist' {
	* generate unique Weight and lnP variable for both sp and pt. 
	cap gen `var' = `var'_sp
	replace `var' = `var'_pt if missing(`var')
}

*** save dataset
save "${outpath}soep_pretest_2_MI.dta", replace
* this dataset we only use in 08_robustness_checks.do

***
