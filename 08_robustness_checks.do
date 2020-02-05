*** 8. Robustness Checks ***

/*
Here we test whether our Pareto's coefficients are equal between pt and sp and
stable across p95 and p99. Furthermore, we cross-check our Pareto-estimation
with Stata's paretofit().
*/

********************************************************************************
*	OUTLINE
*	8.1: Hausman (suest + mi dataset)
*	8.2: Visual Assessment of Overlapping CIs
*	8.3: Two-sample t-test with unequal variances (Welch's t-test)
*	8.4: Cross-Check with Stata's paretofit()
********************************************************************************



********************************************************************************
*
*	8.1 Hausman-test suest (seemingly unrelated regressions)
*
********************************************************************************
/*
based on: https://stats.idre.ucla.edu/stata/code/comparing-regression-coefficients-across-groups-using-suest/
and https://www.stata.com/support/faqs/statistics/combine-results-with-multiply-imputed-data/#suest

Note:
To cross-check our calculation of CIs in 07_fitting_pareto.do, we again estimate 
CIs with Stata's command 'mi estimate'. This is only for our internal robustness check.
We test H0: alpha_sp - alpha_pt = 0
*/

use "${outpath}soep_pretest_2_MI.dta", clear

* define mysuest program
cap program drop mysuest
program mysuest, eclass properties(mi)
		// USAGE: mysuest "NameOfFirstReg" "FirstReg" "NameOfSecondReg" "SecondReg"

        version 14.2
		args data1 model1 data2 model2

        qui `model1'
        estimates store `data1'
        qui `model2'
        estimates store `data2'
        suest `data1' `data2'
        estimates drop `data1' `data2'
        
        ereturn local title "Seemingly unrelated estimation"
end

*** Test mi estimate
scalar sc_lb95 = (sc_thres_1_p95 + sc_thres_2_p95 + sc_thres_3_p95 + sc_thres_4_p95 + sc_thres_5_p95)/m
scalar sc_lb99 = (sc_thres_1_p99 + sc_thres_2_p99 + sc_thres_3_p99 + sc_thres_4_p99 + sc_thres_5_p99)/m

* "esampvaryok" necessary due to variation in observations across imputations.
mi estimate, esampvaryok: reg lnP_sp ln_nw if(nw >= sc_lb95 & D_pt == 0) [iw=W]
mi estimate, esampvaryok: reg lnP_pt ln_nw if(nw >= sc_lb95 & D_pt == 1) [iw=W]

mi estimate, esampvaryok: reg lnP_sp ln_nw if(nw >= sc_lb99 & D_pt == 0) [iw=W]
mi estimate, esampvaryok: reg lnP_pt ln_nw if(nw >= sc_lb99 & D_pt == 1) [iw=W]

*** Run mi estimate with mysuest
* Lower bound p95
mi estimate, esampvaryok: mysuest "soep" "reg lnP_sp ln_nw if(nw >= sc_lb95 & D_pt == 0) [iw=W]" ///
		"pretest" "reg lnP_pt ln_nw if(nw >= sc_lb95 & D_pt == 1) [iw=W]"
qui matrix mat_tmp =  e(b_mi)
scalar sc_alpha_sp_95_mean = mat_tmp[1,1]
scalar sc_alpha_pt_95_mean = mat_tmp[1,4]

qui matrix mat_tmp = r(table)
scalar sc_alpha_sp_95_ci_low = mat_tmp[5,1]
scalar sc_alpha_sp_95_ci_upp = mat_tmp[6,1]

scalar sc_alpha_pt_95_ci_low = mat_tmp[5,4]
scalar sc_alpha_pt_95_ci_upp = mat_tmp[6,4]

mi estimate, vartable nocitable // Variance following Rubin's rule

* Lower bound p99
mi estimate, esampvaryok: mysuest "soep" "reg lnP_sp ln_nw if(nw >= sc_lb99 & D_pt == 0) [iw=W]" ///
		"pretest" "reg lnP_pt ln_nw if(nw >= sc_lb99 & D_pt == 1) [iw=W]"
		* note: only iweight seems to work in mi estimate + suest, otherwise error. 
qui matrix mat_tmp =  e(b_mi)
scalar sc_alpha_sp_99_mean = mat_tmp[1,1]
scalar sc_alpha_pt_99_mean = mat_tmp[1,4]
		
qui matrix mat_tmp = r(table)
scalar sc_alpha_sp_99_ci_low = mat_tmp[5,1]
scalar sc_alpha_sp_99_ci_upp = mat_tmp[6,1]

scalar sc_alpha_pt_99_ci_low = mat_tmp[5,4]
scalar sc_alpha_pt_99_ci_upp = mat_tmp[6,4]

mi estimate, vartable nocitable // Variance following Rubin's rule

foreach type in pt sp {
	foreach pct in 95 99 {
		di in red "alpha (`type', `pct', mean): `=sc_alpha_`type'_`pct'_mean'"
		foreach ci in low upp {
			di in red "alpha (`type', `pct', `ci'): `=sc_alpha_`type'_`pct'_ci_`ci''"
		}
	}
}

*** Run mi estimate with mysuest and test difference accross coefficients

foreach pct in 95 99 {
	di in red 80 * "*"
	di in red "* Test of coeffecients (`pct'th percentile)"
	mi estimate (diff: [soep_mean]ln_nw - [pretest_mean]ln_nw), esampvaryok nocoef: ///
		mysuest "soep" "reg lnP_sp ln_nw if(nw >= sc_lb`pct' & D_pt == 0) [iw=W]"   ///
				"pretest" "reg lnP_pt ln_nw if(nw >= sc_lb`pct' & D_pt == 1) [iw=W]"
	mi testtransform diff
	scalar sc_pval_`pct' = r(p)
	scalar sc_F_`pct' = r(F)	
}
/*
RESULTS:	
at p95 and p99: pvalue = 0.000. 
-> 	Reject H0: alpha_pt == alpha_sp at 99% confidence level for p95 and p99
*/


*** SOEP vs SOEP and Pretest vs Pretest

* SOEP vs SOEP
mi estimate, esampvaryok: mysuest "soep95" "reg lnP_sp ln_nw if(nw >= sc_lb95 & D_pt == 0) [iw=W]" ///
		"soep99" "reg lnP_sp ln_nw if(nw >= sc_lb99 & D_pt == 0) [iw=W]"

* Pretest vs Pretest 
mi estimate, esampvaryok: mysuest "pretest95" "reg lnP_pt ln_nw if(nw >= sc_lb95 & D_pt == 1) [iw=W]" ///
		"pretest99" "reg lnP_pt ln_nw if(nw >= sc_lb99 & D_pt == 1) [iw=W]"

*** Run mi estimate with mysuest and test difference accross coefficients

foreach data in sp pt {
	if "`data'" == "sp" {
		local D_pt_loc = 0
	}
	if "`data'" == "pt" {
		local D_pt_loc = 1
	}
	di in red 80 * "*"
	di in red "* Test of coeffecients (`pct'th percentile)"
	mi estimate (diff: [reg_`data'_95_mean]ln_nw - [reg`data'_99_mean]ln_nw), esampvaryok nocoef: 	///
		mysuest "reg_`data'_95" "reg lnP_`data' ln_nw if(nw >= sc_lb95 & D_pt == `D_pt_loc') [iw=W]"      	///
				"reg`data'_99" "reg lnP_`data' ln_nw if(nw >= sc_lb99 & D_pt == `D_pt_loc') [iw=W]"
	mi testtransform diff
	scalar sc_pval_`pct' = r(p)
	scalar sc_F_`pct' = r(F)	
}

/*
RESULTS:	SOEP vs SOEP:
			->	diff: .1155249 
			->	pvalue = 0.0448: Evidence to reject H_0 that
				alpha_sp95 == alpha_sp99 at 95% confidence level
			
			Pretest vs Pretest:
			->	diff:  .1244399
			->	pvalue = 0.000: Strong evidence to reject H_0 that 
				alpha_pt95 == alpha_pt99 at 99% confidence level
*/

********************************************************************************
*
*	8.2. Visual Assessment of Overlapping CIs 
*		 (in paper: Figure 4.1)
*
********************************************************************************

use "${outpath}soep_pretest_2.dta", clear

set graph off

local pctile p95 p99

/*
template of dataset used to plot CIs:
	/ sc_cilo_p95_sp / sc_cihi_p95_sp / sc_avgalpha_p95_sp /
	/ sc_cilo_p99_sp / sc_cihi_p99_sp / sc_avgalpha_p99_sp /
	/ sc_cilo_p95_pt / sc_cihi_p95_pt / sc_avgalpha_p95_pt /
	/ sc_cilo_p99_pt / sc_cihi_p99_pt / sc_avgalpha_p99_pt /
*/

* save CIs to matrix A
matrix A = J(4, 3 ,.)
matrix A[1,1] = `=sc_cilo_p95_sp'
matrix A[2,1] = `=sc_cilo_p99_sp'
matrix A[3,1] = `=sc_cilo_p95_pt'
matrix A[4,1] = `=sc_cilo_p99_pt'
matrix A[1,2] = `=sc_cihi_p95_sp'
matrix A[2,2] = `=sc_cihi_p99_sp'
matrix A[3,2] = `=sc_cihi_p95_pt'
matrix A[4,2] = `=sc_cihi_p99_pt'
matrix A[1,3] = `=sc_avgalpha_p95_sp'
matrix A[2,3] = `=sc_avgalpha_p99_sp'
matrix A[3,3] = `=sc_avgalpha_p95_pt'
matrix A[4,3] = `=sc_avgalpha_p99_pt'
matrix list A

* transform matrix to Stata-file and process
drop _all
svmat double A
matrix drop A

ren (A1 A2 A3) (cilo cihi avgalpha)
gen X=_n
gen pctile = "p95"
gen data = "SOEP"
replace data = "Pretest" if inlist(X,3,4)
replace pctile = "p99" if inlist(X,2,4)
replace avgalpha = round(avgalpha,0.01)

* plot CIs
foreach pct of local pctile {

	if "`pct'" == "p95" {
	local yscale "yscale(range(0 4))"
	}
	if "`pct'" == "p99" {
	local yscale "yscale(range(1 5))"
	}

	twoway rcap cihi cilo X if pctile=="`pct'", horizontal lstyle(ci) ||			///
		scatter X avgalpha 	if pctile=="`pct'", mstyle(p1) mlabsize(medium) mlabel(avgalpha) mlabposition(6) || ///
		scatter X avgalpha 	if pctile=="`pct'", mstyle(p1) mlabsize(medium) mlabel(data) mlabposition(12)		///
		legend(order(1 "CI (95%)" 2 "average of Pareto's Alpha")) ylabel(, nogrid)	///
		xtitle("alpha") ytitle("") `yscale' xscale(range(0.75 2.25)) 				///
		xlabel(, grid) yscale(off) scheme(s2mono) graphregion(color(white))

	graph export "${graphs}08_overlapping_CIs_`pct'.pdf", replace
}

set graph on


********************************************************************************
*		
*	8.3: Two-sample t-test with unequal variances (Welch's t-test)
*		
********************************************************************************

/*
NOTE:
Here we test first by-hand and then using Statas built-in ttesti if the Pareto's
coefficients are stable:
alpha(sp, p95) == alpha(pt, p95)	and		alpha(sp, p95) == alpha(sp, p99)
alpha(sp, p99) == alpha(pt, p99)	and		alpha(pt, p95) == alpha(pt, p99)

RESULTS:
All tests deliver p-values of 0.0000. The results achieved "by-hand" match those 
by Stata's ttesti function.
*/

* Preparation: Calculate average number of observations for degrees of freedom		
qui scalar sc_n_p95_sp = round(1/5 * (sc_n_1_p95_sp + sc_n_2_p95_sp + sc_n_3_p95_sp + sc_n_4_p95_sp + sc_n_5_p95_sp))
qui scalar sc_n_p95_pt = round(1/5 * (sc_n_1_p95_pt + sc_n_2_p95_pt + sc_n_3_p95_pt + sc_n_4_p95_pt + sc_n_5_p95_pt))
qui scalar sc_n_p99_sp = round(1/5 * (sc_n_1_p99_sp + sc_n_2_p99_sp + sc_n_3_p99_sp + sc_n_4_p99_sp + sc_n_5_p99_sp))
qui scalar sc_n_p99_pt = round(1/5 * (sc_n_1_p99_pt + sc_n_2_p99_pt + sc_n_3_p99_pt + sc_n_4_p99_pt + sc_n_5_p99_pt))


* "By hand* calculations of Two Sample test statistic and degrees of freedom (within dataset / accross thresholds)
foreach data in sp pt {
	scalar sc_2Sttest_`data' = (sc_avgalpha_p95_`data' - sc_avgalpha_p99_`data') /  sqrt( (sc_tv_p95_`data')/(sc_n_p95_`data') + (sc_tv_p99_`data')/(sc_n_p99_`data') )
	scalar sc_2Sttest_df_`data' = (((sc_tv_p95_`data' / sc_n_p95_`data') + (sc_tv_p99_`data' / sc_n_p99_`data'))^2) / ///
				((sc_tv_p95_`data'^2 / (sc_n_p95_`data'^2 * (sc_n_p95_`data' - 1))  ) + (sc_tv_p99_`data'^2 / (sc_n_p99_`data'^2 * (sc_n_p99_`data' - 1))))
	di "sc_2Sttest_ (`data'):"sc_2Sttest_`data'
	di "sc_2Sttest_df_(`data'): "  sc_2Sttest_df_`data'
}

* "By hand* calculations of Two Sample test statistic and degrees of freedom (within thresholds / accross datasets)
foreach pct in p95 p99 {
	scalar sc_2Sttest_`pct' = (sc_avgalpha_`pct'_sp - sc_avgalpha_`pct'_pt) /  sqrt( (sc_tv_`pct'_sp)/(sc_n_`pct'_sp) + (sc_tv_`pct'_pt)/(sc_n_`pct'_pt) )
	scalar sc_2Sttest_df_`pct' = (((sc_tv_`pct'_sp / sc_n_`pct'_sp) + (sc_tv_`pct'_pt / sc_n_`pct'_pt))^2) / ///
				((sc_tv_`pct'_sp^2 / (sc_n_`pct'_sp^2 * (sc_n_`pct'_sp - 1))  ) + (sc_tv_`pct'_pt^2 / (sc_n_`pct'_pt^2 * (sc_n_`pct'_pt - 1))))
	di "sc_2Sttest_ (`pct')): "sc_2Sttest_`pct'
	di "sc_2Sttest_df_(`pct'): " sc_2Sttest_df_`pct'
}

* Same calculations as above using stata buitin "ttesti, ..., unequal" (within dataset / accross thresholds)
foreach data in sp pt {
	foreach pct in p95 p99 {
		* generate std deviation from total variance for each pct and data
		scalar sc_sd_`pct'_`data' = sqrt(`=sc_tv_`pct'_`data'')
	}
	
	* two sample ttest with unequal variance
	ttesti `=sc_n_p95_`data'' `=sc_avgalpha_p95_`data'' `=sc_sd_p95_`data'' `=sc_n_p99_`data'' `=sc_avgalpha_p99_`data'' `=sc_sd_p99_`data'',  unequal
}

* Same calculations as above using stata buitin "ttesti, ..., unequal" (within thresholds / accross datasets)
foreach pct in p95 p99 {
	ttesti `=sc_n_`pct'_sp' `=sc_avgalpha_`pct'_sp' `=sc_sd_`pct'_sp' `=sc_n_`pct'_pt' `=sc_avgalpha_`pct'_pt' `=sc_sd_`pct'_pt',  unequal
}


********************************************************************************
*		
*	8.4: Cross-Check with Stata's paretofit() 
*		
********************************************************************************

/*
NOTE:
Here we run the Stata's command paretofit() over p95 and p99 to cross-check our
estimations in 07_fitting_pareto.do. 

RESULTS:
The results are similar to ours.
*/

use "${outpath}soep_pretest_2.dta", clear

* define threshold levels as percentiles of SOEP
local pctile p95 p99

* define datasets: SOEP (sp), Pretest (pt)
local data pt sp

foreach dat of local data {
	forval imp=1(1)`=m' {
		foreach pc of local pctile {
			* thresholds from SOEP
			qui sum _`imp'_nw if D_pt==0, d
			scalar sc_thres_`imp'_`pc'_sp = r(`pc')
			
			di in red 80 * "="
			di in red "Data: `dat'. Imputation: `imp'. Percentile: `pc'. Lower bound: `=sc_thres_`imp'_`pc''"
			di in red 80 * "="
			
			if "`dat'" == "sp" {
				paretofit _`imp'_nw if D_pt == 0 [fweight=round(W_`dat')], x0(`=sc_thres_`imp'_`pc'_sp') robust
				scalar sc_alpha_pfit_`imp'_`pc'_`dat' = e(ba)
			}
			if "`dat'" == "pt" {
				qui paretofit _`imp'_nw if D_pt == 1 [fweight=round(W_`dat')], x0(`=sc_thres_`imp'_`pc'_sp') robust
				scalar sc_alpha_pfit_`imp'_`pc'_`dat' = e(ba)
			}
		}
	}
}

foreach dat in pt sp {
	foreach pc of local pctile {
		di in red  "alpha (pareto fit, `pc', `dat', mean): " (sc_alpha_pfit_1_`pc'_`dat' + sc_alpha_pfit_2_`pc'_`dat' + /// 
		sc_alpha_pfit_3_`pc'_`dat' + sc_alpha_pfit_4_`pc'_`dat' +  sc_alpha_pfit_5_`pc'_`dat') / m
	}
}

***
