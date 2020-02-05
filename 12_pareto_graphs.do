
*** 12. Pareto Graphs ***

/*
Here we produce: Figure 3.1, Figure 3.2, Figure 4.2, Figure 4.3, Figure A.1-6
*/

********************************************************************************
*	OUTLINE
*	12.1 Empirical Pareto-Distribution
*	12.2 Empirical Pareto-Distribution with Fitted Line
*	12.3 Theoretical Pareto-Distribution
********************************************************************************

set graph off
use "${outpath}soep_pretest_3.dta", clear

*** Preparation for Plots ***

* define threshold levels as percentiles of SOEP
local pctile p95 p99

* define datasets: SOEP (_sp) and Pretest (_pt)
local data pt sp

* m: number of imputations
scalar m = 5

* create png folder
capture mkdir ${graphs}png

*** prepare scalars, vars in logs, mean of cilo, cihi
forval imp=1(1)`=m' {
	foreach pct of local pctile {
		
		** threshold scalar in logs
		scalar sc_lnthres_`imp'_`pct' = ln(sc_thres_`imp'_`pct')
		if "${show_all}" == "TRUE" {
			di in red "threshold (`pct', imp`imp'): 	`=sc_thres_`imp'_`pct''"
			di in red "ln(threshold) (`pct', imp`imp'): 	`=sc_lnthres_`imp'_`pct''"
		}
		
		** nw_hat in logs
		* with Pretest Alpha
		gen lnnw_hat_`imp'_`pct'_cilo = ln(nw_hat_`imp'_`pct'_cilo)
		gen lnnw_hat_`imp'_`pct'_cihi = ln(nw_hat_`imp'_`pct'_cihi)
		* with SOEP Alpha
		gen lnnw_hat2_`imp'_`pct'_cilo = ln(nw_hat2_`imp'_`pct'_cilo)
		gen lnnw_hat2_`imp'_`pct'_cihi = ln(nw_hat2_`imp'_`pct'_cihi)

		** mean of cilo, ciho
		* with Pretest Alpha
		egen lnnw_hat_`imp'_`pct'_mean = rowmean(lnnw_hat_`imp'_`pct'_cilo lnnw_hat_`imp'_`pct'_cihi)
		* with SOEP Alpha
		egen lnnw_hat2_`imp'_`pct'_mean = rowmean(lnnw_hat2_`imp'_`pct'_cilo lnnw_hat2_`imp'_`pct'_cihi)

		** regression line (reglin) based on estimated alphas (y=-a*x+c)
		*  (cilopt, cihipt: alpha from Pretest, cilosp, cihisp: alpha from SOEP)
		gen reglin_`imp'_`pct'_cilopt = -sc_cilo_`pct'_pt*ln_nw_`imp'_ + sc_cons_`imp'_`pct'_sp if ln_nw_`imp'_>= sc_lnthres_`imp'_`pct' & D_pt==1
		gen reglin_`imp'_`pct'_cihipt = -sc_cihi_`pct'_pt*ln_nw_`imp'_ + sc_cons_`imp'_`pct'_sp if ln_nw_`imp'_>= sc_lnthres_`imp'_`pct' & D_pt==1
		* mean of cilo, cihi (alpha pretest)
		egen reglin_`imp'_`pct'_meanpt = rowmean(reglin_`imp'_`pct'_cilopt reglin_`imp'_`pct'_cihipt)
		gen reglin_`imp'_`pct'_cilosp = -sc_cilo_`pct'_sp*ln_nw_`imp'_ + sc_cons_`imp'_`pct'_sp if ln_nw_`imp'_>= sc_lnthres_`imp'_`pct' & D_pt==0
		gen reglin_`imp'_`pct'_cihisp = -sc_cihi_`pct'_sp*ln_nw_`imp'_ + sc_cons_`imp'_`pct'_sp if ln_nw_`imp'_>= sc_lnthres_`imp'_`pct' & D_pt==0
		* mean of cilo, cihi (alpha soep)
		egen reglin_`imp'_`pct'_meansp = rowmean(reglin_`imp'_`pct'_meanpt reglin_`imp'_`pct'_meanpt)

	}
}

*** define dataset-specific locals with graph cosmetics
local symbol_opt "jitter(1.25) msize(small) msymbol(oh) mlwidth(thin)"
local line_opt " lpattern(_--) lwidth(thin)"

* loop over three datasets
foreach dat of local data {
	
	if "`dat'" == "sp" {
		local g_title "SOEP"
		local g_nw_opt "D_pt==0"
		local g_src "SOEP 2012 (v33.1)"
		local g_note " with applied household weights."
		local g_note2 "         "
		local g_col "blue%80"
	}
	if "`dat'" == "pt" {
		local g_title "Pretest"
		local g_nw_opt "D_pt==1"
		local g_src "Pretest 2017"
		local g_note " with own calculated weights."
		local g_note2 "         "
		local g_col "cranberry%80"

	}
	****************************************************************************
	*
	* 12.1 	Empirical Pareto-Distribution
	*		(in paper: Figure 3.2, Figure A.1-2)
	*
	****************************************************************************

	/* 
	In total 20x graphs: 
		- sp, pt
		- threshold p95, p99
		- 5 imputations + 5 nw_hat
	*/

	****************************************************************************
	*		12.2 a) Pareto Distribution: Scatterplot
	****************************************************************************
	
	forval imp=1(1)`=m' {
		* Plot only 1 net wealth implicat
		twoway ///
		(scatter lnP_`dat'_`imp'_ ln_nw_`imp'_ if `g_nw_opt', sort(ln_nw_`imp'_) `symbol_opt' mlcolor(gs5%50) ///
		ylab(, grid) ytitle("ln(P)") xlab(, grid) xtitle("ln(nw)") scheme(s2mono) graphregion(color(white)) ///
		/* title("Pareto Distribution of the `g_title'") */ 	///
		legend(row(1) order(1 "ln(nw{sub:`imp'})")) 	///
		note("Source: `g_src'." "Note: implicate `imp' of the net wealth variable (nw) displayed,`g_note'"))
		
		*graph export "${graphs}12_1_pareto_distrib_scatter_`dat'_imp`imp'.pdf", replace
		graph export "${graphs}png/12_1_pareto_distrib_scatter_`dat'_imp`imp'.png", replace
	}

	****************************************************************************
	*		
	*	12.2 	Empirical Pareto-Distribution with Fitted Line, 
	*			Actual vs. Predicted, SOEP with Pretest Alpha (12.2a) and with 
	*			SOEP Alpha (12.2b)
	*			(in paper: Figure 4.2, Figure 4.3, Figure A.3-6)
	*		
	****************************************************************************

	/* 
	In total 20x graphs: 
		- sp, pt
		- threshold p95, p99
		- 5 imputations + 5 nw_hat
	*/
	
	if "`dat'" == "sp" {
		
		forval imp=1(1)`=m' {
			foreach pct of local pctile {
				
				* prepare abbrev./rounded values for Note of graphs
				local thres_pct = substr("`pct'",2,3)
				local alpha_cilo_pt = round(`=sc_cilo_`pct'_pt', .01)
				local alpha_cihi_pt = round(`=sc_cihi_`pct'_pt', .11)
				local alpha_mean_pt = round(`=sc_avgalpha_`pct'_pt', .01)
				local alpha_cilo_sp = round(`=sc_cilo_`pct'_sp', .01)
				local alpha_cihi_sp = round(`=sc_cihi_`pct'_sp', .01)				
				local alpha_mean_sp = round(`=sc_avgalpha_`pct'_sp', .01)
				local sc_avgrrmse_`imp'_`pct' = (sc_cilorrmse_`imp'_`pct' + sc_cihirrmse_`imp'_`pct')/2
				local sc_avgrrmse2_`imp'_`pct' = (sc_cilorrmse2_`imp'_`pct' + sc_cihirrmse2_`imp'_`pct')/2
				
				di in red "+++++ actual vs. predicted"
				di in red "+++++ dataset: `dat', imp`imp', threshold: `pct', threshold value: `=sc_thres_`imp'_`pct''"
				di in red "+++++ Pretest alpha (CIlo, CIhi): `alpha_cilo_pt' to `alpha_cihi_pt'"
				di in red "+++++ SOEP alpha (CIlo, CIhi): `alpha_cilo_sp' to `alpha_cihi_sp'"

				*** 12.2a Pareto Distribution with Pretest Alpha
				twoway ///
				(scatter lnP_sp_`imp'_ ln_nw_`imp'_ if `g_nw_opt', sort(ln_nw_`imp'_) `symbol_opt' mcolor(gs4%50) msize(vsmall)) ///
				(line lnP_sp_`imp'_ lnnw_hat_`imp'_`pct'_cilo if ln_nw_`imp'_ >= `=sc_lnthres_`imp'_`pct'' & `g_nw_opt', lcolor(gold%20)) 	///
				(scatter lnP_sp_`imp'_ lnnw_hat_`imp'_`pct'_mean if ln_nw_`imp'_ >= `=sc_lnthres_`imp'_`pct'' & `g_nw_opt', msymbol(o) mcolor(red%20) msize(vsmall))	///
				(line lnP_sp_`imp'_ lnnw_hat_`imp'_`pct'_cihi if ln_nw_`imp'_ >= `=sc_lnthres_`imp'_`pct'' & `g_nw_opt', lcolor(gold%20)	///
				ytitle("ln(P)") xlab(, grid) xtitle("ln(nw)") scheme(s2mono) graphregion(color(white)) ///
				/* title("Pareto Distribution of `g_title'") subtitle("actual vs. predicted data") */ ///
				legend(row(1) order(1 "actual ln({it:nw}{sub:{it:`imp'}})" 3 "predicted ln({it:nw}{sub:{it:`imp'}})" 2 "CI (95%)")) ///
				note("Source: `g_src'." "Note: actual data: net wealth (nw) implicate `imp' (grey),`g_note'" "Threshold at `thres_pct'% percentile of the SOEP (`=sc_thres_`imp'_`pct'' Euros), predicted data: nw above threshold," "predicted with {it:alpha}{sub:{it:Pretest}} between `alpha_cilo_pt' to `alpha_cihi_pt' (yellow), mean of alpha = `alpha_mean_pt' (red), RRMSE = `sc_avgrrmse_`imp'_`pct''%."))

				*graph export "${graphs}12_2a_actuals_vs_predicted_`imp'_`pct'_withpretestalpha.pdf", replace
				graph export "${graphs}png/12_2a_actuals_vs_predicted_`imp'_`pct'_withpretestalpha.png", replace

				*** 12.2b Pareto Distribution with SOEP Alpha
				twoway ///
				(scatter lnP_sp_`imp'_ ln_nw_`imp'_ if `g_nw_opt', sort(ln_nw_`imp'_) `symbol_opt' mcolor(gs4%50) msize(vsmall)) ///
				(line lnP_sp_`imp'_ lnnw_hat2_`imp'_`pct'_cilo if ln_nw_`imp'_ >= `=sc_lnthres_`imp'_`pct'' & `g_nw_opt', lcolor(gold%20)) ///
				(scatter lnP_sp_`imp'_ lnnw_hat2_`imp'_`pct'_mean if ln_nw_`imp'_ >= `=sc_lnthres_`imp'_`pct'' & `g_nw_opt', msymbol(o) mcolor(red%20) msize(vsmall))  ///
				(line lnP_sp_`imp'_ lnnw_hat2_`imp'_`pct'_cihi if ln_nw_`imp'_ >= `=sc_lnthres_`imp'_`pct'' & `g_nw_opt', lcolor(gold%20)  ///
				ytitle("ln(P)") xlab(, grid) xtitle("ln(nw)") scheme(s2mono) graphregion(color(white)) ///
				/* title("Pareto Distribution of `g_title'") subtitle("actual vs. predicted data") */ ///
				legend(row(1) order(1 "actual ln({it:nw}{sub:{it:`imp'}})" 3 "predicted ln({it:nw}{sub:{it:`imp'}})" 2 "CI (95%)")) ///
				note("Source: `g_src'." "Note: actual data: net wealth (nw) implicate `imp' (grey),`g_note'" "Threshold at `thres_pct'% percentile of the SOEP (`=sc_thres_`imp'_`pct'' Euros), predicted data: nw above threshold," "predicted with {it:alpha}{sub:{it:SOEP}} between `alpha_cilo_sp' to `alpha_cihi_sp' (yellow), mean of alpha = `alpha_mean_sp' (red), RRMSE = `sc_avgrrmse2_`imp'_`pct''%."))

				*graph export "${graphs}12_2b_actuals_vs_predicted_`imp'_`pct'_withsoepalpha.pdf", replace
				graph export "${graphs}png/12_2b_actuals_vs_predicted_`imp'_`pct'_withsoepalpha.png", replace

			}
		}
	}
}

********************************************************************************
*
* 12.3 	Theoretical Pareto-Distribution
*		(in paper: Figure 3.1)
*
********************************************************************************

local a1 = 2
local a2 = 3
local a3 = 5
local rl = 1000
local rh = 3000

local ylow  = 1000

local yhig1 = 1000000
local yhig2 = 100000
local yhig3 = 10000 * 1.6
local plotLowBound = 0.01

*** generate obs
drop _all
set obs 2000
egen x = seq(), from(0)
gen x1 = x + `ylow'

gen y1 = (`a1' * `ylow'^`a1') / (x1^(`a1' + 1))
gen y2 = (`a2' * `ylow'^`a2') / (x1^(`a2' + 1))
gen y3 = (`a3' * `ylow'^`a3') / (x1^(`a3' + 1))

gen y4 = (1 - (x1 / `ylow')^(-`a1'))
gen y5 = (1 - (x1 / `ylow')^(-`a2'))
gen y6 = (1 - (x1 / `ylow')^(-`a3'))

*** plot 

#delimit ;
graph twoway (line y1 x1 if x < `rh', lw(medthick))
             (line y2 x1 if x < `rh', lw(medthick))
             (line y3 x1 if x < `rh', lw(medthick)),
  xscale(range(`ylow3'))
  xlabel(,angle(90))

  bgcolor(white) graphregion(color(white)) scheme(s2mono)  
  xline(`ylow', lpattern(dash) lcolor(grey))

  legend(order( 1 "{&alpha} = `a1'" 2 "{&alpha} = `a2'" 3 "{&alpha} = `a3'") rows(1))
  //title("Probability Distribution Function")
  xtitle("y")
  ytitle("f(y|{&alpha}, threshold=`ylow')")
  
  saving(paretoCDF, replace);
#delimit cr 
graph export ${graphs}/th_paretodiag01a.pdf, replace

#delimit ;
graph twoway (line y4 x1, lw(medthick))
             (line y5 x1, lw(medthick))
             (line y6 x1, lw(medthick)), 
  xscale(range(`ylow3'))
  xlabel(,angle(90))
  
  xline(`ylow', lpattern(dash) lcolor(grey))
  bgcolor(white) graphregion(color(white)) scheme(s2mono)                                           
  legend(order( 1 "{&alpha} = `a1'" 2 "{&alpha} = `a2'" 3 "{&alpha} = `a3'") rows(1))
  //title("Cumulative Distribution Function")
  xtitle("y")
  ytitle("F(y|{&alpha}, threshold=`ylow')")
  
  saving(paretoPDF, replace);
#delimit cr
graph export ${graphs}/th_paretodiag01b.pdf, replace


gr combine paretoCDF.gph paretoPDF.gph,  xsize(7) ///
  graphregion(color(white))

graph export ${graphs}/04_paretoDistGraphs.pdf, replace


/* Plot Pareto distributions */

local a1 = 2
local a2 = 3
local a3 = 5
local rl = 1000
local rh = 2000
local ylow  = 1000

local yhig1 = 1000000
local yhig2 = 100000
local yhig3 = 10000 * 1.6
local plotLowBound = 0.01

drop _all
set obs 10000
egen x = seq(), from(0)

gen x1 = `ylow' + (x * 100)
gen x2 = `ylow' + (x * 10)
gen x3 = x + `ylow'

gen y1 = (1 - (1 - (`ylow' / x1)^(`a1')))
gen y2 = (1 - (1 - (`ylow' / x2)^(`a2')))
gen y3 = (1 - (1 - (`ylow' / x3)^(`a3')))


*** Plot CCDF linear and CCDF log-log

* Pareto Diagram (log-log)

#delimit ;
graph twoway (line y1 x1 if y1 > `plotLowBound', lw(medthick)) 
             (line y2 x2 if y2 > `plotLowBound', lw(medthick))  
	     (line y3 x3 if y3 > `plotLowBound', lw(medthick)),  
  xscale(range(`ylow1'))
  ylabel(0.05 0.2(0.2)1 1)
  xlabel(,angle(90))
  xline(`ylow', lpattern(dash) lcolor(grey))

  bgcolor(white) graphregion(color(white))  scheme(s2mono)                                                                  
  legend(order( 1 "{&alpha} = `a1'" 2 "{&alpha} = `a2'" 3 "{&alpha} = `a3'") rows(1))
  //title("Pareto Diagram (log-log)")
  xtitle("y")
  ytitle("1 - F(y|{&alpha}, threshold=`ylow')")

  saving(paretoDiagram01, replace);
#delimit cr
graph export ${graphs}/th_paretodiag02a.pdf, replace


* Pareto Diagram (log-log)

#delimit ;
graph twoway (line y1 x1 if y1 > `plotLowBound', lw(medthick)) 
             (line y2 x2 if y2 > `plotLowBound', lw(medthick))  
	     (line y3 x3 if y3 > `plotLowBound', lw(medthick)),  
  xscale(log r(`ylow2')) yscale(log)
  ylabel(0.05 0.2(0.2)1 1)
  xlabel(,angle(90))
  
  xline(`ylow', lpattern(dash) lcolor(grey))

  bgcolor(white) graphregion(color(white)) scheme(s2mono)                                                                    
  legend(order( 1 "{&alpha} = `a1'" 2 "{&alpha} = `a2'" 3 "{&alpha} = `a3'") rows(1))
  //title("Pareto Diagram (log-log)")
  xtitle("y")
  ytitle("1 - F(y|{&alpha}, threshold=`ylow')")

  saving(paretoDiagram02, replace);
#delimit cr
graph export ${graphs}/th_paretodiag02b.pdf, replace


gr combine paretoDiagram01.gph paretoDiagram02.gph,   xsize(7) ///
  graphregion(color(white))

graph export ${graphs}/04_paretoDiagram.pdf, replace

set graph on


di in red "+++++ END +++++"

***
