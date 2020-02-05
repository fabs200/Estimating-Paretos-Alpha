*** 10. threshold alpha relation ***

/* 
To check the stability of alpha, we estimate each 10.000 euros the corresponding
alpha and plot the resulting alpha-threshold-combination.

Here we produce: Figure 4.4
*/

set type double, permanently
set matsize 11000
set graph on

* m: number of imputations
scalar m = 5

* define datasets: SOEP (sp), Pretest (pt)
local data sp pt

foreach dat of local data {

	forval imp=1(1)5 {

		set graph off

		use "${outpath}soep_pretest_2.dta", clear
	
		* sort
		sort _`imp'_nw

		* calculate the percentiles
		_pctile _`imp'_nw if  _`imp'_nw>0, percentiles(1 99 99.5 99.9)
		local p1	= r(r1)
		local p99	= r(r2)
		local p99_5	= r(r3)
		local p99_9	= r(r4)
		
		scalar nw_min = `p1'
		scalar nw_max = `p99' /*`r(max)': not possible to regress very few obs */

		scalar steps = 10000 		/* <--- define steps, 10.000 euros */

		* define matrix size
		scalar mat_size = ((nw_max - nw_min) / steps) +1
		* define matrix
		matrix P = J(mat_size,5,.)
		matrix list P

		di in red "++++++++++ matsize: `=mat_size', steps: `=steps'" _newline ///
			"           nw_min: `=nw_min', nw_max: `=nw_max'"

		************************************************************************
		* 10.1 Calculate Threshold-Alpha-Combination
		************************************************************************
		di in red "++++ `dat'"

		if "`dat'" == "sp" {
			local sp_opt "& D_pt==0"
			local sp_src "SOEP 2012 (v33.1)"
			local sp_w_info "with applied household weights."
		}
		if "`dat'" == "pt" {
			local pt_opt "& D_pt==1"
			local pt_src "Pretest 2017"
			local pt_w_info "own weights applied."
		}
					
		* loop over multiple thresholds
		scalar a = 1
		forval i = `=nw_min' (`=steps') `=nw_max' {

			* threshold value
			scalar sc_thres = `i'

			qui reg lnP_`dat'_`imp'_ ln_nw_`imp'_ if (_`imp'_nw >= sc_thres ``dat'_opt') [fw=round(W_`dat')]
				matrix P[a,1] = -_b[ln_nw]
				matrix P[a,2] = sc_thres
				matrix P[a,3] = e(r2)
				matrix P[a,4] = e(r2_a)
				matrix P[a,5] = e(N)
				scalar a = a + 1
			di in red "+++ threshold: `= sc_thres', alpha: `=_b[ln_nw]'"

		}
		
		* transform matrix to Stata-file and process
		matrix list P
		drop _all
		svmat double P
		matrix drop P

		* rename and sort
		ren (P1 P2 P3 P4 P5) (alpha thres r2 adjR2 N)
		gen _imp_ = `imp'
		gsort -thres
		drop if thres == .
		
		save "${outpath}alpha_threshold_`dat'_`imp'.dta", replace

	}

	****************************************************************************
	* 10.2 Graph: Threshold-Alpha-Combination with all 5 nw in one Graph
	*	   (see Figure 4.4)
	****************************************************************************

	use "${outpath}alpha_threshold_`dat'_5.dta", clear
	
	forval i=4(-1)1 {
		append using "${outpath}alpha_threshold_`dat'_`i'.dta"
		append using "${outpath}alpha_threshold_`dat'_`i'.dta"
		append using "${outpath}alpha_threshold_`dat'_`i'.dta"
		append using "${outpath}alpha_threshold_`dat'_`i'.dta"
	}

	gen thres_mio = thres / 10^6

	twoway                      ///
		(line  alpha thres_mio if _imp_==1, sort(_imp_ thres_mio)) ///
		(line  alpha thres_mio if _imp_==2, sort(_imp_ thres_mio)) ///
		(line  alpha thres_mio if _imp_==3, sort(_imp_ thres_mio)) ///
		(line  alpha thres_mio if _imp_==4, sort(_imp_ thres_mio)) ///
		(line  alpha thres_mio if _imp_==5, sort(_imp_ thres_mio)), ///
		legend(row(1) order(1 "{it:nw}{sub:{it:1}}" 2 "{it:nw}{sub:{it:2}}" 3 "{it:nw}{sub:{it:3}}" 4 "{it:nw}{sub:{it:4}}" 5 "{it:nw}{sub:{it:5}}")) ///
		ytitle("Pareto's Alpha") xlab(, grid) xtitle("Threshold in mio. Euro")	///
		scheme(s2mono) graphregion(color(white)) ylabel(, angle(horizontal)) 	///
		title("Pareto-Threshold-Combinations") subtitle("``dat'_src'")			///
		note("Source: ``dat'_src'." "Note: data bases on all `=m' net wealth variables. We ran OLS every `=steps' Euros" "``dat'_w_info'")
	
	graph export "${graphs}10_alpha_thres_combi_`dat'_nw_combined.pdf", replace

}

***

