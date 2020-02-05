*** 7. Fitting Pareto ***

/*
Here, the Pareto's coefficients are estimated on the basis of the SOEP and 
the Pretest seperately according to section 3.2 in the paper. 
Since the nw variable contains 5 values per HH, we estimate 5 alphas per dataset.
We use as thresholds the 95th and 99th percentile (p95 and p99) of the SOEP. 

In total, we estimate 20 alphas: 
SOEP: 		5x for threshold at p95, 5x for threshold at p99
Pretest:	5x for threshold at p95, 5x for threshold at p99
	
For each estimation, we save following scalars:
- Pareto's coefficients (alphas)
- standard deviation (sd)
- total number of observations (N)
- constant (cons)
- R-squared (r2)
- resid. sum of sq. (rss)

After the estimation, we follow Rubin (1987) by adjusting for the increased 
uncertainty about the true parameter of alpha due to multiple imputation. This 
enables us to compose confidence intervals and to test the equality of the estimated 
parameters in the presence of multiple imputed data.
*/

********************************************************************************
*	OUTLINE
*	7.1.	Estimating Pareto's coefficients
*	7.2.	Calculating CI of Pareto's coefficients according to Rubin (1987)
********************************************************************************


use "${outpath}soep_pretest_2.dta", clear

* m: number of imputations
scalar m = 5

* define threshold levels as percentiles of SOEP
local pctile p95 p99

* define percentile of alpha CI
scalar sc_alphaCI = .975

* global: display all scalars below
global show_all "TRUE"

********************************************************************************
*		
*	7.1. Estimating Pareto's Coefficients
*		
********************************************************************************

forval imp=1(1)`=m' {

	* Estimate Threshold Values (only for SOEP)
	qui sum _`imp'_nw if D_pt==0 [fw=round(W_sp)], d
	foreach pc of local pctile {
		scalar sc_thres_`imp'_`pc' = r(`pc')
	}
	
	* Estimate Pareto's Coefficients (SOEP and Pretest, with threshold from SOEP)
	foreach dat in sp pt {

		foreach pct of local pctile {
	
			if "`dat'" == "sp" {
				local D_type 0
			}
			if "`dat'" == "pt" {
				local D_type 1
			}
			
				* run regression
				reg lnP_`dat'_`imp'_ ln_nw_`imp'_ if (_`imp'_nw > sc_thres_`imp'_`pct' & D_pt==`D_type') [fw=round(W_`dat')]
					
					* name regression
					estimates title: reg lnP_`dat'_`imp'_ ~ ln_nw_`imp'_ (thres: `pct')
					
					* save regression
					est store reg_`imp'_`pct'_`dat'

					* save scalars
					scalar sc_alpha_`imp'_`pct'_`dat' = -_b[ln_nw] /* -alpha */
					scalar sc_sd_`imp'_`pct'_`dat' = _se[ln_nw]
					scalar sc_N_`imp'_`pct'_`dat' = e(N)
					scalar sc_rss_`imp'_`pct'_`dat' = e(rss) /* resid. sum of sq. needed for CI calc. */
					scalar sc_r2_`imp'_`pct'_`dat' = e(r2)
					scalar sc_cons_`imp'_`pct'_`dat' = _b[_cons]

				* display results
				if "$show_all" == "TRUE" {
					if "`dat'" == "sp" {
						local data "SOEP"
					}
					else {
						local data "Pretest"
					}
					
					di in red "+++++++++++++++++++++++"
					di "dataset:		`data'"
					di "threshold pctile: 	`pct'"
					di "threshold value: 	sc_thres_`imp'_`pct' = `=sc_thres_`imp'_`pct''"
					di "alpha:		sc_alpha_`imp'_`pct'_`dat' = `=sc_alpha_`imp'_`pct'_`dat''"
					di "sd:    		sc_sd_`imp'_`pct'_`dat' = `=sc_sd_`imp'_`pct'_`dat''"	
					di "N: 			sc_N_`imp'_`pct'_`dat' = `=sc_N_`imp'_`pct'_`dat''"
					di "sse:		sc_sse_`imp'_`pct'_`dat' = `=sc_rss_`imp'_`pct'_`dat''"
					di "r2:			sc_r2_`imp'_`pct'_`dat' = `=sc_r2_`imp'_`pct'_`dat''"
				}
		}
	}
}

********************************************************************************
*		
*	7.2.	Calculating CI of alpha according to Rubin (1987)
*		
********************************************************************************

/* 
NOTE:
One estimate of alpha for each dataset and treshold is desirable. 
We follow Rubin's (1987) instructions of dealing with multiply imputed data. By
this, we derive one alpha per dataset and threshold. In a last step, we calculate 
the CIs for each of the four alphas.

	Content:
		1. Estimating Alpha (Average of each five alphas)
		2. AWV: Average-Within-Variance
			2.1. Weighted mean of logged NW (needed for denominator of variance of alpha)
			2.2. Denominator of variance of alpha
			2.3. Variance of Alpha
			2.3. WV: Within-Variance
			2.4. AWV: Average-Within-Variance
		3. BV: Between Variance
		4. TV: Total Variance 
		5. Confidence Intervals
*/

foreach dat in sp pt {

	if "`dat'" == "sp" {
		local dataset "D_pt==0"
	}
	
	if "`dat'" == "pt" {
		local dataset "D_pt==1"
	}

/*
NOTE:
We refrain from looping over percentiles to increase the readability of the code.
*/

** 1. Estimate of Alpha (Mean of each 5 alphas per dataset and percentile) 
*	  (see formula 3.7)

		* threshold: p95
	qui scalar sc_avgalpha_p95_`dat' = 1/m*(sc_alpha_1_p95_`dat' + sc_alpha_2_p95_`dat' + sc_alpha_3_p95_`dat' + sc_alpha_4_p95_`dat' + sc_alpha_5_p95_`dat')
	di " avg. alpha (p95, `dat'): `=sc_avgalpha_p95_`dat''"

		* threshold: p99
	qui scalar sc_avgalpha_p99_`dat' = 1/m*(sc_alpha_1_p99_`dat' + sc_alpha_2_p99_`dat' + sc_alpha_3_p99_`dat' + sc_alpha_4_p99_`dat' + sc_alpha_5_p99_`dat')
	di " avg. alpha (p99, `dat'): `=sc_avgalpha_p99_`dat''"
	
	
** 2. AWV: Average-Within-Variance (see formula 3.15)

	forval imp=1(1)`=m' {
	
		** 2.1. Generate weighted mean of logged net wealth 
		**	   (needed for the denominator of the Variance of alpha) 
		**	   (see formula 3.8 and 3.11)
		
			* threshold: p95
		qui sum ln_nw_`imp'_ [fw=round(W_`dat')] if _`imp'_nw > sc_thres_`imp'_p95 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_sumwlnnw_`dat'_`imp'_p95 = r(sum)
		qui sum ln_nw_`imp'_ if _`imp'_nw > sc_thres_`imp'_p95 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_n_`imp'_p95_`dat' = r(N)
		di " sc_sumwlnnw_`dat'_`imp'_p95 = `=sc_sumwlnnw_`dat'_`imp'_p95'"
		di " sc_n_`imp'_p95_`dat' = `=sc_n_`imp'_p95_`dat''"
		
		qui sum W_`dat' if _`imp'_nw > sc_thres_`imp'_p95 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_wmeanlnnw_p95_`dat'_`imp' = sc_sumwlnnw_`dat'_`imp'_p95 / r(sum)
		di " sc_wmeanlnnw_p95_`dat'_`imp' = `=sc_wmeanlnnw_p95_`dat'_`imp''"

			* threshold: p99
		qui sum ln_nw_`imp'_ [fw=round(W_`dat')] if _`imp'_nw > sc_thres_`imp'_p99 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_sumwlnnw_`dat'_`imp'_p99 = r(sum)
		qui sum ln_nw_`imp'_ if _`imp'_nw > sc_thres_`imp'_p99 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_n_`imp'_p99_`dat' = r(N)
		di " sc_sumwlnnw_`dat'_`imp'_p99 = `=sc_sumwlnnw_`dat'_`imp'_p99'"
		di " sc_n_`imp'_p99_`dat' = `=sc_n_`imp'_p99_`dat''"

		qui sum W_`dat' if _`imp'_nw > sc_thres_`imp'_p99 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_wmeanlnnw_p99_`dat'_`imp' = sc_sumwlnnw_`dat'_`imp'_p99 / r(sum)
		di " sc_wmeanlnnw_p99_`dat'_`imp' = `=sc_wmeanlnnw_p99_`dat'_`imp''"

		
		** 2.2. Denominator of the variance of Alpha (see formula 3.9 and 3.12)
		**		(Sum of the squared residuals of the weighted log. net wealth 
		**		from their weighted means)
		
			* threshold: p95
		gen sq2devmeanlnnw_`dat'_`imp'_p95 = (W_`dat' * ln_nw_`imp'_ - sc_wmeanlnnw_p95_`dat'_`imp')^2
		qui sum sq2devmeanlnnw_`dat'_`imp'_p95 if _`imp'_nw > sc_thres_`imp'_p95 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_sumsq2devlnnw_p95_`dat'_`imp' = r(sum)
		di " sc_sumsq2devlnnw_p95_`dat'_`imp' = `=sc_sumsq2devlnnw_p95_`dat'_`imp''"

			* threshold: p99
		gen sq2devmeanlnnw_`dat'_`imp'_p99 = (W_`dat' * ln_nw_`imp'_ - sc_wmeanlnnw_p99_`dat'_`imp')^2
		qui sum sq2devmeanlnnw_`dat'_`imp'_p99 if _`imp'_nw > sc_thres_`imp'_p99 & `dataset' & ln_nw_`imp'_!=.
		qui scalar sc_sumsq2devlnnw_p99_`dat'_`imp' = r(sum)
		di " sc_sumsq2devlnnw_p99_`dat'_`imp' = `=sc_sumsq2devlnnw_p99_`dat'_`imp''"

		drop sq2devmeanlnnw_`dat'_`imp'_p95   sq2devmeanlnnw_`dat'_`imp'_p99
		
		
		** 2.3 	Variance of Alpha, according to the variance of the parameter alpha 
		**		in the Simple Linear Regression (SLR). (see formula 3.10 and 3.13)
		**		Let y=c0+alpha*x+epsilon be the SLR model function. Following, 
		**		Var(alpha) = Sum(e_i^2)/ sum(x_i-xbar)^2
		**		nominator Var(alpha): rss taken from the estimation above (see formula 3.9 and 3.12)
		**		denominator: calculated in section 2.2 in this do-file
		
			* threshold: p95
		qui scalar sc_v_`dat'_`imp'_p95 = (sc_rss_`imp'_p95_`dat' / sc_sumsq2devlnnw_p95_`dat'_`imp')
		di " sc_v_`dat'_`imp'_p95 = `=sc_v_`dat'_`imp'_p95'"	
		
			* threshold: p99
		qui scalar sc_v_`dat'_`imp'_p99 =(sc_rss_`imp'_p99_`dat' / sc_sumsq2devlnnw_p99_`dat'_`imp')
		di " sc_v_`dat'_`imp'_p99 = `=sc_v_`dat'_`imp'_p99'"					
		
		
		** 2.4. WV: Within-Variance (= Variance of Alpha / n) (see formula 3.14)
		
			* threshold: p95
		qui scalar sc_wv_`dat'_`imp'_p95 = (sc_v_`dat'_`imp'_p95 / sc_n_`imp'_p95_`dat')
		di " sc_wv_`dat'_`imp'_p95 = `=sc_wv_`dat'_`imp'_p95'"								

			* threshold: p99
		qui scalar sc_wv_`dat'_`imp'_p99 = (sc_v_`dat'_`imp'_p99 / sc_n_`imp'_p99_`dat')
		di " sc_wv_`dat'_`imp'_p99 = `=sc_wv_`dat'_`imp'_p99'"								

		di "++++++++"
	}


		** 2.5. AWV: Average Within-Variance (= 1/m * sum(WV)) (see formula 3.15)

			* threshold: p95
		qui scalar sc_awv_p95_`dat' = 1/m*((sc_wv_`dat'_1_p95 + sc_wv_`dat'_2_p95 + sc_wv_`dat'_3_p95 + sc_wv_`dat'_4_p95 + sc_wv_`dat'_5_p95))
		di " avg. within variance (p95, `dat'): `=sc_awv_p95_`dat''" 

			* threshold: p99
		qui scalar sc_awv_p99_`dat' = 1/m*((sc_wv_`dat'_1_p99 + sc_wv_`dat'_2_p99 + sc_wv_`dat'_3_p99 + sc_wv_`dat'_4_p99 + sc_wv_`dat'_5_p99))
		di " avg. within variance (p99, `dat'): `=sc_awv_p99_`dat''"


	** 3. BV: Between-Variance of Alphas (= sum((alpha_j - mean_alpha)^2)) (see formula 3.16)

		* threshold: p95
	qui scalar sc_bv_p95_`dat' = 0
		forval i=1(1)`=m' {
			qui scalar sc_bv_p95_`dat' = sc_bv_p95_`dat' + (sc_alpha_`i'_p95_`dat' - sc_avgalpha_p95_`dat')^2
		}
	di " between variance (p95, `dat'): `=sc_bv_p95_`dat''"

		* threshold: p99
	qui scalar sc_bv_p99_`dat' = 0
		forval i=1(1)`=m' {
			qui scalar sc_bv_p99_`dat' = sc_bv_p99_`dat' + (sc_alpha_`i'_p99_`dat' - sc_avgalpha_p99_`dat')^2
		}
	di " between variance (p99, `dat'): `=sc_bv_p99_`dat''"


	** 4. TV: total variance (= AWV + (1 + m^(-1)) * BV) (compare to formula 3.17)	

		* threshold: p99
	qui scalar sc_tv_p95_`dat' = sc_awv_p95_`dat' + (1+m^-1) * sc_bv_p95_`dat'
	di " total variance (p95, `dat'): `=sc_tv_p95_`dat''"

		* threshold: p99
	qui scalar sc_tv_p99_`dat' = sc_awv_p99_`dat' + (1+m^-1) * sc_bv_p99_`dat'
	di " total variance (p99, `dat'): `=sc_tv_p99_`dat''"

	
	** 5 Confidence Intervals: lower and upper bound (compare to formulas 3.18 and 3.19)

		* threshold: p95
	qui scalar sc_cilo_p95_`dat' = sc_avgalpha_p95_`dat' - invnorm(sc_alphaCI) * sqrt(sc_tv_p95_`dat')
	qui scalar sc_cihi_p95_`dat' = sc_avgalpha_p95_`dat' + invnorm(sc_alphaCI) * sqrt(sc_tv_p95_`dat')
	di " CI low (p95, `dat'): `=sc_cilo_p95_`dat''"
	di " CI high (p95, `dat'): `=sc_cihi_p95_`dat''"

		* threshold: p99
	qui scalar sc_cilo_p99_`dat' = sc_avgalpha_p99_`dat' - invnorm(sc_alphaCI) * sqrt(sc_tv_p99_`dat')
	qui scalar sc_cihi_p99_`dat' = sc_avgalpha_p99_`dat' + invnorm(sc_alphaCI) * sqrt(sc_tv_p99_`dat')
	di " CI low (p99, `dat'): `=sc_cilo_p99_`dat''"
	di " CI high (p99, `dat'): `=sc_cihi_p99_`dat''"

}

***
