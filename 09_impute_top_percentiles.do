*** 9. Impute Top Percentiles ***

/*		
Here we impute the Top-Tails of the SOEP above p95, p99 with a) alpha_pt and for
robustness check with b) alpha_sp. Using the imputed values (for a) nw_hat, for 
b) nw_hat2), we then calculate the top percentiles and the total sum of wealth 
for the top 5%, top 1% and top 0.1%. Lastly, we use Van der Wijk's Law to 
calculate the average wealth and the total sum of wealth of those top %.

Impute Top-Tails:	nw_hat = ((1-F)^(-1/alpha)) * threshold

Note:
The above formula is valid for a distribution that is indeed Pareto-distributed 
accross all observations. However, we assume that the Pareto-Distribution is only
valid for obs. with a nw larger than the thresholds p95/p99. With this in mind, 
we construct a "synthetic" CDF (F2) by expanding the span where we assume the 
data to be pareto distributed (ie from 95th to 100th and from 99th to 100th pct).
For that, we use the formula: 

	F2 = (F - p) * ( 1 / (1 - p)), where p = {0.95; 0.99}
	
Had we not take this measure, our estimates with suffer a extreme jump at the 
threshold value, since we would be estimating from the 95th percentile onwards 
based on the original dataset, which is only pareto distributed after this pct. 
We would be, therefore, estimating, say, the 95th pct of the 99th pct. Further,
our estimated top-tails would be implausibly high.
*/



********************************************************************************
*	OUTLINE
*	9.1. Impute Top-Tails of SOEP and calculate some Goodness of Fit Measures
*	9.2. Calculate Top Percentiles and Total Sum based on imputed Top-Tails
*	9.3. Calculate the Mean of the Scalars over all five Implicates
*	9.4. Estimating the Mean and Total Sum above threshold with Van der Wijk's Law
******************************************************************************** 


use "${outpath}soep_pretest_2.dta", clear

scalar m=5

********************************************************************************
*
* 9.1. Impute Top-Tails of SOEP and calculate some Goodness of Fit Measures
* 
********************************************************************************


forvalue imp=1(1)`=m' {

	foreach pct in p95 p99 {
		if "`pct'" == "p95" {
			local p = 0.95
		}
		if "`pct'" == "p99" {
			local p = 0.99
		}
		gen F2_sp_`imp'_`pct' = (F_sp_`imp'_ - `p') * ( 1 / (1 - `p'))

		foreach citype in cilo cihi {
	
			* display alphas
			di in red "sc_`citype'_`pct'_pt = `=sc_`citype'_`pct'_pt'"
	
		**** 9.1a) with Alpha Pretest: nw_hat
			gen nw_hat_`imp'_`pct'_`citype' 		= _`imp'_nw 	if D_pt == 0 & _`imp'_nw <= sc_thres_`imp'_`pct'
			gen D_nw_hat_`imp'_`pct'_`citype' 		= 1 			if D_pt == 0
			replace D_nw_hat_`imp'_`pct'_`citype' 	= 0 			if D_pt == 0 & _`imp'_nw <= sc_thres_`imp'_`pct'	
			//replace nw_hat_`imp'_`pct'_`citype' 	= ((1-F_sp_`imp'_)^(-1/sc_`citype'_`pct'_pt))*sc_thres_`imp'_`pct' if D_pt == 0 & _`imp'_nw > sc_thres_`imp'_`pct'
			replace nw_hat_`imp'_`pct'_`citype' 	= ((1-F2_sp_`imp'_`pct')^(-1/sc_`citype'_`pct'_pt))*sc_thres_`imp'_`pct' if D_pt == 0 & _`imp'_nw > sc_thres_`imp'_`pct'

			* RMSE (rmse)
			gen trash = (_`imp'_nw - nw_hat_`imp'_`pct'_`citype')^2	if D_nw_hat_`imp'_`pct'_`citype' == 1
			qui sum trash [fw=round(W_sp)] 	if D_nw_hat_`imp'_`pct'_`citype' == 1
			scalar sc_`citype'rmse_`imp'_`pct' = round(sqrt((1/r(N))*r(sum)), .01)
			di in red "sc_`citype'rmse_`imp'_`pct' = `=sc_`citype'rmse_`imp'_`pct''"
			drop trash
			
			* RRMSE (rrmse) (values close to 0: good imp.)
			qui sum _`imp'_nw  [fw=round(W_sp)]  if D_nw_hat_`imp'_`pct'_`citype' == 1
			scalar sc_`citype'rrmse_`imp'_`pct' = round((sc_`citype'rmse_`imp'_`pct' / r(sum))*100, .001)
			di in red "sc_`citype'rrmse_`imp'_`pct' = `=sc_`citype'rrmse_`imp'_`pct''"

			
		**** 9.1b) with Alpha SOEP: nw_hat2
			gen nw_hat2_`imp'_`pct'_`citype' = _`imp'_nw if D_pt == 0 & _`imp'_nw <= sc_thres_`imp'_`pct'
			gen D_nw_hat2_`imp'_`pct'_`citype' = 1 if D_pt == 0
			replace D_nw_hat2_`imp'_`pct'_`citype' = 0 if D_pt == 0 & _`imp'_nw <= sc_thres_`imp'_`pct'	
			//replace nw_hat2_`imp'_`pct'_`citype' = ((1-F_sp_`imp'_)^(-1/sc_`citype'_`pct'_sp))*sc_thres_`imp'_`pct' if D_pt == 0 & _`imp'_nw > sc_thres_`imp'_`pct'
			replace nw_hat2_`imp'_`pct'_`citype' = ((1-F2_sp_`imp'_`pct')^(-1/sc_`citype'_`pct'_sp))*sc_thres_`imp'_`pct' if D_pt == 0 & _`imp'_nw > sc_thres_`imp'_`pct'
			
			* RMSE (rmse2)
			gen trash = (_`imp'_nw - nw_hat2_`imp'_`pct'_`citype')^2  if D_nw_hat2_`imp'_`pct'_`citype' == 1
			qui sum trash [fw=round(W_sp)] if D_nw_hat2_`imp'_`pct'_`citype' == 1
			scalar sc_`citype'rmse2_`imp'_`pct' = round(sqrt((1/r(N))*r(sum)), .01)
			di in red "sc_`citype'rmse2_`imp'_`pct' = `=sc_`citype'rmse2_`imp'_`pct''"
			drop trash
			
			* RRMSE (rrmse2) (values close to 0: good imp.)
			qui sum _`imp'_nw  [fw=round(W_sp)]  if D_nw_hat2_`imp'_`pct'_`citype' == 1
			scalar sc_`citype'rrmse2_`imp'_`pct' = round((sc_`citype'rmse2_`imp'_`pct' / r(sum))*100, .001)
			di in red "sc_`citype'rrmse2_`imp'_`pct' = `=sc_`citype'rrmse2_`imp'_`pct''"
			
		}
	}
	
	****************************************************************************
	*
	* 9.2. Calculate Top Percentiles and Total Sum based on imputed Top-Tails
	*
	****************************************************************************
	
	
	* Get the values of selected percentiles of original SOEP (needed below):
	_pctile _`imp'_nw [fw = round(W_sp)] if D_pt == 0, percentiles(95 99 99.9)
	scalar sc_pcttop50_sp_`imp' = r(r1)
	scalar sc_pcttop10_sp_`imp' = r(r2)
	scalar sc_pcttop01_sp_`imp' = r(r3)

	* Total Sum Wealth (tw) above threshold, original SOEP:
	foreach top in top50 top10 top01 {	
		sum _`imp'_nw [fw = round(W_sp)] if D_pt == 0 & (_`imp'_nw > sc_pct`top'_sp_`imp')
		scalar sc_tw`top'_sp_`imp' = r(sum)
		di in red "sc_tw`top'_sp_`imp' = `=sc_tw`top'_sp_`imp''"
	}

	
**** 9.2a) with Alpha Pretest, (p95, p99)
	
	** Using nw_hat from 9.1a) estimate Top Percentiles (I) / Total Sum Wealth (I)
	** with Alpha Pretest, p95:
	
	* (I) Estimate Top Percentiles
	foreach citype in cilo cihi {		
		_pctile nw_hat_`imp'_p95_`citype' [fw = round(W_sp)] if D_pt == 0, percentiles(95 99 99.9)
		scalar sc_pcttop50_`imp'_p95_`citype' = r(r1)
		scalar sc_pcttop10_`imp'_p95_`citype' = r(r2)
		scalar sc_pcttop01_`imp'_p95_`citype' = r(r3)
		di in red "sc_pcttop50_`imp'_p95_`citype' `=sc_pcttop50_`imp'_p95_`citype''"
		di in red "sc_pcttop10_`imp'_p95_`citype' `=sc_pcttop10_`imp'_p95_`citype''"
		di in red "sc_pcttop01_`imp'_p95_`citype' `=sc_pcttop01_`imp'_p95_`citype''"

		* (II) Sum nw_hat from 9.1a) with threshold p95 to get r(sum) = Total Sum Wealth (tw):
		foreach top in top50 top10 top01 {		
			sum nw_hat_`imp'_p95_`citype' [fw = round(W_sp)] if D_pt == 0 & (nw_hat_`imp'_p95_`citype' > sc_pct`top'_`imp'_p95_`citype')
			scalar sc_tw`top'_`imp'_p95_`citype' = r(sum)
			di in red "sc_tw`top'_`imp'_p95_`citype' = `=sc_tw`top'_`imp'_p95_`citype''"
		}
	}

	** Using nw_hat from 9.1a) estimate Top Percentiles (I) / Total Sum Wealth (I)
	** with Alpha Pretest, p99:
	
	* (I) Estimate Top Percentiles
	foreach citype in cilo cihi {
	
		_pctile nw_hat_`imp'_p99_`citype' [fw = round(W_sp)] if D_pt == 0, percentiles(99 99.9)
		scalar sc_pcttop10_`imp'_p99_`citype' = r(r1)
		scalar sc_pcttop01_`imp'_p99_`citype' = r(r2)
		di in red "sc_pcttop10_`imp'_p99_`citype' `=sc_pcttop10_`imp'_p99_`citype''"
		di in red "sc_pcttop01_`imp'_p99_`citype' `=sc_pcttop01_`imp'_p99_`citype''"
	
		* (II) Sum nw_hat from 9.1b) with threshold p99 to get r(sum) = Total Sum Wealth (tw):
		foreach top in top10 top01 {
			sum nw_hat_`imp'_p99_`citype' [fw = round(W_sp)] if D_pt == 0 & (nw_hat_`imp'_p99_`citype' > sc_pct`top'_`imp'_p99_`citype')
			scalar sc_tw`top'_`imp'_p99_`citype' = r(sum)
			di in red "sc_tw`top'_`imp'_p99_`citype' = `=sc_tw`top'_`imp'_p99_`citype''"
		}
	}

**** 9.2b) with Alpha SOEP, (p95, p99)
	* Use alpha with threshold p95 - percentiles:
	foreach citype in cilo cihi {
		_pctile nw_hat2_`imp'_p95_`citype' [fw = round(W_sp)] if D_pt == 0, percentiles(95 99 99.9)
		scalar sc_pct2top50_`imp'_p95_`citype' = r(r1)
		scalar sc_pct2top10_`imp'_p95_`citype' = r(r2)
		scalar sc_pct2top01_`imp'_p95_`citype' = r(r3)
		di in red "sc_pct2top50_`imp'_p95_`citype' = `=sc_pct2top50_`imp'_p95_`citype''"
		di in red "sc_pct2top10_`imp'_p95_`citype' = `=sc_pct2top10_`imp'_p95_`citype''"
		di in red "sc_pct2top01_`imp'_p95_`citype' = `=sc_pct2top01_`imp'_p95_`citype''"

		* Use alpha with threshold p95 - total wealth:
		foreach top in top50 top10 top01 {
			sum nw_hat2_`imp'_p95_`citype' [fw = round(W_sp)] if D_pt == 0 & (nw_hat2_`imp'_p95_`citype' > sc_pct2`top'_`imp'_p95_`citype')
			scalar sc_tw2`top'_`imp'_p95_`citype' = r(sum)
			di in red "sc_tw2`top'_`imp'_p95_`citype' = `=sc_tw2`top'_`imp'_p95_`citype''"
		}	
	}

	* Use alpha with threshold p99 - percentiles:
	foreach citype in cilo cihi {
		_pctile nw_hat2_`imp'_p99_`citype' [fw = round(W_sp)] if D_pt == 0, percentiles(99 99.9)
		scalar sc_pct2top10_`imp'_p99_`citype' = r(r1)
		scalar sc_pct2top01_`imp'_p99_`citype' = r(r2)
		di in red "sc_pct2top10_`imp'_p99_`citype' = `=sc_pct2top10_`imp'_p99_`citype''"
		di in red "sc_pct2top01_`imp'_p99_`citype' = `=sc_pct2top01_`imp'_p99_`citype''"

		
		* Use alpha with threshold p99 - total wealth:
		foreach top in top10 top01 {
			sum nw_hat2_`imp'_p99_`citype' [fw = round(W_sp)] if D_pt == 0 & (nw_hat2_`imp'_p99_`citype' > sc_pct2`top'_`imp'_p99_`citype')
			scalar sc_tw2`top'_`imp'_p99_`citype' = r(sum)
			di in red "sc_tw2`top'_`imp'_p99_`citype' = `=sc_tw2`top'_`imp'_p99_`citype''"
		}
	}
}


** check nw_hat, 1 example comparison
	if "$show_all" == "TRUE" {
		* check nw_hat based on Pretest alpha
	sum nw_hat_1_p95_cihi nw_hat_1_p95_cilo _1_nw  if D_nw_hat_1_p95_cihi ==1
	sum nw_hat_1_p99_cihi nw_hat_1_p99_cilo _1_nw  if D_nw_hat_1_p99_cihi ==1
		* check nw_hat2 based on SOEP alpha
	sum nw_hat2_1_p95_cihi nw_hat2_1_p95_cilo _1_nw  if D_nw_hat2_1_p95_cihi ==1
	sum nw_hat2_1_p99_cihi nw_hat2_1_p99_cilo _1_nw  if D_nw_hat2_1_p99_cihi ==1
}

********************************************************************************
*
*	9.3 Calculate the Mean of the Scalars over all five Implicates
*
********************************************************************************


***** Average original SOEP values across all imp.: (I) Top Percentiles and (II) Total Sum Wealths
foreach top in top50 top10 top01 {
	* (I) Top Percentiles (pct)
	scalar sc_pct`top'_sp = (sc_pct`top'_sp_1 + sc_pct`top'_sp_2 + sc_pct`top'_sp_3 + sc_pct`top'_sp_4 + sc_pct`top'_sp_5)/m 
	* (II) Total Sum Wealths (tw)
	scalar sc_tw`top'_sp = (sc_tw`top'_sp_1 + sc_tw`top'_sp_2 + sc_tw`top'_sp_3 + sc_tw`top'_sp_4 + sc_tw`top'_sp_5)/m 
}
	
**** 9.3a) with Alpha Pretest
foreach top in top50 top10 top01 {
	** Average of Imputations across all imp.: (I) Top Percentiles  and (II) Total Sum Wealths
	* (Pretest alpha, threshold p95)
	foreach citype in cilo cihi {	
		* (I) top percentiles (pct) in mio.
		scalar sc_pct`top'_p95_`citype'	= (sc_pct`top'_1_p95_`citype' + sc_pct`top'_2_p95_`citype' + sc_pct`top'_3_p95_`citype' + sc_pct`top'_4_p95_`citype' + sc_pct`top'_5_p95_`citype')/m 
		scalar sc_pct`top'_p95_`citype'	= sc_pct`top'_p95_`citype' / (10^6)
		di in red "sc_pct`top'_p95_`citype' in mio. = `=sc_pct`top'_p95_`citype''"		
		* (II) total wealth (tw) in billions
		scalar sc_tw`top'_p95_`citype'	= (sc_tw`top'_1_p95_`citype' + sc_tw`top'_2_p95_`citype' + sc_tw`top'_3_p95_`citype' + sc_tw`top'_4_p95_`citype' + sc_tw`top'_5_p95_`citype')/m 
		scalar sc_tw`top'_p95_`citype'	= sc_tw`top'_p95_`citype' / (10^9)
		di in red "sc_tw`top'_p95_`citype' in mio. = `=sc_tw`top'_p95_`citype''"
	}
}

	** Average of Imputations across all imp.: (I) Top Percentiles  and (II) Total Sum Wealths
	* (Pretest alpha, threshold p99)
	foreach top in top10 top01 {
		foreach citype in cilo cihi {
		* (I) top percentiles (pct) in mio.
		scalar sc_pct`top'_p99_`citype'	= (sc_pct`top'_1_p99_`citype' + sc_pct`top'_2_p99_`citype' + sc_pct`top'_3_p99_`citype' + sc_pct`top'_4_p99_`citype' + sc_pct`top'_5_p99_`citype')/m 
		scalar sc_pct`top'_p99_`citype' = sc_pct`top'_p99_`citype' / (10^6)
		di in red "sc_pct`top'_p99_`citype' in mio. = `=sc_pct`top'_p99_`citype''"
		* (II) total wealth (tw) in billions
		scalar sc_tw`top'_p99_`citype'	= (sc_tw`top'_1_p99_`citype' + sc_tw`top'_2_p99_`citype' + sc_tw`top'_3_p99_`citype' + sc_tw`top'_4_p99_`citype' + sc_tw`top'_5_p99_`citype')/m 
		scalar sc_tw`top'_p99_`citype'	= sc_tw`top'_p99_`citype' / (10^9)
		di in red "sc_tw`top'_p99_`citype' in mio. = `=sc_tw`top'_p99_`citype''"	
	}	
}

**** 9.3b) with Alpha SOEP
foreach top in top50 top10 top01 {
	** Average of Imputations across all imp.: (I) Top Percentiles  and (II) Total Sum Wealths
	* (SOEP alpha, threshold p95)
	foreach citype in cilo cihi {	
		* (I) top percentiles (pct2) in mio.
		scalar sc_pct2`top'_p95_`citype'	= (sc_pct2`top'_1_p95_`citype' + sc_pct2`top'_2_p95_`citype' + sc_pct2`top'_3_p95_`citype' + sc_pct2`top'_4_p95_`citype' + sc_pct2`top'_5_p95_`citype')/m 
		scalar sc_pct2`top'_p95_`citype'	= sc_pct2`top'_p95_`citype' / (10^6)
		di in red "sc_pct2`top'_p95_`citype' in mio. = `=sc_pct2`top'_p95_`citype''"		
		* (II) total wealth (tw2) in billions
		scalar sc_tw2`top'_p95_`citype'	= (sc_tw2`top'_1_p95_`citype' + sc_tw2`top'_2_p95_`citype' + sc_tw2`top'_3_p95_`citype' + sc_tw2`top'_4_p95_`citype' + sc_tw2`top'_5_p95_`citype')/m 
		scalar sc_tw2`top'_p95_`citype'	= sc_tw2`top'_p95_`citype' / (10^9)
		di in red "sc_tw2`top'_p95_`citype' in mio. = `=sc_tw2`top'_p95_`citype''"
	}
}

	** Average of Imputations across all imp.: (I) Top Percentiles  and (II) Total Sum Wealths
	* (SOEP alpha, threshold p99)
	foreach top in top10 top01 {
		foreach citype in cilo cihi {
		* (I) top percentiles (pct2) in mio.
		scalar sc_pct2`top'_p99_`citype' = (sc_pct2`top'_1_p99_`citype' + sc_pct2`top'_2_p99_`citype' + sc_pct2`top'_3_p99_`citype' + sc_pct2`top'_4_p99_`citype' + sc_pct2`top'_5_p99_`citype')/m 
		scalar sc_pct2`top'_p99_`citype' = sc_pct2`top'_p99_`citype' / (10^6)
		di in red "sc_pct2`top'_p99_`citype' in mio. = `=sc_pct2`top'_p99_`citype''"
		* (II) total wealth (tw2) in billions
		scalar sc_tw2`top'_p99_`citype'	= (sc_tw2`top'_1_p99_`citype' + sc_tw2`top'_2_p99_`citype' + sc_tw2`top'_3_p99_`citype' + sc_tw2`top'_4_p99_`citype' + sc_tw2`top'_5_p99_`citype')/m 
		scalar sc_tw2`top'_p99_`citype'	= sc_tw2`top'_p99_`citype' / (10^9)
		di in red "sc_tw2`top'_p99_`citype' in mio. = `=sc_tw2`top'_p99_`citype''"	
	}	
}


********************************************************************************
*
* 9.4. Estimating the Mean and Total Sum above threshold with Van der Wijk's Law
*
********************************************************************************

/*
Note:
Calculate the Mean above Threshold according to Van der Wijk's Law (vdw) with 
the Pretest Alpha and Thresholds

First, we calculated vdw-Measures: (mean, total sum wealth)
	twvdw = total sum wealth Van der Wijk's Law
	avgvdw = mean Van der Wijk's Law
		
However, it turns out that the Pretest delivers alphas which are smaller 
than 1. Thus, the quotient (alpha/(alpha-1)) becomes negative and we end 
up in negative vdw-measures. Finally, we leave this section out in our paper.
*/


* Calculate the Mean above Threshold according to Van der Wijk's Law (vdw)
*	with the Pretest Alpha and Threshold defined at p95

forvalue imp=1(1)`=m' {

	foreach top in top50 top10 top01 {
		
		foreach citype in cilo cihi {
		
			** p95 **
			* mean above `top' pct (avg)
			scalar sc_avgvdw`top'_`imp'_p95_`citype'	= (sc_`citype'_p95_pt / (sc_`citype'_p95_pt - 1)) * sc_pct`top'_`imp'_p95_`citype'
			di in red "sc_avgvdw`top'_`imp'_p95_`citype' = `=sc_avgvdw`top'_`imp'_p95_`citype''"
			
			* total sum above `top' pct (twvdw)
				if "`top'" == "top50" {
				* Top 5%
					scalar sc_twvdw`top'_`imp'_p95_`citype'	= 0.05 * (sc_`citype'_p95_pt / (sc_`citype'_p95_pt - 1)) * sc_pct`top'_`imp'_p95_`citype' * sc_N_sp
					di in red "sc_twvdw`top'_`imp'_p95_`citype' = `=sc_twvdw`top'_`imp'_p95_`citype''"
				}
				* Top 1%
				if "`top'" == "top10" {
					scalar sc_twvdw`top'_`imp'_p95_`citype'	= 0.01 * (sc_`citype'_p95_pt / (sc_`citype'_p95_pt - 1)) * sc_pct`top'_`imp'_p95_`citype' * sc_N_sp
					di in red "sc_twvdw`top'_`imp'_p95_`citype' = `=sc_twvdw`top'_`imp'_p95_`citype''"
				}
				* Top 0.1%
				if "`top'" == "top01" {
					scalar sc_twvdw`top'_`imp'_p95_`citype'	= 0.001 * (sc_`citype'_p95_pt / (sc_`citype'_p95_pt - 1)) * sc_pct`top'_`imp'_p95_`citype' * sc_N_sp
					di in red "sc_twvdw`top'_`imp'_p95_`citype' = `=sc_twvdw`top'_`imp'_p95_`citype''"
				}
		}
		
	}
}

* Calculate the Mean above Threshold according to Van der Wijk's Law (vdw)
*	with the Pretest Alpha and Threshold defined at p99

forvalue imp=1(1)`=m' {

	foreach top in top10 top01 {
		
		foreach citype in cilo cihi {
		
			** p99 **
			* mean above `top' pct (avg)
			scalar sc_avgvdw`top'_`imp'_p99_`citype'	= (sc_`citype'_p99_pt / (sc_`citype'_p99_pt - 1)) * sc_pct`top'_`imp'_p99_`citype'
			di in red "sc_avgvdw`top'_`imp'_p99_`citype' = `=sc_avgvdw`top'_`imp'_p99_`citype''"
			
			* total sum above `top' pct (twvdw)
				* Top 1%
				if "`top'" == "top10" {
					scalar sc_twvdw`top'_`imp'_p99_`citype'	= 0.01 * (sc_`citype'_p99_pt / (sc_`citype'_p99_pt - 1)) * sc_pct`top'_`imp'_p99_`citype' * sc_N_sp
					di in red "sc_twvdw`top'_`imp'_p99_`citype' = `=sc_twvdw`top'_`imp'_p99_`citype''"
				}
				* Top 0.1%
				if "`top'" == "top01" {
					scalar sc_twvdw`top'_`imp'_p99_`citype'	= 0.001 * (sc_`citype'_p99_pt / (sc_`citype'_p99_pt - 1)) * sc_pct`top'_`imp'_p99_`citype' * sc_N_sp
					di in red "sc_twvdw`top'_`imp'_p99_`citype' = `=sc_twvdw`top'_`imp'_p99_`citype''"
				}
		}
		
	}
}

* Calculate the Average of above Scalars

foreach top in top50 top10 top01 {
	
	foreach citype in cilo cihi {
	
		** p95 **
		scalar sc_avgvdw`top'_p95_`citype'	= (sc_avgvdw`top'_1_p95_`citype' + sc_avgvdw`top'_2_p95_`citype' + sc_avgvdw`top'_3_p95_`citype' + sc_avgvdw`top'_4_p95_`citype' + sc_avgvdw`top'_5_p95_`citype')/m 
		scalar sc_avgvdw`top'_p95_`citype' = round(sc_avgvdw`top'_p95_`citype' / 1000000, .001)
		di in red "sc_avgvdw`top'_p95_`citype' = `=sc_avgvdw`top'_p95_`citype''"
		** p95 **
		scalar sc_twvdw`top'_p95_`citype'	= (sc_twvdw`top'_1_p95_`citype' + sc_twvdw`top'_2_p95_`citype' + sc_twvdw`top'_3_p95_`citype' + sc_twvdw`top'_4_p95_`citype' + sc_twvdw`top'_5_p95_`citype')/m 
		scalar sc_twvdw`top'_p95_`citype' = round(sc_twvdw`top'_p95_`citype' / 1000000, .001)
		di in red "sc_twvdw`top'_p95_`citype' = `=sc_twvdw`top'_p95_`citype''"
	}
}

foreach top in top10 top01 {
	
	foreach citype in cilo cihi {
	
		** p99 **
		scalar sc_avgvdw`top'_p99_`citype'	= (sc_avgvdw`top'_1_p99_`citype' + sc_avgvdw`top'_2_p99_`citype' + sc_avgvdw`top'_3_p99_`citype' + sc_avgvdw`top'_4_p99_`citype' + sc_avgvdw`top'_5_p99_`citype')/m 
		scalar sc_avgvdw`top'_p99_`citype' = round(sc_avgvdw`top'_p99_`citype' / 1000000, .001)
		di in red "sc_avgvdw`top'_p99_`citype' = `=sc_avgvdw`top'_p99_`citype''"
		** p99 **
		scalar sc_twvdw`top'_p99_`citype'	= (sc_twvdw`top'_1_p99_`citype' + sc_twvdw`top'_2_p99_`citype' + sc_twvdw`top'_3_p99_`citype' + sc_twvdw`top'_4_p99_`citype' + sc_twvdw`top'_5_p99_`citype')/m 
		scalar sc_twvdw`top'_p99_`citype' = round(sc_twvdw`top'_p99_`citype' / 1000000, .001)
		di in red "sc_twvdw`top'_p99_`citype' = `=sc_twvdw`top'_p99_`citype''"
	}
}
	
save "${outpath}soep_pretest_3.dta", replace

***
