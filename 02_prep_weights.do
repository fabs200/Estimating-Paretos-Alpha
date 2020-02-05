***	2. Weights ***

/*
For the SOEP HHs we use the HH weights as loaded in 01_prep_data.do. We do NOT 
have frequency weights for the Pretest dataset. However, due to the probability 
sampling design applied in the Pretest, we need to adjust for the different 
drawing probabilities by reweighting the interviewed HHs 
(see our detailed explanation in section 2.3 in our paper). 
*/

use "${outpath}soep_pretest_0.dta", clear

********************************************************************************
*		
*	2.1 Weights: SOEP (_sp)
*		
********************************************************************************

* check plausibility of weights
forval imp=1(1)5{
	sum _`imp'_nw [fw=round(W_sp,1)] if D_pt == 0
}

qui sum _1_nw [fw=round(W_sp,1)] if D_pt == 0
scalar sc_N_sp = r(N)

********************************************************************************
*		
*	2.2 Re-weighting: Pretest (_pt)
*		
********************************************************************************

gen W_pt = .
scalar sc_strata = 3

/* 	Since the drawn Pretest HHs have different response rates than the sampling 
probabilities, we apply the inverse of the sample share for each stratum. */

* actual shares of Pretest respondents
qui sum D_pt if D_pt==1
scalar sc_N_pt = r(N)
di in red "pretest n = `=sc_N_pt'"

forval i=1(1)3 {
	qui sum D_pt if schicht==`i'
	scalar sc_share_strat`i' = r(N)/sc_N_pt
	di in red "share stratum `i' = `=sc_share_strat`i''" 
}

* re-weighting according to their relative response rates
* Stratum 1
replace W_pt = (1/`=sc_share_strat1') if D_pt==1 & schicht==1
* Stratum 2
replace W_pt = (1/`=sc_share_strat2') if D_pt==1 & schicht==2
* Stratum 3
replace W_pt = (1/`=sc_share_strat3') if D_pt==1 & schicht==3

replace W_pt = W_pt / sc_strata if D_pt==1

* Extrapolating from the sample to the top 1% of the wealth distribution in Germany 
* (on HH-level, 2017: about 41,304,000 HHs * 0.01 = 413,000)
* (compare to "https://de.statista.com/statistik/daten/studie/156950/umfrage/anzahl-der-privathaushalte-in-deutschland-seit-1991/")
qui sum _1_nw if D_pt==1
replace W_pt = W_pt * (413000/r(N))

* check frequency weights
sum _1_nw if D_pt==0 [fw = round(W_sp)]
sum _1_nw if D_pt==1 [fw = round(W_pt)]


save "${outpath}soep_pretest_1.dta", replace


***
