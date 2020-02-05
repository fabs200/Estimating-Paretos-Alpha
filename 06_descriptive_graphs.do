*** 6. graphs ***

/* 
Here we plot histograms and CDFs. For simplicity, we use only the 1st 
implicate (_1_nw_mio) of the SOEP and of the Pretest.

Here we produce: Figure 2.1, Figure A.7
*/

set graph off

use "${outpath}soep_pretest_2.dta", clear

********************************************************************************
*		
*	6.1 histograms: SOEP vs. Pretest
*		(in paper: Figure 2.1)
*		
********************************************************************************

/*
Note: We loop twice, first, we plot the histogram with nw_mio > 0 (s=0) and 
	  second, with nw_mio > 5 (s=5).
*/

forval s=0(5)5 {
	twoway	(histogram _1_nw_mio if _1_nw_mio>`s' & D_pt==0, width(15) color(gray%50) freq) ///
		(histogram _1_nw_mio if _1_nw_mio>`s' & D_pt==1, width(15) color(cranberry%35)		///
        xtitle("Net Wealth (in Mio. Euro)")													///
		legend(label(1 "SOEP 2012") label(2 "Pretest 2017")) 				 				///
		note("Note: Net wealth (nw) in million euros; nw larger than `s'm, implicate 1 plotted.") ///
		ytitle("frequency") ylabel(, format(%9.0fc)) freq									///
		scheme(s2mono) graphregion(color(white)))
	graph export "${graphs}6_1_histogram_nw_mio_larger`s'_sp_pt.pdf", replace
}

********************************************************************************
*		
*	6.2 CDF: SOEP (_sp) and Pretest (_pt)
*		(in paper: Figure A.7)
*		
********************************************************************************

/*
Note: We plot the CDF. We provide the CDFs in the appendix of our paper.
*/

foreach data in sp pt {
	
	if "`data'" == "sp" {
		local cumul_title "CDF of Net Wealth of SOEP"
		local cumul_col "gray"
		local cumul_dat "D_pt==0"
		local cumul_leg "SOEP 2012"
	}
	if "`data'" == "pt" {
		local cumul_title "CDF of Net Wealth of Pretest"
		local cumul_col "cranberry%75"
		local cumul_dat "D_pt==1"
		local cumul_leg "Pretest 2017"
	}

	* prepare CDF variable (nw implicate 1)
	cumul _1_nw_mio [fw=round(W_`data')] if _1_nw_mio>=0, gen(_1_nw_mio_cumul_`data')
	sort _1_nw_mio_cumul_`data'
	
	* plot CDF
	line _1_nw_mio_cumul_`data' _1_nw_mio if `cumul_dat', ylab(, grid) ytitle("")		///
		lpattern(--) lcolor(`cumul_col')xlab(, grid) xtitle("Net Wealth (in Mio. Euro)") ///
		note("Source: `cumul_leg'." "Note: CDF of nw implicate 1 plotted, CDF calculated with nw larger 0.") ///
		/* title("`cumul_title'") saving(${graphs}cumul_`data', replace) */				///
		scheme(s2mono) graphregion(color(white))
	
	graph export "${graphs}6_2_cumul_`data'.pdf", replace
	
	drop _1_nw_mio_cumul_`data'
}

set graph on

***
