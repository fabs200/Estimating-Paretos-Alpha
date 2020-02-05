*** 5. Descriptive Statistics ***

/* 
We create summary statistics for the weighted SOEP and Pretest datasets 
including the respective weights as applied in do-file 02_prep_weights.do.

Here, we produce: table 2.2., table A.1, table A.2, table A.3
*/

********************************************************************************
* OUTLINE
* 5.1. Summary statistics of ALL net wealth implicates of SOEP and Pretest
* 5.2. Summary statistics of ALL net wealth imp., weighted, SOEP and Pretest
* 5.3. Summary statistics only for 1st net wealth imp. of SOEP and Pretest
* 5.4. Summary statistics only for 1st net wealth imp. weighted, SOEP and Pretest
********************************************************************************

use "${outpath}soep_pretest_2.dta", clear

* format variables for better display
doubletofloat _*_nw


********************************************************************************
*
*	5.1. Summary Statistics ALL net wealth SOEP (_sp) and Pretest (_pt)
*		(in paper: table A.1)
*
********************************************************************************

* define display formats
local fmt_count "%10.0fc"
local fmt_stats "%10.2fc"

/* 
Note: the code loops twice: First, i=1, meaning that nw including neg. values.
Second, i=2, meaning that we only include pos. values (see local option_`data').
*/

forval i=1(1)2 {
	
	*local i=1
	* generate used locals for 'estpost summarize'
	local estpost_sp	""
	local estpost_pt	""

	* nw full range (incl. neg. values)
	local option1_sp	"if D_pt == 0"
	local option1_pt	"if D_pt == 1"

	* nw>0 (only pos. values)
	local option2_sp	"if _1_nw_mio > 0 & D_pt == 0"
	local option2_pt	"if _1_nw_mio > 0 & D_pt == 1"

	* sources
	local source_sp		"SOEP 2012 (v33.1)"
	local source_pt		"Pretest 2017"
	
	* generate estpost locals and variables
	foreach data in sp pt {
		
		forval imp=1(1)5 {
			local estpost_`data' "`estpost_`data'' _`imp'_nw_mio_`data'"
			* label
			local label_sp		"SOEP nw imp. `imp'"
			local label_pt		"Pretest nw imp. `imp'"

			di "`estpost_`data''"
			gen _`imp'_nw_mio_`data' = _`imp'_nw_mio `option`i'_`data''
			label variable _`imp'_nw_mio_`data' "`label_`data''"
		}
	
		* sample sizes
		qui sum _1_nw_mio_`data' `option`i'_`data''
		scalar sc_N_i`i'_`data'=r(N)

	}
	
	* generate summary statistics
	estpost summarize `estpost_sp' `estpost_pt', detail

	esttab using "${tables}5_1_summary_statistics`i'_in_mio_noW_ALL.tex", 		///
				replace varwidth(44) nonumber noobs nomtitles label				///
				cells("count(fmt(`fmt_count')) mean(fmt(`fmt_stats')) p50(fmt(`fmt_stats')) p75(fmt(`fmt_stats')) p90(fmt(`fmt_stats')) p99(fmt(`fmt_stats')) min(fmt(`fmt_stats')) max(fmt(`fmt_stats'))" ) ///
				title(Descriptive Statistics of `source_sp' and `source_pt')	///
				addnote("Source: `source_sp', `source_pt'." "Note: Net wealth (nw) imputed, in mio. Euro, unweighted, for simplicity rounded.")


	****************************************************************************
	*
	* 5.2. Summary statistics of ALL net wealth imp., weighted, SOEP and Pretest
	* (in paper: table A.2)
	* (in paper: table A.3)
	*
	****************************************************************************

	
	local W_sp_info		"SOEP: with applied frequency weights on household-level (N=`=sc_N_i`i'_sp')"
	local W_pt_info		"Pretest: with own re-weighting scheme (N=`=sc_N_i`i'_pt')"
		
	foreach data in sp pt {
	
		estpost summarize `estpost_`data'' [fw=round(W_`data')] `option`i'_`data'', detail 

		esttab using "${tables}5_2_summary_statistics`i'_in_mio_`data'_W_ALL.tex", 	///
				replace varwidth(44) nonumber noobs nomtitles label					///
				cells("count(fmt(`fmt_count')) mean(fmt(`fmt_stats')) p50(fmt(`fmt_stats')) p75(fmt(`fmt_stats')) p90(fmt(`fmt_stats')) p99(fmt(`fmt_stats')) min(fmt(`fmt_stats')) max(fmt(`fmt_stats'))" ) ///
				title(Descriptive Statistics of `source_`data'')					///
				addnote("Source: `source_`data''." "Note: Net wealth (nw), imputed, weighted and displayed in mio. Euro, for simplicity rounded." "`W_`data'_info'.")
				
	}

	drop `estpost_sp' `estpost_pt'

}


********************************************************************************
*
*	5.3. Summary statistics only for 1st net wealth implicate of SOEP (_sp) and 
*		 Pretest (_pt) (in paper: table 2.2.)
*
********************************************************************************

/*
Note: Due to technical reasons, we redo section 5.1 but only for the 1st 
	  implicate of SOEP, Pretest. (indicated as "_imp1only")
*/

* define display formats
local fmt_count "%10.0fc"
local fmt_stats "%10.2fc"

forval i=1(1)2 {
	
	*local i=1
	* generate used locals for 'estpost summarize'
	local estpost_sp	""
	local estpost_pt	""

	* nw full range
	local option1_sp	"if D_pt == 0"
	local option1_pt	"if D_pt == 1"

	* nw>0
	local option2_sp	"if _1_nw_mio > 0 & D_pt == 0"
	local option2_pt	"if _1_nw_mio > 0 & D_pt == 1"

	* sources
	local source_sp		"SOEP 2012 (v33.1)"
	local source_pt		"Pretest 2017"
	
	* generate estpost locals and variables
	foreach data in sp pt {
		
		forval imp=1(1)1 {
			local estpost_`data' "`estpost_`data'' _`imp'_nw_mio_`data'"
			* label
			local label_sp		"SOEP nw imp. `imp'"
			local label_pt		"Pretest nw imp. `imp'"

			di "`estpost_`data''"
			gen _`imp'_nw_mio_`data' = _`imp'_nw_mio `option`i'_`data''
			label variable _`imp'_nw_mio_`data' "`label_`data''"
		}
	
		* sample sizes
		qui sum _1_nw_mio_`data' `option`i'_`data''
		scalar sc_N_i`i'_`data'=r(N)

	}
	
	* generate summary statistics
	estpost summarize `estpost_sp' `estpost_pt', detail

	esttab using "${tables}5_3_summary_statistics`i'_in_mio_noW_imp1only.tex", 	///
				replace varwidth(44) nonumber noobs nomtitles label				///
				cells("count(fmt(`fmt_count')) mean(fmt(`fmt_stats')) p50(fmt(`fmt_stats')) p75(fmt(`fmt_stats')) p90(fmt(`fmt_stats')) p99(fmt(`fmt_stats')) min(fmt(`fmt_stats')) max(fmt(`fmt_stats'))" ) ///
				title(Descriptive Statistics of `source_sp' and `source_pt')	///
				addnote("Source: `source_sp', `source_pt'." "Note: Net Wealth (nw) implicate 1, in mio. Euro, unweighted, for simplicity rounded.")

	****************************************************************************
	*
	* 	5.4. Summary statistics of 1st net wealth implicate, weighted, SOEP, Pretest
	*
	****************************************************************************
	/*
	Note: Due to technical reasons, we redo section 5.2 but only for the 1st 
		  implicate of SOEP, Pretest. (indicated as "_imp1only")
	*/


	
	local W_sp_info		"SOEP: with applied frequency weights on household-level (N=`=sc_N_i`i'_sp')"
	local W_pt_info		"Pretest: with own re-weighting scheme (N=`=sc_N_i`i'_pt')"
		
	foreach data in sp pt {
	
		estpost summarize `estpost_`data'' [fw=round(W_`data')] `option`i'_`data'', detail 

		esttab using "${tables}5_4_summary_statistics`i'_in_mio_`data'_W_imp1only.tex", 	///
				replace varwidth(44) nonumber noobs nomtitles label						///
				cells("count(fmt(`fmt_count')) mean(fmt(`fmt_stats')) p50(fmt(`fmt_stats')) p75(fmt(`fmt_stats')) p90(fmt(`fmt_stats')) p99(fmt(`fmt_stats')) min(fmt(`fmt_stats')) max(fmt(`fmt_stats'))" ) ///
				/* title(Descriptive Statistics of `source_`data'') */						///
				addnote("Source: `source_`data''." "Note: Net wealth (nw), imputed, weighted and displayed in mio. Euro, for simplicity rounded." "`W_`data'_info'.")
				
	}

	drop `estpost_sp' `estpost_pt'

}


***
