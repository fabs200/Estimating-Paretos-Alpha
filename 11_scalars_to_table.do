*** 11. Scalars to table ***

/*
Here we save all calculated scalars to tables which are automatically save as
.tex-files.

Here we produce: table 4.1, table 4.2, table 4.3, table 4.4

*/

********************************************************************************
* OUTLINE
* Table 11.1: Alphas Result (table 4.1 and table 4.2)
* Table 11.2: Calculation Percentiles (table 4.3)
* Table 11.3: Calculation Total Sums Wealth (table 4.4)
********************************************************************************

clear
set matsize 11000
scalar m=5
local pctile p95 p99
local data sp pt

********************************************************************************
* Table 11.1: Alphas Result
* (in paper: table 4.1 and table 4.2)
********************************************************************************

foreach pct of local pctile {

	/*
	table TEMPLATE (each p95, p99)
					/ alpha_sp	/ N_sp	/ alpha_pt	/ N_pt	/ thres_val	/ pval-ci	/
	1	nw1			/			/		/			/		/			/			/
	2	sd1			/			/		/			/		/			/			/
	3	nw2			/			/		/			/		/			/			/			
	4	sd2			/			/		/			/		/			/			/
	5	nw3			/			/		/			/		/			/			/			
	6	sd3			/			/		/			/		/			/			/
	7	nw4			/			/		/			/		/			/			/			
	8	sd4			/			/		/			/		/			/			/
	9	nw5			/			/		/			/		/			/			/			
	10	sd5			/			/		/			/		/			/			/
	11	avgalpha	/			/		/			/		/			/sc_pval_pct/
	12	sd_avgalpha	/			/		/			/		/			/			/
	13	ci lb		/			/		/			/		/			/			/
	14	ci up		/			/		/			/		/			/			/
	
	NOTE: Additionally, testing CI_p95_sp = CI_p99_sp delivers the p-value: xxx 
		and CI_p95_pt = CI_p99_pt delivers the p-value: xxx.
	*/
	
	* define matrix: scalar results
	matrix R = J(14, 7 ,.)
	matrix colnames R = var alpha_sp N_sp alpha_pt N_pt thres_val_`pct' p-val
	matrix rownames R = nw1 sd1 nw2 sd2 nw3 sd3 nw4 sd4 nw5 sd5 avgalpha tv cilo cihi
	*matrix list R

	* pct without 'p'
	scalar PC = substr("`pct'",2,3)
	
	scalar b = 2
	scalar c = 3
	foreach dat of local data {
		scalar a = 1
		forval imp=1(1)`=m' {
		
			* display scalars
			if "$show_all" == "TRUE" {
			di in red "scalar sc_alpha_`imp'_`pct'_`dat' = `=sc_alpha_`imp'_`pct'_`dat''"
			di in red "matrix: R[`=a',`=b']"
			di in red "scalar sc_N_`imp'_`pct'_`dat' = `=sc_N_`imp'_`pct'_`dat''"
			di in red "matrix: R[`=a',`=c']"
			}

			* write scalars into matrix R
			matrix R[`=a',`=b'] = round(`=sc_alpha_`imp'_`pct'_`dat'',.001)
			matrix R[11,`=b'] 	= round(`=sc_avgalpha_`pct'_`dat'',.001)
			matrix R[`=a',`=c'] = `=sc_N_`imp'_`pct'_`dat''
			matrix R[`=a',6] 	= `=sc_thres_`imp'_`pct''

			scalar a = `=a' + 1

			matrix R[`=a',`=b'] = round(`=sc_sd_`imp'_`pct'_`dat'',.0001)
			matrix R[12,`=b'] 	= round(`=sc_tv_`pct'_`dat'',.0001)

			matrix R[11,7] 		= round(`=sc_pval_`pct'',.00001)
			matrix R[13,`=b'] 	= round(`=sc_cilo_`pct'_`dat'',.001)
			matrix R[14,`=b'] 	= round(`=sc_cihi_`pct'_`dat'',.001)

			scalar a = `=a' + 1
		}
		
		scalar b = `=b' + 2
		scalar c = `=c' + 2

	}

	* transform matrix to Stata-file and process
	matrix list R
	drop _all
	svmat double R
	matrix drop R

	* rename
	ren (R2 R3 R4 R5 R6 R7) (alpha_sp N_sp alpha_pt N_pt thres_val_`pct' pval)

	save "${outpath}alpha_results_`pct'_temp.dta", replace
	clear

	use "${outpath}alpha_results_`pct'_temp.dta", clear

	* label
	cap qui label drop alpha_results
	scalar sc_alphaCI_left = 1-`=sc_alphaCI'
	label define alpha_results 1 "ln(nw_{1})" 2 "" 3 "ln(nw_{2})" 4 "" 5 "ln(nw_{3})" 6 "" 7 "ln(nw_{4})" 8 "" 9 "ln(nw_{5})" 10 "" 11 "mean of $\alpha$" 12 "" 13 "CI lower (`=sc_alphaCI_left')" 14 "CI upper (`=sc_alphaCI')"
	replace R1 = _n
	label values R1 alpha_results
	label variable R1 ""

	* prepare sd in parenthesis
	foreach dat of local data {
		tostring alpha_`dat', replace force
		replace alpha_`dat' = "(" + alpha_`dat' + ")" if !mod(_n,2) & _n!=14
	}

	* check
	if "$show_all" == "TRUE" {
	l
	}
	
	*** Latex table ***
	listtex * using "${tables}11_1_alpha_results_`pct'.tex", replace type rstyle(tabular) missnum() ///
			head("\begin{table} \centering \begin{tabular}{ccccccc}" `"\hline \textit{} &	\textit{$\alpha_{SOEP}$} &	\textit{$\N_{SOEP}$} &	\textit{$\alpha_{Pretest}$} &	\textit{$\N_{Pretest}$} & \textit{threshold} &	\textit{p-value} \\ \hline"') ///
			foot("\hline" "\multicolumn{7}{l}{%" "\begin{minipage}{10cm}%"  "\vspace{.1cm} Source: SOEP (2012, v33.1) and Pretest 2017. \\Note: p-value represents the test result of xxx nw: net wealth, sd in parenthesis, household weights applied to SOEP, own calculated weights applied to the Pretest \end{minipage}%" "} \\  \end{tabular} \caption{Estimation Results of Pareto's Coefficients (Threshold at `pct')} \label{tab:alpha_results_`pct'} \end{table}")

}

********************************************************************************
* Table 11.2: Calculation Percentiles
* (in paper: table 4.3)
********************************************************************************

/*
table TEMPLATE (each p95, p99)
	
		/			/  		pct2 with SOEP alpha		/ 	pct with pretest alpha			/
	1	/			/ 	cilo	/ 	mean	/ 	cihi	/ 	cilo	/ 	mean	/ 	cihi	/
	2	/ Top 5%	/			/			/			/	(x)		/	 (x)	/	 (x)	/
	3	/ Top 1%	/			/			/			/			/			/			/
	4	/ Top 0.1%	/			/			/			/			/			/			/

			(Note: no values for Top 5% if p99: (x))
*/

foreach pct of local pctile {
	
	* calculate means between cilo and cihi
	if "`pct'" == "p95" {
	scalar sc_pcttop50_`pct'_mean = (sc_pcttop50_`pct'_cilo + sc_pcttop50_`pct'_cihi)/2
	scalar sc_pct2top50_`pct'_mean = (sc_pct2top50_`pct'_cilo + sc_pct2top50_`pct'_cihi)/2
	}
	scalar sc_pcttop10_`pct'_mean = (sc_pcttop10_`pct'_cilo + sc_pcttop10_`pct'_cihi)/2
	scalar sc_pct2top10_`pct'_mean = (sc_pct2top10_`pct'_cilo + sc_pct2top10_`pct'_cihi)/2
	
	scalar sc_pcttop01_`pct'_mean = (sc_pcttop01_`pct'_cilo + sc_pcttop01_`pct'_cihi)/2
	scalar sc_pct2top01_`pct'_mean = (sc_pct2top01_`pct'_cilo + sc_pct2top01_`pct'_cihi)/2

	* define matrix: scalar results
	capture drop matrix P
	matrix P = J(3, 7,.)
	matrix colnames P = top cilo_sp_alph mean_sp_alph cihi_sp_alph cilo_pt_alph mean_pt_alph cihi_pt_alph
	matrix rownames P = top50 top10 top01
	matrix list P

	if "`pct'" == "p95" {
	matrix P[1,2] = round(`=sc_pct2top50_`pct'_cihi',.001)
	matrix P[1,3] = round(`=sc_pct2top50_`pct'_mean',.001)
	matrix P[1,4] = round(`=sc_pct2top50_`pct'_cilo',.001)
	matrix P[1,5] = round(`=sc_pcttop50_`pct'_cihi',.001)
	matrix P[1,6] = round(`=sc_pcttop50_`pct'_mean',.001)
	matrix P[1,7] = round(`=sc_pcttop50_`pct'_cilo',.001)
	}
	matrix P[2,2] = round(`=sc_pct2top10_`pct'_cihi',.001)
	matrix P[2,3] = round(`=sc_pct2top10_`pct'_mean',.001)
	matrix P[2,4] = round(`=sc_pct2top10_`pct'_cilo',.001)
	matrix P[2,5] = round(`=sc_pcttop10_`pct'_cihi',.001)
	matrix P[2,6] = round(`=sc_pcttop10_`pct'_mean',.001)
	matrix P[2,7] = round(`=sc_pcttop10_`pct'_cilo',.001)
	matrix P[3,2] = round(`=sc_pct2top01_`pct'_cihi',.001)
	matrix P[3,3] = round(`=sc_pct2top01_`pct'_mean',.001)
	matrix P[3,4] = round(`=sc_pct2top01_`pct'_cilo',.001)
	matrix P[3,5] = round(`=sc_pcttop01_`pct'_cihi',.001)
	matrix P[3,6] = round(`=sc_pcttop01_`pct'_mean',.001)
	matrix P[3,7] = round(`=sc_pcttop01_`pct'_cilo',.001)

	matrix list P

	* transform matrix to Stata-file and process
	matrix list P
	drop _all
	svmat double P
	matrix drop P

	* rename
	ren (P2 P3 P4 P5 P6 P7) (cilo_sp_alph mean_sp_alph cihi_sp_alph cilo_pt_alph mean_pt_alph cihi_pt_alph)

	* label
	cap qui label drop label_predict_pcts
	label define label_predict_pcts 1 "Top 5\%" 2 "Top 1\%" 3 "Top 0.1\%"
	replace P1 = _n
	label values P1 label_predict_pcts
	label variable P1 ""

	* check
	if "$show_all" == "TRUE" {
	l
	}

	*** Latex table ***
	listtex * using "${tables}11_2_predict_top_pct_`pct'.tex", replace type rstyle(tabular) missnum() ///
			head("\begin{table} \centering \begin{tabular}{c|ccc|ccc}" `"\hline \textit{} &	\textit{$\CI_{low}$} & \textit{mean} &	\textit{$\CI_{high}$} &	\textit{$\CI_{low}$} & \textit{mean} &	\textit{$\CI_{high}$} \\ \hline"') ///
			foot("\hline \multicolumn{7}{l}{%" "\begin{minipage}{10cm}%"  "\vspace{.1cm} Source: SOEP (2012, v33.1) and Pretest 2017. \\Note: All values in mio. euros, 95\% confidence intervals, calculated according to Rubin (1987), estimation based on our imputation, threshold at `pct'. \end{minipage}%" "} \\  \end{tabular} \caption{Estimation of Selected Percentiles (Threshold at `pct')} \label{tab:predict_pcts_`pct'} \end{table}")

}

********************************************************************************
* Table 11.3: Calculation Total Sums Wealth
* (in paper: table 4.4)
********************************************************************************

/*
table TEMPLATE: (each p95, p99)
	
		/			/ 		tw2 with SOEP alpha			/ 		tw with pretest alpha		/
	1	/			/ 	cilo	/ 	mean	/ 	cihi	/ 	cilo	/ 	mean	/ 	cihi	/
	2	/ Top 5%	/			/			/			/	(x)		/	(x)		/	 (x)	/
	3	/ Top 1%	/			/			/			/			/			/			/
	4	/ Top 0.1%	/			/			/			/			/			/			/

			(Note: no values for Top 5% if p99: (x))
*/

foreach pct of local pctile {
	
	* calculate means between cilo and cihi
	if "`pct'" == "p95" {
	scalar sc_twtop50_`pct'_mean = (sc_twtop50_`pct'_cilo + sc_twtop50_`pct'_cihi)/2
	scalar sc_tw2top50_`pct'_mean = (sc_tw2top50_`pct'_cilo + sc_tw2top50_`pct'_cihi)/2	
	}
	scalar sc_twtop10_`pct'_mean = (sc_twtop10_`pct'_cilo + sc_twtop10_`pct'_cihi)/2
	scalar sc_tw2top10_`pct'_mean = (sc_tw2top10_`pct'_cilo + sc_tw2top10_`pct'_cihi)/2
	
	scalar sc_twtop01_`pct'_mean = (sc_twtop01_`pct'_cilo + sc_twtop01_`pct'_cihi)/2
	scalar sc_tw2top01_`pct'_mean = (sc_tw2top01_`pct'_cilo + sc_tw2top01_`pct'_cihi)/2

	* define matrix: scalar results
	capture drop matrix T
	matrix T = J(3, 7,.)
	matrix colnames T = top cilo_sp_alph mean_sp_alph cihi_sp_alph cilo_pt_alph mean_pt_alph cihi_pt_alph
	matrix rownames T = top50 top10 top01
	matrix list T

	if "`pct'" == "p95" {
	matrix T[1,2] = round(`=sc_tw2top50_`pct'_cihi',.001)
	matrix T[1,3] = round(`=sc_tw2top50_`pct'_mean',.001)
	matrix T[1,4] = round(`=sc_tw2top50_`pct'_cilo',.001)
	matrix T[1,5] = round(`=sc_twtop50_`pct'_cihi',.001)
	matrix T[1,6] = round(`=sc_twtop50_`pct'_mean',.001)
	matrix T[1,7] = round(`=sc_twtop50_`pct'_cilo',.001)
	}
	matrix T[2,2] = round(`=sc_tw2top10_`pct'_cihi',.001)
	matrix T[2,3] = round(`=sc_tw2top10_`pct'_mean',.001)
	matrix T[2,4] = round(`=sc_tw2top10_`pct'_cilo',.001)
	matrix T[2,5] = round(`=sc_twtop10_`pct'_cihi',.001)
	matrix T[2,6] = round(`=sc_twtop10_`pct'_mean',.001)
	matrix T[2,7] = round(`=sc_twtop10_`pct'_cilo',.001)
	matrix T[3,2] = round(`=sc_tw2top01_`pct'_cihi',.001)
	matrix T[3,3] = round(`=sc_tw2top01_`pct'_mean',.001)
	matrix T[3,4] = round(`=sc_tw2top01_`pct'_cilo',.001)
	matrix T[3,5] = round(`=sc_twtop01_`pct'_cihi',.001)
	matrix T[3,6] = round(`=sc_twtop01_`pct'_mean',.001)
	matrix T[3,7] = round(`=sc_twtop01_`pct'_cilo',.001)

	matrix list T

	* transform matrix to Stata-file and Process
	matrix list T
	drop _all
	svmat double T
	matrix drop T

	* rename
	ren (T2 T3 T4 T5 T6 T7) (cilo_sp_alph mean_sp_alph cihi_sp_alph cilo_pt_alph mean_pt_alph cihi_pt_alph)

	* label
	cap qui label drop label_predict_tws
	label define label_predict_tws 1 "Top 5\%" 2 "Top 1\%" 3 "Top 0.1\%"
	replace T1 = _n
	label values T1 label_predict_tws
	label variable T1 ""

	* check
	if "$show_all" == "TRUE" {
	l
	}

	*** Latex table ***
	listtex * using "${tables}11_3_predict_top_tw_`pct'.tex", replace type rstyle(tabular) missnum() ///
			head("\begin{table} \centering \begin{tabular}{c|ccc|ccc}" `"\hline \textit{} &	\textit{$\CI_{low}$} & \textit{mean} &	\textit{$\CI_{high}$} &	\textit{$\CI_{low}$} & \textit{mean} &	\textit{$\CI_{high}$} \\ \hline"') ///
			foot("\hline \multicolumn{7}{l}{%" "\begin{minipage}{10cm}%"  "\vspace{.1cm} Source: SOEP (2012, v33.1) and Pretest 2017. \\Note: All values in billions euros, 95\% confidence intervals, calculated according to Rubin (1987), estimation based on our imputation, threshold at `pct'. \end{minipage}%" "} \\  \end{tabular} \caption{Calculation of Total Sum of Net Wealth (Threshold at `pct')} \label{tab:predict_tws_`pct'} \end{table}")

}



		
		
/*
*** NOTE for Latex table: ***
** add packages:
	\usepackage[euler]{textgreek}
	\usepackage{graphicx}
	\usepackage{makecell}

** to place table accurately on page: 
	\resizebox{.5\width}{!}{\input{Seminararbeit/img/alpha_results_p99_edited.tex}}
*/


***
