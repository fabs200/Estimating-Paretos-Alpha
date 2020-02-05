********************************************************************************
* Master
********************************************************************************

clear all
set more off
set graph off

* define project-path
if "`c(username)'" == "Fabian" {
	global mypath	"/Users/Fabian/OneDrive/Studium/Seminar/"
	global soep 	"/Users/Fabian/Documents/DATA/STATA/SOEP_v33.1/SOEP_wide/"
	global pretest	"/Users/Fabian/Documents/DATA/STATA/pretest/"
}
else if "`c(username)'" == "avila" {
	global mypath	"/home/avila/Documents/Projects/WS1819_FU/MeasuringInequality/stata/"
	global soep		"/data/DatasetsSOEP/SOEP_v33.1/SOEP_wide/"
	global pretest	"/data/DatasetsSOEP/pretest/"
}
else if "`c(username)'" == "Tobias?" {
	global mypath	"path"
	global soep 	"path"
	global pretest	"path"
}
else if "`c(username)'" == "Sebastian?" {
	global mypath	"path"
	global soep 	"path"
	global pretest	"path"
}

* globals of further paths
global do		"${mypath}do/"
global outpath	"${mypath}outpath/"
global graphs	"${outpath}graphs/"
global tables	"${outpath}tables/"

if "`c(username)'" == "Fabian" {
	global outpath	"/Users/Fabian/Documents/DATA/STATA/"
}

* create folders if not existing
foreach dir in outpath graphs tables {
	if "`dir'" == "outpath" {
		capture confirm file "${mypath}`dir'"
		if _rc mkdir "${mypath}`dir'/"
	}
	else { 
		capture confirm file "${outpath}`dir'"
		if _rc mkdir "${outpath}`dir'/"
	}
}

* display all generated scalars
global show_all "TRUE"

***


*** I: Data Preparation
	* 1. SOEP, Pretest data preparation
	do "${do}01_prep_data.do"

	* 2. preparation of weights (Pretest)
	do "${do}02_prep_weights.do"

	* 3. preparation of variables used for pareto distribution
	do "${do}03_prep_pareto.do"

	* 4. preparation of mi impute dataset
	do "${do}04_prep_mi_impute_dataset.do"

	
*** II: Descriptive Statistics and Graphs
	* 5. descriptive statistics
	do "${do}05_descriptive_stats.do"

	* 6. additional graphs
	do "${do}06_descriptive_graphs.do"

	
*** III: Analysis
	* 7. fitting pareto distribution and predict top percentiles of SOEP
	do "${do}07_fitting_pareto.do"

	* 8. Robustness Checks
	do "${do}08_robustness_checks.do"

	* 9. Impute Top-Tails of SOEP
	do "${do}09_impute_top_percentiles.do"
	
	* 10. threshold-alhpa-relation
	do "${do}10_threshold_alpha.do"

	
*** IV: Data Export
	* 11. Exporting scalars as tables
	do "${do}11_scalars_to_table.do"

	* 12. Pareto Graphs
	do "${do}12_pareto_graphs.do"


set graph on

***
