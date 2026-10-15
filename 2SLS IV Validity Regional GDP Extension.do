/****************************************************************************************
 0. Setup
****************************************************************************************/
clear all                  // Remove all data and matrices from memory
set matsize 5000           // Increase maximum matrix size (needed for xtivreg2 with many FEs)
set more off               // Prevent Stata from pausing output

***Install all necessary packages***
******************************************
cap ssc install xtabond2 
cap ssc install xtivreg2 
cap ssc install spmat 
cap ado install st0292, from("http://www.stata-journal.com/software/sj13-2/") replace
cap ssc install spmack
cap ssc install ivreg2, replace
cap ssc install ranktest, replace

***Set Base Directory***
******************************************
global project "C:/Users/sohom/OneDrive/Desktop/Courses/Foundations of Econometrics & ML/Group Project/replication_files_ddcg"  /* Set base directory                         */
global limit=25                  // Evaluate effects 25 years after transition    
local repsBS=100                 // Number of bootstrap repetitions

/****************************************************************************************
 1. DEFINE PROGRAMS & FUNCTIONS THAT WILL BE CALLED UPON IN THE PROGRAM
****************************************************************************************/


/* The purpose of these programs are to compute metrics of interest using 
the regression results*/


/*Vareffects function for posting outputs of modified columns 
in second-stage regression*/
**************************************************************************************

capture program drop vareffects_with_yreg
program define vareffects_with_yreg, eclass

quietly: nlcom (effect1: _b[shortrun]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post

quietly: nlcom (effect2: _b[effect1]*_b[lag1]+_b[shortrun]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post

quietly: nlcom (effect3: _b[effect2]*_b[lag1]+_b[effect1]*_b[lag2]+_b[shortrun]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post
	  
quietly: nlcom (effect4: _b[effect3]*_b[lag1]+_b[effect2]*_b[lag2]+_b[effect1]*_b[lag3]+_b[shortrun]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post	  

forvalues j=5(1)$limit{	  
local j1=`j'-1
local j2=`j'-2
local j3=`j'-3
local j4=`j'-4

quietly: nlcom (effect`j': _b[effect`j1']*_b[lag1]+_b[effect`j2']*_b[lag2]+_b[effect`j3']*_b[lag3]+_b[effect`j4']*_b[lag4]+_b[shortrun]) ///
	  (effect`j1': _b[effect`j1']) ///
	  (effect`j2': _b[effect`j2']) ///
	  (effect`j3': _b[effect`j3']) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post	  	  

}

quietly: nlcom (effect$limit: _b[effect$limit]) ///
	  (longrun: _b[shortrun]/(1-_b[lag1]-_b[lag2]-_b[lag3]-_b[lag4])) ///
      (shortrun: _b[shortrun]) ///
	  (persistence: _b[lag1]+_b[lag2]+_b[lag3]+_b[lag4]) ///
	  (cumulativereg: _b[yreg_lag1]+_b[yreg_lag2]+_b[yreg_lag3]+_b[yreg_lag4]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post
ereturn display
end


/*Vareffects function for posting outputs of modified columns 
in first-stage regression*/
****************************************************************************************

capture program drop vareffects_with_yreg_d
program define vareffects_with_yreg_d, eclass

quietly: nlcom (persistence: _b[lag1]+_b[lag2]+_b[lag3]+_b[lag4]) ///
	  (cumulativereg: _b[yreg_lag1]+_b[yreg_lag2]+_b[yreg_lag3]+_b[yreg_lag4]) ///
	  (cumulativeinstr: _b[instr_lag1]+_b[instr_lag2]+_b[instr_lag3]+_b[instr_lag4]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (yreg_lag1: _b[yreg_lag1]) ///
	  (yreg_lag2: _b[yreg_lag2]) ///
	  (yreg_lag3: _b[yreg_lag3]) ///
	  (yreg_lag4: _b[yreg_lag4]) ///
	  , post
ereturn display
end

/*Vareffects function for posting outputs of instrument coefficients in 
authhors' original first-stage regression*/
*************************************************************************************

capture program drop vareffects_original_d
program define vareffects_original_d, eclass

quietly: nlcom (persistence: _b[lag1]+_b[lag2]+_b[lag3]+_b[lag4]) ///
	  (cumulativeinstr: _b[instr_lag1]+_b[instr_lag2]+_b[instr_lag3]+_b[instr_lag4]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
ereturn display
end



/****************************************************************************************
 2. TOP PANEL: 2SLS ESTIMATES — GLOBAL RESULTS (WITH YREG EXTENSION)
****************************************************************************************/

use "$project/DDCGdata_final.dta", clear
capture drop instrument
gen instrument=demreg

*Create Spatial Lags for Columns 8 and 9*
*********************************************************************

/* Purpose:
   Construct time-varying spatial lag matrices for:
       y_w      = spatial lag of GDP
       dem_w    = spatial lag of democracy
       yreg_w   = spatial lag of regional GDP  
*/

forvalues j=1960(1)2010{
preserve
keep if year==`j'
keep if cen_lon!=.&cen_lat!=.&y!=.&dem!=.
keep wbcode2 cen_lon cen_lat y yreg dem year 
sort wbcode2

spmat idistance dmat cen_lon cen_lat, id(wbcode2) dfunction(dhaversine) norm(row) replace
spmat lag  y_w dmat y
spmat lag  dem_w dmat dem
spmat lag  yreg_w dmat yreg

keep y_w dem_w yreg_w wbcode2 year
tempfile m
save `m', replace
restore

merge 1:1 wbcode2 year using `m', update
drop _m

}

sort wbcode2 year


/**********************************************************
Top panel: 2SLS estimates — GLOBAL Results (With Extension)
of Columns of Interest (2 & 5-9)
***********************************************************/


*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*, fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7, fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy*, fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy*, fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy*, fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy*, fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add


***********************************************/
*Export Regression Tables as csv files*
***********************************************/

#delimit ;
esttab e2 e5 e6 e7 e8 e9  using "$project/results/Table2SLSTop.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))   
keep(dem) order(dem)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add e5_add e6_add e7_add e8_add e9_add using 
"$project/results/Table2SLSTop_Add.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg) 
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr

/**************************************************************************
First Stage Results of Columns of Interest 2 & 5-9
(Reproduce Original Authors' Results but display Instrument Coefficients)
***************************************************************************/

*Column 2*
**********************************
quietly: xtreg dem l(1/4).y l(1/4).instrument yy* if samp2==1, fe cluster(wbcode2) r 
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_original_d
estimates store e2fs_o_add

*Column 5*
**********************************
quietly: xtreg dem l(1/4).y l(1/4).instrument yy*  rtrend* if samp5==1, fe cluster(wbcode2) r 
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_original_d
estimates store e5fs_o_add

*Column 6*
**********************************
quietly: xtreg dem l(1/4).y l(1/4).instrument l(1/4).(yreg tradewbreg) yy* if samp6==1, fe cluster(wbcode2) r 
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_original_d
estimates store e6fs_o_add

*Column 7*
**********************************
quietly: xtreg dem l(1/4).y l(1/4).instrument l(1/4).(yreg tradewbreg unrestreg) yy* if samp7==1, fe cluster(wbcode2) r
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_original_d
estimates store e7fs_o_add

*Column 8*
**********************************
quietly: xtreg dem l(1/4).y l(1/4).instrument l(1/4).y_w yy*  if samp8==1, fe cluster(wbcode2) r 
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_original_d
estimates store e8fs_o_add

*Column 9*
**********************************
quietly: xtreg dem l(1/4).y l(1/4).instrument l(1/5).y_w l(1/4).dem_w yy*  if samp9==1, fe cluster(wbcode2) r 
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_original_d
estimates store e9fs_o_add

**************************************************************
*Export Authors' First-Stage Regression IV Tables as csv file*
**************************************************************

#delimit ;
esttab e2fs_o_add e5fs_o_add e6fs_o_add e7fs_o_add e8fs_o_add e9fs_o_add using "$project/results/TableFS_Original_Add.csv",
varlabels(persistence "Persistence of GDP process" cumulativeinstr "Original Instrument") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(persistence cumulativeinstr) 
order(persistence cumulativeinstr)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr

/**************************************************************
First Stage Results of Our Extension with Regional GDP Controls
***************************************************************/

*Column 2*
******************************************************

quietly: xtreg dem l(1/4).y l(1/4).yreg l(1/4).instrument yy* if samp2==1, fe cluster(wbcode2) r 
estimates store e2fs
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_with_yreg_d
estimates store e2fs_add

*Column 5*
******************************************************

quietly: xtreg dem l(1/4).y l(1/4).yreg l(1/4).instrument yy*  rtrend* if samp5==1, fe cluster(wbcode2) r 
estimates store e5fs
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_with_yreg_d
estimates store e5fs_add

*Column 6*
******************************************************

quietly: xtreg dem l(1/4).y l(1/4).yreg l(1/4).instrument l(1/4).(tradewbreg) yy* if samp6==1, fe cluster(wbcode2) r 
estimates store e6fs
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_with_yreg_d
estimates store e6fs_add

*Column 7*
******************************************************

quietly: xtreg dem l(1/4).y l(1/4).yreg l(1/4).instrument l(1/4).(tradewbreg unrestreg) yy* if samp7==1, fe cluster(wbcode2) r
estimates store e7fs
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_with_yreg_d
estimates store e7fs_add

*Column 8*
******************************************************

quietly: xtreg dem l(1/4).y l(1/4).yreg l(1/4).instrument l(1/4).yreg_w yy*  if samp8==1, fe cluster(wbcode2) r 
estimates store e8fs
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_with_yreg_d
estimates store e8fs_add

*Column 9*
******************************************************

quietly: xtreg dem l(1/4).y l(1/4).yreg l(1/4).instrument l(1/5).yreg_w l(1/4).dem_w yy*  if samp9==1, fe cluster(wbcode2) r 
estimates store e9fs
nlcom (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]) (instr_lag1: _b[L.instrument]) (instr_lag2: _b[L2.instrument]) (instr_lag3: _b[L3.instrument]) (instr_lag4: _b[L4.instrument]), post
vareffects_with_yreg_d
estimates store e9fs_add


*************************************************************
*Export Extension First-Stage Regression Tables as csv files*
*************************************************************

#delimit ;
esttab e2fs e5fs e6fs e7fs e8fs e9fs using "$project/results/TableFS.csv",
varlabels(L.instrument "Democracy wave t-1"  L2.instrument "Democracy wave t-2" L3.instrument "Democracy wave t-3" L4.instrument "Democracy wave t-4") 
cells(b(star fmt(%9.3f)) se(par))    
keep(L.instrument L2.instrument L3.instrument L4.instrument) 
order(L.instrument L2.instrument L3.instrument L4.instrument)
starlevels(* 0.10 ** 0.05 *** 0.01 **** 0.001)   
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none) ; 
#delimit cr

#delimit ;
esttab e2fs_add e5fs_add e6fs_add e7fs_add e8fs_add e9fs_add using "$project/results/TableFS_Add.csv",
varlabels(persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect" cumulativeinstr "Original Instrument") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(persistence cumulativereg cumulativeinstr) 
order(persistence cumulativereg cumulativeinstr)
starlevels(* 0.10 ** 0.05 *** 0.01 **** 0.001)  
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr



*******************************************************************
/******************************************************************
DIAGNOSTIC: Re-run 2SLS Extension Results by REGION
******************************************************************/
*******************************************************************


*******************************************************************
***AFRICA (EXCL. MENA)***
*******************************************************************


*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "AFR", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_AFR
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_AFR


*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "AFR", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_AFR
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_AFR

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "AFR", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_AFR
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_AFR


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "AFR", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_AFR
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_AFR


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "AFR", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_AFR
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_AFR


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "AFR", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_AFR
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_AFR


***********************************************************
*Export AFRICA (EXCL. MENA) Regression Tables as csv files*
***********************************************************

#delimit ;
esttab e2_AFR e5_AFR e6_AFR e7_AFR e8_AFR e9_AFR  using "$project/results/Table2SLSTop_AFR.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))     
keep(dem) order(dem)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_AFR e5_add_AFR e6_add_AFR e7_add_AFR e8_add_AFR e9_add_AFR using "$project/results/Table2SLSTop_Add_AFR.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg) 
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***EAST ASIA AND PACIFIC***
*******************************************************************

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "EAP", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_EAP
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_EAP


*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "EAP", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_EAP
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_EAP

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "EAP", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_EAP
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_EAP


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "EAP", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_EAP
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_EAP


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "EAP", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_EAP
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_EAP


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "EAP", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_EAP
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_EAP


*************************************************************
*Export EAST ASIA AND PACIFIC Regression Tables as csv files*
*************************************************************

#delimit ;
esttab e2_EAP e5_EAP e6_EAP e7_EAP e8_EAP e9_EAP  using "$project/results/Table2SLSTop_EAP.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))  
keep(dem) order(dem)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_EAP e5_add_EAP e6_add_EAP e7_add_EAP e8_add_EAP e9_add_EAP using 
"$project/results/Table2SLSTop_Add_EAP.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***EUROPE AND CENTRAL ASIA***
*******************************************************************

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "ECA", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_ECA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_ECA

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "ECA", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_ECA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_ECA

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "ECA", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_ECA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_ECA


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "ECA", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_ECA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_ECA


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "ECA", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_ECA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_ECA


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "ECA", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_ECA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_ECA


***************************************************************
*Export EUROPE AND CENTRAL ASIA Regression Tables as csv files*
***************************************************************

#delimit ;
esttab e2_ECA e5_ECA e6_ECA e7_ECA e8_ECA e9_ECA  using "$project/results/Table2SLSTop_ECA.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))   
keep(dem) order(dem)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_ECA e5_add_ECA e6_add_ECA e7_add_ECA e8_add_ECA e9_add_ECA using "$project/results/Table2SLSTop_Add_ECA.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***INDUSTRIALISED NATIONS***
*******************************************************************

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "INL", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_INL
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_INL

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "INL", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_INL
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_INL

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "INL", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_INL
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_INL


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "INL", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_INL
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_INL


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "INL", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_INL
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_INL


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "INL", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_INL
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_INL


**************************************************************
*Export INDUSTRIALISED NATIONS Regression Tables as csv files*
**************************************************************

#delimit ;
esttab e2_INL e5_INL e6_INL e7_INL e8_INL e9_INL  using "$project/results/Table2SLSTop_INL.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))    
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_INL e5_add_INL e6_add_INL e7_add_INL e8_add_INL e9_add_INL using 
"$project/results/Table2SLSTop_Add_INL.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***LATIN AMERICA & CARRIBEAN***
*******************************************************************

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_LAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_LAC

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_LAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_LAC

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_LAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_LAC

*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_LAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_LAC

*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_LAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_LAC

*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_LAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_LAC


*****************************************************************
*Export LATIN AMERICA & CARRIBEAN Regression Tables as csv files*
*****************************************************************

#delimit ;
esttab e2_LAC e5_LAC e6_LAC e7_LAC e8_LAC e9_LAC  using "$project/results/Table2SLSTop_LAC.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))    
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_LAC e5_add_LAC e6_add_LAC e7_add_LAC e8_add_LAC e9_add_LAC using "$project/results/Table2SLSTop_Add_LAC.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***MIDDLE EAST AND NORTH AFRICA***
*******************************************************************

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "MNA", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_MNA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_MNA

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "MNA", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_MNA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_MNA

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "MNA", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_MNA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_MNA


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "MNA", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_MNA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_MNA


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "MNA", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_MNA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_MNA


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "MNA", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_MNA
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_MNA


********************************************************************
*Export MIDDLE EAST AND NORTH AFRICA Regression Tables as csv files*
********************************************************************

#delimit ;
esttab e2_MNA e5_MNA e6_MNA e7_MNA e8_MNA e9_MNA  using "$project/results/Table2SLSTop_MNA.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))     
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_MNA e5_add_MNA e6_add_MNA e7_add_MNA e8_add_MNA e9_add_MNA using 
"$project/results/Table2SLSTop_Add_MNA.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr


*******************************************************************
***SOUTH ASIA***
*******************************************************************


*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region == "SAS", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_SAS
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_SAS

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region == "SAS", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_SAS
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_SAS

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region == "SAS", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_SAS
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_SAS


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region == "SAS", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_SAS
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_SAS


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region == "SAS", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_SAS
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_SAS


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region == "SAS", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_SAS
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_SAS


**************************************************
*Export SOUTH ASIA Regression Tables as csv files*
**************************************************

#delimit ;
esttab e2_SAS e5_SAS e6_SAS e7_SAS e8_SAS e9_SAS  using "$project/results/Table2SLSTop_SAS.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))     
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_SAS e5_add_SAS e6_add_SAS e7_add_SAS e8_add_SAS e9_add_SAS using "$project/results/Table2SLSTop_Add_SAS.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***GLOBAL EXCLUDING LAC***
*******************************************************************

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if region != "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_noLAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_noLAC

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if region != "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_noLAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_noLAC

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if region != "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_noLAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_noLAC


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if region != "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_noLAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_noLAC


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if region != "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_noLAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_noLAC


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if region != "LAC", fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_noLAC
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_noLAC


********************************************************
*Export GLOBAL EXCL. LAC Regression Tables as csv files*
********************************************************

#delimit ;
esttab e2_noLAC e5_noLAC e6_noLAC e7_noLAC e8_noLAC e9_noLAC  using "$project/results/Table2SLSTop_noLAC.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))     
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_noLAC e5_add_noLAC e6_add_noLAC e7_add_noLAC e8_add_noLAC e9_add_noLAC using "$project/results/Table2SLSTop_Add_noLAC.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


***ONLY EU REGRESSION (SIMPLE TIME-INVARYING DUMMY)***
*******************************************************************

***Create EU Dummy (EU=1 if country ever joins EU, EU = 0 if country never joins EU)***
capture drop EU
gen byte EU = 0


foreach C in "Slovenia" "Slovakia" "Romania" "Poland" "Latvia" "Lithuania" "Hungary" "Croatia" "Estonia" "Czech Republic" "Bulgaria" "Sweden" "Portugal" "Netherlands" "Luxembourg" "Italy" "Ireland" "Greece" "France" "Finland" "Spain" "Denmark" "Germany" "Cyprus" "Belgium" "Austria" "United Kingdom" {
    replace EU = 1 if country == "`C'" 
}

/* SPECIAL NOTES:
	1. We include United Kingdom here as dataset is from 1960-2010 while BREXIT was in 2020
	2. We exclude Malta as it is not present in the dataset
	3. Here, dummy does not segregate by when a country joins the EU, so Croatia included even
	   though they joined in 2013 and dataset ends in 2010
*/



*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* if EU == 1, fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_onlyEU
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_onlyEU

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy*  rtrend2-rtrend7 if EU == 1, fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_onlyEU
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_onlyEU

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* if EU == 1, fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_onlyEU
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_onlyEU


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* if EU == 1, fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_onlyEU
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_onlyEU


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w))  yy* if EU == 1, fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_onlyEU
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_onlyEU


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* if EU == 1, fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_onlyEU
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_onlyEU


*********************************************************************************
*Export ONLY EU REGRESSION (TIME-INVARYING DUMMY) Regression Tables as csv files*
*********************************************************************************

#delimit ;
esttab e2_onlyEU e5_onlyEU e6_onlyEU e7_onlyEU e8_onlyEU e9_onlyEU  using "$project/results/Table2SLSTop_onlyEU.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))     
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_onlyEU e5_add_onlyEU e6_add_onlyEU e7_add_onlyEU e8_add_onlyEU e9_add_onlyEU using "$project/results/Table2SLSTop_Add_onlyEU.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr


*******************************************************************
***GLOBAL REGRESSION WITH EU DUMMY CONTROLS***
*******************************************************************

/***Add EU Time-Varying Dummy (EU_t = 1 for country 'i' in all years 't' during which
	country is EU Member, otherwise EU_t = 0)***/

capture drop EU_t
gen byte EU_t = 0

* Original EU (Founding Members in 1958)
foreach c in "Belgium" "France" "Germany" "Italy" "Luxembourg" "Netherlands" {
    replace EU_t= 1 if country == "`c'" & year >= 1960  // Dataset starts in 1960
}

* 1973 enlargement
foreach c in "Denmark" "Ireland" "United Kingdom" {
    replace EU_t = 1 if country == "`c'" & year >= 1973
}

* 1981 enlargement
replace EU_t = 1 if country == "Greece" & year >= 1981

* 1986 enlargement
foreach c in "Spain" "Portugal" {
    replace EU_t = 1 if country == "`c'" & year >= 1986
}

* 1995 enlargement
foreach c in "Austria" "Finland" "Sweden" {
    replace EU_t = 1 if country == "`c'" & year >= 1995
}

* 2004 enlargement
foreach c in "Cyprus" "Czech Republic" "Estonia" "Hungary" "Latvia" "Lithuania" ///
             "Poland" "Slovakia" "Slovenia" {
    replace EU_t = 1 if country == "`c'" & year >= 2004
}

* 2007 enlargement
foreach c in "Bulgaria" "Romania" {
    replace EU_t = 1 if country == "`c'" & year >= 2007
}

/* SPECIAL NOTES:
	1. Croatia joined the EU in 2013 while the dataset ends at 2010, so Croatia excluded
	2. We exclude Malta as it is not present in the dataset
	3. BREXIT not reflected as it occurred in 2020 while the dataset ends at 2010
*/

*Column 2 (4 lags of GDP & Instrument)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* EU_t, fe cluster(wbcode2) r partial(yy*)
capture drop samp2
gen samp2=e(sample)
estimates store e2_EUControl
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e2_add_EUControl

*Column 5 (Add region-specific trends)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem=l(1/4).instrument) yy* EU_t rtrend2-rtrend7, fe cluster(wbcode2) r partial(yy*)
capture drop samp5
gen samp5=e(sample)
estimates store e5_EUControl
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e5_add_EUControl

*Column 6 (Add controls for regional trade)
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg=l(1/4).(instrument tradewbreg)) yy* EU_t, fe cluster(wbcode2) r partial(yy*)
capture drop samp6
gen samp6=e(sample)
estimates store e6_EUControl
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e6_add_EUControl


*Column 7 (Add control for regional unrest)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem tradewbreg unrestreg=l(1/4).(instrument tradewbreg unrestreg)) yy* EU_t, fe cluster(wbcode2) r partial(yy*)
capture drop samp7
gen samp7=e(sample)
estimates store e7_EUControl
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e7_add_EUControl


*Column 8 (Spatial Lags of GDP)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem yreg_w=l(1/4).(instrument yreg_w)) yy* EU_t, fe cluster(wbcode2) r partial(yy*)
capture drop samp8
gen samp8=e(sample)
estimates store e8_EUControl
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e8_add_EUControl


*Column 9 (Spatial Lags of GDP & Democracy)*
**********************************************************

quietly: xtivreg2 y l(1/4).y l(1/4).yreg (dem l(0/4).yreg_w dem_w=l(1/4).(instrument dem_w) l(1/5).yreg_w)  yy* EU_t, fe cluster(wbcode2) r partial(yy*)
capture drop samp9
gen samp9=e(sample)
estimates store e9_EUControl
nlcom (shortrun: _b[dem])  (lag1: _b[L.y])  (lag2: _b[L2.y])  (lag3: _b[L3.y])  (lag4: _b[L4.y]) (yreg_lag1: _b[L.yreg]) (yreg_lag2: _b[L2.yreg]) (yreg_lag3: _b[L3.yreg]) (yreg_lag4: _b[L4.yreg]), post
vareffects_with_yreg
estimates store e9_add_EUControl


********************************************************************************
*Export GLOBAL REGRESSION WITH EU DUMMY CONTROLS Regression Tables as csv files*
********************************************************************************

#delimit ;
esttab e2_EUControl e5_EUControl e6_EUControl e7_EUControl e8_EUControl e9_EUControl  using "$project/results/Table2SLSTop_EUControl.csv",
varlabels(dem "Democracy" ) 
cells(b(star fmt(%9.3f)) se(par)) 
stats( N N_g widstat jp, fmt(%7.0f %7.0f %7.1f %7.3f) labels( "Observations" "Countries in sample"  "Exc. Instruments F-stat." "Hansen p-value" ))     
keep(dem) order(dem) star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers 
collabels(none) ; 
#delimit cr

#delimit ;
esttab e2_add_EUControl e5_add_EUControl e6_add_EUControl e7_add_EUControl e8_add_EUControl e9_add_EUControl using "$project/results/Table2SLSTop_Add_EUControl.csv",
varlabels(longrun "Long-run effect of democracy"  effect$limit "Effect of democracy after $limit years"  persistence "Persistence of GDP process" cumulativereg "Cumulative Regional GDP Effect") 
cells(b(star fmt(%9.3f)) se(par)) 
keep(longrun effect$limit persistence cumulativereg) 
order(longrun effect$limit persistence cumulativereg)
star(* 0.10 ** 0.05 *** 0.01 **** 0.001)
nolabel replace mlabels("(2)" "(5)" "(6)" "(7)" "(8)" "(9)") nonumbers collabels(none); 
#delimit cr



/***IMPORTANT NOTE: When closing Stata console after running this do-file, do not save 
the workspace, as otherwise it may modify the original .dta file and cause issues
when running the do-file again***/
