*! version 1.1  Split-Panel Jackknife (SPJ) - Variable Name Fixed
program define xtlp, eclass sortpreserve
    version 14.0
    
    syntax varlist(min=2 numeric fv ts) [if] [in], [ FE TFE ] Method(string) [ Hor(numlist integer) YTRansf(string) SHock(numlist integer) Graph]

    * 1. Setup and Checks
    * -------------------    
    qui xtset
    local idvar `r(panelvar)'
    local timevar `r(timevar)'
    
    if "`idvar'" == "" {
        di as error "Error: You must xtset your data before using xtspj."
        exit 459
    }

    marksample touse
    markout `touse' `varlist'
	
	if "`method'" == "" {
        di as error "Error: option method() is required."
        exit 198
    }

    * 2. Pass data to Mata for processing
    * -----------------------------------
    gettoken depvar indepvars : varlist
	
	_fv_check_depvar `depvar'
	fvexpand `indepvars' if `touse'
	
	local clean_names "`r(varlist)'"
	
	* fe_type
	* --------
	local fe_type = .
    if "`tfe'" != "" {
        local fe_type = 2
    }
    if "`fe'"  != "" {
        if `fe_type' != . {
			di as error "Error: options fe and tfe cannot be specified together."
			exit 198
		}
		local fe_type = 1
    }
    if `fe_type' == . {
        local fe_type = 1
    }
	
	* method
	* -------
	local method = lower("`method'")
    local method_code = 0
    if "`method'" == "fe" {
        local method_code = 1
    }
    else if "`method'" == "spj" {
        local method_code = 2
    }
    else {
        di as error "Error: option method() must be either 'fe' or 'spj'."
        exit 198
    }
	
	* horizon
	* --------
	loc hors : subinstr local hor "," " ", all
	loc nh = wordcount("`hors'")
	
	if `nh' > 1 {
		* Case 1: User entered two numbers (e.g., "0 5" or "1 10")
		tokenize "`hors'"
		loc hs `1'   // Extract the first number as the start point
		loc hor `2'  // Extract the second number as the end point
	
		* Check: No more than 2 arguments allowed
		if `nh' > 2 {
			di as error "Error: Too many arguments in hor(). Please enter 'H' or 'Start End'."
			exit 198
		}
		* Core Constraint: Start horizon must be 0 or 1
		if `hs'!=0 & `hs'!=1 {
			di as error "Error: Start Horizon must be 0 or 1."
			exit 198
		}
		* Check: End horizon must be greater than start horizon
		if `hor' <= `hs' {
			di as error "Error: End horizon must be greater than start horizon."
			exit 198
		}
	
		loc hran `hs'/`hor'
	}
	else if `nh' == 1 {
		* Case 3: User entered only one number (e.g., "8")
		* Default start point set to 0 (i.e., 0 to 8)
		* Note: If you want single-digit input to start from 1 by default, change 0 to 1 below
		loc hs = 0 
		loc hor `hors'
		loc hran `hs'/`hor'
	}
	else if `nh' == 0 {
		* Case 2: User did not provide input (Default settings)
		* Default set to 0 to 0
		loc hs = 0
		loc hor = 0
		loc hran `hs'/`hor'
	}
	
	* dep variables
	* --------------
	if `hor' > 0 {
		* levels
		if "`ytransf'"==""|"`ytransf'"=="level" {
			local ytransf "level"
			forvalues h = `hran' {
				loc hstr = `h' - `hs'
				loc m = `h'
				tempvar y_h`hstr'
				qui gen `y_h`hstr'' = f`h'.`depvar'		
				loc trn`hstr' "`depvar'_h(`m')"
			}
			loc y y_h
		}
		* differences
		else if "`ytransf'"=="diff" {
			forvalues h = `hran' {
				loc hstr = `h' - `hs'
				loc m = `h'
				tempvar dy_h`hstr'
				qui gen `dy_h`hstr'' = f`h'.`depvar' - l.f`h'.`depvar' 
				loc trn`hstr' "D.`depvar'_h(`m')"
			}
			loc y dy_h
		}
		* Cumulative differences
		else if "`ytransf'"=="cmltdiff" {
			forvalues h = `hran' {
				loc hstr = `h' - `hs'
				loc m = `h'
				tempvar cdy_h`hstr'
				qui gen `cdy_h`hstr'' = f`h'.`depvar' - l.`depvar' 
				loc trn`hstr' "C.D.`depvar'_h(`m')"
			}
			loc y cdy_h
		}
		* Cumulative sum
		else if "`ytransf'"=="cmltsum" {
			forvalues h = `hran' {
				loc hstr = `h' - `hs'
				loc m = `h'
				tempvar cy_h`hstr'
				quietly gen double `cy_h`hstr'' = `depvar'
				if `h' > 0 {
					forvalues k = 1/`h' {
						quietly replace `cy_h`hstr'' = `cy_h`hstr'' + f`k'.`depvar'
					}
				}
				loc trn`hstr' "C.`depvar'_h(`m')"
			}
			loc y cy_h
		}
		else {
		    di as error "Error: option ytransf() must be level, diff, cmltdiff, or cmltsum."
			exit 198
		}
	}
	
	* shock variables
	* --------------
	local nshock = 1
    if "`shock'" != "" {
	    local nshock_nwords : word count `shock'
        if `nshock_nwords' != 1 {
            di as error "Error: option shock() must contain exactly one integer."
            exit 198
        }
        local nshock : word 1 of `shock'
		capture confirm integer number `nshock'
        if _rc {
            di as error "Error: shock() must be an integer."
            exit 198
        }
        local K : word count `clean_names'
        if (`nshock' < 1 | `nshock' > `K') {
            di as error "Error: shock() must be between 1 and `K' (number of covariates)."
            exit 198
        }
    }
	if `nshock' == 1 {
	    gettoken shockvar restvar : indepvars
	}
	else {
	    local shocklist ""
		forvalues i = 1/`nshock' {
			local shocklist `shocklist' `i'
		}
		local shockvarlist ""
		foreach j of local shocklist {
			local thisvar : word `j' of `clean_names'
			local shockvarlist `shockvarlist' `thisvar'
		}
		local Nshock : word count `shocklist'
	}
	
	* 3. Display
	* -----------
	if `method_code' == 1 {
		if `fe_type' == 1 {
			local title_txt "xtlp - FE - Individual Fixed Effects"
		}
		else {
			local title_txt "xtlp - FE - Two-way Fixed Effects (Individual + Time)"
		}
	}
	else {
		if `fe_type' == 1 {
			local title_txt "xtlp - SPJ - Individual Fixed Effects"
		}
		else {
			local title_txt "xtlp - SPJ - Two-way Fixed Effects (Individual + Time)"
		}
	}
	di _n as txt "`title_txt'"
	
	* 4. Estimate a single horizon or the full horizon range
	* -------------------------------------------------------
	if `hor' == 0 {
		
		mata: lp_work("`depvar'", "`indepvars'", "`touse'", "`idvar'", "`timevar'", `fe_type', `method_code')
		
		local N_val = scalar(N)
		di _n as txt _col(53) as txt "Number of obs =" ///
					 _col(67) as res %10.0fc `N_val'
		
// 		local depvar_trans_h "`trn`hstr''"
// 		local clean_names : subinstr local clean_names "`depvar'" "`depvar_trans_h'", word all
		
		matrix colnames b = `clean_names'
		matrix colnames V = `clean_names'
		matrix rownames V = `clean_names'
	
		ereturn post b V, esample(`touse')
	
		ereturn scalar N    = scalar(N)
		ereturn scalar N_g  = scalar(N_g)
		ereturn scalar df_r = scalar(df_r)
	
		ereturn local depvar   "depvar"
		ereturn local indepvars "`indepvars'"
		ereturn local cmd      "xtlp"
		ereturn local properties "b V"
		
		ereturn display
		
	}
	else if (`hor' != 0 & `nshock' == 1) {
		* plot data
		cap drop _birf _seirf _birf_lo _birf_up
		tempvar _t _zero birf seirf birf_up birf_lo 
		
		if `hs'<=0 loc h1 = `hor'+ 1 -`hs'
		else 	   loc h1 = `hor'
		if `hs'<=0 qui gen `_t' =_n-1+`hs'
		else  	   qui gen `_t' =_n
		qui gen `_zero' = 0
		
		qui gen `birf'    = 0 if _n<=`h1'
		qui gen `seirf'   = 0 if _n<=`h1'
		qui gen `birf_up' = 0 if _n<=`h1'
		qui gen `birf_lo' = 0 if _n<=`h1'
		
		* estimation		
		forval h=`hran' {
			if `hs'<=0 loc k=`h'+ 1 - `hs'
			else loc k=`h'
			loc hstr = `h' - `hs'
			
			local depvar_transf ``y'`hstr''
			tempvar touse_current
            qui gen byte `touse_current' = `touse'
            qui markout `touse_current' `depvar_transf' `indepvars'
            
			mata: lp_work("`depvar_transf'", "`indepvars'", "`touse_current'", "`idvar'", "`timevar'", `fe_type', `method_code')
			
			matrix colnames b = `clean_names'
			matrix colnames V = `clean_names'
			matrix rownames V = `clean_names'
			
			matrix b`h' = b
			matrix V`h' = V
			ereturn matrix b`h' b`h'
			ereturn matrix V`h' V`h'
			
			local coef = b[1, 1]
			local var  = V[1, 1]
			local se   = sqrt(`var')
			local ub = `coef' + 1.96 * `se'
			local lb = `coef' - 1.96 * `se'
			quietly {
				replace `birf'    = `coef' if _n == `k'
				replace `seirf'   = `se'  if _n == `k'
				replace `birf_up' = `ub'  if _n == `k'
				replace `birf_lo' = `lb'  if _n == `k'
			}
		}
		
		mkmat `birf'    if _n<=`h1', mat(BIRF)
		mkmat `seirf'   if _n<=`h1', mat(SEIRF)
		mkmat `birf_lo' if _n<=`h1', mat(SEIRF_LO)
		mkmat `birf_up' if _n<=`h1', mat(SEIRF_UP)
		
		mat IRF = BIRF , SEIRF , SEIRF_LO , SEIRF_UP
		matrix colnames IRF = "IRF" "Std.Err." "IRF LOW" "IRF UP"
		
		loc rows ""
		forval i=`hran' {
			loc rows `rows' `i'
		}
		matrix rownames IRF = `rows'		
		matlist IRF, noheader format(%9.5f) title("Impulse Response Function") lines(oneline)
		ereturn matrix irf IRF
		
		* graph
		loc mod = mod(`hor'-`hs',2)
		if `hor'-`hs'>12 & `mod'==0 loc p 2
		else if `hor'-`hs'>12 & `mod'==1 loc p 3
		else loc p 1
		
		if "`graph'"!="" {
		    if "`method'" == "fe" {
				loc lcolor blue
			}
			else if "`method'" == "spj" {
			    loc lcolor red
			}
			qui twoway (rarea `birf_up' `birf_lo' `_t', fcolor(`lcolor'%15) lc(`lcolor'%7)) ///
			(line `_zero' `_t', lcolor(gs5) lpattern(dash)) ///
			(line `birf' `_t', lcolor(`lcolor') lpattern(solid)) if _n<=`h1', ///
			legend(`off' order(3 "IRF of `y' (`depvar') to shock (`shockvar'), method(`method')") position(6)) tlabel(`hs'(`p')`hor') xtitle("Horizon") ///
			name("IRF_`method'", replace)
		}
	}
	else {
		* plot data
		cap drop _birf _seirf _birf_lo _birf_up
		tempvar _t _zero
		
		if `hs'<=0 loc h1 = `hor'+ 1 -`hs'
		else 	   loc h1 = `hor'
		if `hs'<=0 qui gen `_t' =_n-1+`hs'
		else  	   qui gen `_t' =_n
		qui gen `_zero' = 0
		
		tempname IRF_all
		
		local sidx = 0
		foreach idx of local shocklist {
			local ++sidx
			tempvar birf`sidx' seirf`sidx' birf_up`sidx' birf_lo`sidx'
			quietly gen `birf`sidx''    = 0 if _n <= `h1'
			quietly gen `seirf`sidx''   = 0 if _n <= `h1'
			quietly gen `birf_up`sidx'' = 0 if _n <= `h1'
			quietly gen `birf_lo`sidx'' = 0 if _n <= `h1'
		}
		
		local row_names ""
		
		* estimation		
		forval h=`hran' {
		    local row_names "`row_names' `h'"
			
			if `hs'<=0 loc k=`h'+ 1 - `hs'
			else loc k=`h'
			loc hstr = `h' - `hs'
			
			local depvar_transf ``y'`hstr''
			tempvar touse_current
            qui gen byte `touse_current' = `touse'
            qui markout `touse_current' `depvar_transf' `indepvars'
            
			mata: lp_work("`depvar_transf'", "`indepvars'", "`touse_current'", "`idvar'", "`timevar'", `fe_type', `method_code')
			
			matrix colnames b = `clean_names'
			matrix colnames V = `clean_names'
			matrix rownames V = `clean_names'
			
			matrix b`h' = b
			matrix V`h' = V
			ereturn matrix b`h' b`h'
			ereturn matrix V`h' V`h'
			
			local sidx = 0
			foreach idx of local shocklist {
				local ++sidx
				tempname coef var se ub lb
				scalar `coef' = b[1, `idx']
				scalar `var'  = V[`idx', `idx']
				scalar `se'   = sqrt(`var')
				scalar `ub'   = `coef' + 1.96 * `se'
				scalar `lb'   = `coef' - 1.96 * `se'
				quietly {
					replace `birf`sidx''    = `coef' if _n == `k'
					replace `seirf`sidx''   = `se'   if _n == `k'
					replace `birf_up`sidx'' = `ub'   if _n == `k'
					replace `birf_lo`sidx'' = `lb'   if _n == `k'
				}
			}
		}
		
		local sidx = 0
		local first_matrix = 1
		foreach idx of local shocklist {
			local ++sidx
			local shockvar : word `idx' of `clean_names'
			
			tempname BIRF SEIRF SEIRF_LO SEIRF_UP IRF_s
			
			mkmat `birf`sidx''    if _n <= `h1', mat(`BIRF')
			mkmat `seirf`sidx''   if _n <= `h1', mat(`SEIRF')
			mkmat `birf_lo`sidx'' if _n <= `h1', mat(`SEIRF_LO')
			mkmat `birf_up`sidx'' if _n <= `h1', mat(`SEIRF_UP')
			
			matrix `IRF_s' = `BIRF', `SEIRF', `SEIRF_LO', `SEIRF_UP'
			local cname1 "IRF_`sidx'"
			local cname2 "Std.Err._`sidx'"
			local cname3 "IRF LOW_`sidx'"
			local cname4 "IRF UP_`sidx'"
			matrix colnames `IRF_s' = "`cname1'" "`cname2'" "`cname3'" "`cname4'"
			matrix rownames `IRF_s' = `row_names'
			
			if `first_matrix' == 1 {
                matrix `IRF_all' = `IRF_s'
                local first_matrix = 0
            }
            else {
                matrix `IRF_all' = `IRF_all', `IRF_s'
            }
	
			matlist `IRF_s', noheader format(%9.5f) ///
				title("Impulse Response Function for shock #`sidx' (`shockvar')") ///
				lines(oneline)
			
			* graph
			if "`graph'" != "" {
				local mod = mod(`hor' - `hs', 2)
				if `hor' - `hs' > 12 & `mod' == 0 local p = 2
				else if `hor' - `hs' > 12 & `mod' == 1 local p = 3
				else local p = 1
				
				if "`method'" == "fe" {
					loc lcolor blue
				}
				else if "`method'" == "spj" {
					loc lcolor red
				}
				quietly twoway ///
					(rarea `birf_up`sidx'' `birf_lo`sidx'' `_t', fcolor(`lcolor'%15) lc(`lcolor'%7)) ///
					(line `_zero' `_t', lcolor(gs5) lpattern(dash)) ///
					(line `birf`sidx'' `_t', lcolor(`lcolor') lpattern(solid)) ///
					if _n <= `h1', ///
					legend(`off' order(3 "IRF of `y' (`depvar') to shock (`shockvar'), method(`method')") position(6)) ///
					tlabel(`hs'(`p')`hor') ///
					xtitle("Horizon") ///
					name("IRF_`method'_`sidx'", replace)
			}
		}
		matrix rownames `IRF_all' = `row_names'
		ereturn matrix irf `IRF_all'
	}
end

* -----------------------------------------------------------------------------
* MATA CODE BLOCK
* -----------------------------------------------------------------------------
version 14.0
set matalnum on
mata:

void lp_work(string scalar depvar, ///
			 string scalar indepvars, ///
			 string scalar touse, ///
             string scalar idvar, ///
			 string scalar timevar, ///
             real   scalar fe_type, ///
             real   scalar method_code)
{
    real matrix Y, X, info, V_est, b_est
    real colvector ID, Time
    real scalar N, K, N_g
    string rowvector xnames
    
    // Load Data
    Y = st_data(., depvar, touse)
    X = st_data(., indepvars, touse)
    ID = st_data(., idvar, touse)
    Time = st_data(., timevar, touse)
    
    // Get variable names for labeling
    xnames = tokens(indepvars)
    
    // Panel Info setup (sort by ID, Time)
    // We need to sort X, Y, ID, Time together to ensure structure
    real matrix ALL
    ALL = Y, ID, Time, X
    _sort(ALL, (2,3)) // Sort by ID then Time
    
    Y = ALL[., 1]
    ID = ALL[., 2]
    Time = ALL[., 3]
    X = ALL[., 4::cols(ALL)]
    
    N = rows(Y)
    K = cols(X)
    
    // Panel Setup info: [start_index, end_index] for each individual
    info = panelsetup(ID, 1)
    N_g = rows(info) // Number of groups (individuals)

    // -------------------------------------------------------
    // Step 1: Define Split Points (The Logic provided)
    // -------------------------------------------------------
    real colvector T_a_idx, T_b_idx, cut_i
    T_a_idx = J(N, 1, 0)
    T_b_idx = J(N, 1, 0)
    cut_i   = J(N_g, 1, .)
	
    real scalar i, start, end_t, Ti, cut_idx
   
    for (i=1; i<=N_g; i++) {
        start = info[i, 1]
        end_t = info[i, 2]
        Ti = end_t - start + 1
      
        cut_idx = floor((Ti + 1) / 2)
		cut_i[i] = cut_idx
		
        // Mark Part A (1 to cut)
        if (cut_idx >= 1) {
            T_a_idx[|start \ (start + cut_idx - 1)|] = J(cut_idx, 1, 1)
        }
        // Mark Part B (cut+1 to Ti)
        if (cut_idx < Ti) {
            T_b_idx[|(start + cut_idx) \ end_t|] = J(Ti - cut_idx, 1, 1)
        }
    }

    // -------------------------------------------------------
    // Step 2: Within Transformation (Demeaning)
    // -------------------------------------------------------
	// Full sample
	real matrix YX, YX_dm
    real colvector Y_dot
    real matrix    X_dot
    YX    = Y, X
    YX_dm = twoway_demean(YX, ID, Time, info, fe_type)
    Y_dot = YX_dm[., 1]
    X_dot = YX_dm[., 2..(K+1)]
	
	// split
	real colvector selA, selB
    selA = selectindex(T_a_idx :!= 0)
    selB = selectindex(T_b_idx :!= 0)

    // A sample
    real colvector Y_A, ID_A, Time_A
    real matrix   X_A
    real matrix   infoA
    Y_A    = Y[selA, .]
    X_A    = X[selA, .]
    ID_A   = ID[selA, .]
    Time_A = Time[selA, .]
    infoA  = panelsetup(ID_A, 1)
	
	real matrix YX_a, YX_a_dm
	real colvector Y_a_dot
    real matrix   X_a_dot
    YX_a    = Y_A, X_A
    YX_a_dm = twoway_demean(YX_a, ID_A, Time_A, infoA, fe_type)
    Y_a_dot = YX_a_dm[., 1]
    X_a_dot = YX_a_dm[., 2..(cols(YX_a_dm))]

    // B sample
    real colvector Y_B, ID_B, Time_B
    real matrix   X_B
    real matrix   infoB
    Y_B    = Y[selB, .]
    X_B    = X[selB, .]
    ID_B   = ID[selB, .]
    Time_B = Time[selB, .]
    infoB  = panelsetup(ID_B, 1)
	
	real matrix YX_b, YX_b_dm
    real colvector Y_b_dot
    real matrix    X_b_dot
    YX_b    = Y_B, X_B
    YX_b_dm = twoway_demean(YX_b, ID_B, Time_B, infoB, fe_type)
    Y_b_dot = YX_b_dm[., 1]
    X_b_dot = YX_b_dm[., 2..(cols(YX_b_dm))]
	
	// prepare for d_dot
	real matrix X_dot_a_full, X_dot_b_full
	real scalar NT
	NT = rows(X_dot)
	
	X_dot_a_full = J(NT, K, .)
	X_dot_b_full = J(NT, K, .)
	
	X_dot_a_full[selA, .] = X_a_dot
	X_dot_b_full[selB, .] = X_b_dot

	// -------------------------------------------------------
    // Step 3: Estimate Coefficients (OLS on demeaned data)
    // -------------------------------------------------------
    real colvector b_full, b_a, b_b
    real matrix XX_inv_full, XX_inv_a, XX_inv_b
  
    // Full
    XX_inv_full = cholinv(cross(X_dot, X_dot))
    b_full = XX_inv_full * cross(X_dot, Y_dot)
  
    // Part A
    XX_inv_a = cholinv(cross(X_a_dot, X_a_dot))
    b_a = XX_inv_a * cross(X_a_dot, Y_a_dot)
  
    // Part B
    XX_inv_b = cholinv(cross(X_b_dot, X_b_dot))
    b_b = XX_inv_b * cross(X_b_dot, Y_b_dot)
  
    // -------------------------------------------------------
    // Step 4: Estimator
    // -------------------------------------------------------
	if (method_code == 1) {
	    b_est = b_full
	}
	else if (method_code == 2) {
	    b_est = 2 * b_full - 0.5 * (b_a + b_b)
	}
  
    // -------------------------------------------------------
    // Step 5: Variance Calculation
    // -------------------------------------------------------
	real matrix X_mat, X_a_mat, X_b_mat, X_sub, d_dot
	X_mat  = X_dot
	X_a_mat= X_dot_a_full
	X_b_mat= X_dot_b_full
	
	X_sub  = J(NT, K, .)
	d_dot  = J(NT, K, .)
	for (i=1; i<=N_g; i++) {
		real scalar start1, end_t1, Ti1, ci
		start1 = info[i,1]
		end_t1 = info[i,2]
		Ti1    = end_t1 - start1 + 1
		ci    = cut_i[i]
		if (ci != .) {
			if (ci > 0) {
				X_sub[| start1,1 \ start1+ci-1,K |] = ///
					X_a_mat[| start1,1 \ start1+ci-1,K |]
			}
			if (ci < Ti1) {
				X_sub[| start1+ci,1 \ end_t1,K |] = ///
					X_b_mat[| start1+ci,1 \ end_t1,K |]
			}
		}
	}
	if (method_code == 1) {
	    d_dot = X_mat
	}
	else if (method_code == 2) {
	    d_dot = 2:*X_mat :- X_sub
	}
	
	real colvector e
	e = Y_dot - X_dot*b_est      // NT x 1

	real matrix W_N
	W_N = J(K, K, 0)
	for (i = 1; i <= N_g; i++) {
		real scalar start2, end_t2, Ti2
		start2 = info[i,1]
		end_t2 = info[i,2]
		Ti2    = end_t2 - start2 + 1
	
		real matrix d_i
		d_i = d_dot[| start2, 1 \ end_t2, K |]      // Ti x K
	
		real colvector e_i
		e_i = e[| start2, 1 \ end_t2, 1 |]          // Ti x 1
	
		real colvector g_i
		g_i = d_i' * e_i                        // K x 1
	
		W_N = W_N + g_i * g_i'                  // K x K
	}

	// Sandwich: var0_mat = XX_inv * W * XX_inv
// 	real scalar adj
// 	adj = (N_g / (N_g - 1.0)) * ((NT - 1.0) / (NT - K))
	real matrix W
// 	W = W_N * adj
	W = W_N
	real matrix XX_inv
	XX_inv = cholinv(cross(X_dot, X_dot))
	V_est = XX_inv * W * XX_inv
	real colvector se
	se = sqrt(diagonal(V))
	
    // -------------------------------------------------------
    // Step 6: Post Results to Stata
    // -------------------------------------------------------
	st_matrix("b", b_est')    // 1 x K
	st_matrix("V", V_est)     // K x K
	
	st_numscalar("N",   N)
	st_numscalar("N_g", N_g)
	st_numscalar("df_r", N_g - 1)
}

// function: twoway demean
real matrix twoway_demean(real matrix Z, ///
						  real colvector ID, ///
                          real colvector Time, ///
						  real matrix info, ///
                          real scalar fe_type)
{
    real scalar i, start, end_t
    real matrix Z_dm

    Z_dm = Z

    // -------- 1. individual FE (within-id) --------
    for (i = 1; i <= rows(info); i++) {
        start = info[i,1]
        end_t = info[i,2]

        real matrix Zi
        real rowvector mean_i

        Zi     = Z[|start,1 \ end_t,.|]
        mean_i = mean(Zi)
        Z_dm[|start,1 \ end_t,.|] = Zi :- mean_i
    }

    if (fe_type == 1) {
		return(Z_dm)
	}
	else {
		// -------- 2. time FE (within-time) on top of FE --------
		real colvector tuniq
		tuniq = uniqrows(Time)
		real scalar Nt, j, tval
		Nt = rows(tuniq)
		
		for (j = 1; j <= Nt; j++) {
			tval = tuniq[j]
		
			real colvector idx_t
			idx_t = selectindex(Time :== tval)
			if (rows(idx_t) == 0) continue
		
			real matrix Zt
			Zt = Z_dm[idx_t, .]
			real rowvector mean_t
			mean_t = mean(Zt)
		
			Z_dm[idx_t, .] = Zt :- mean_t
		}
		
		return(Z_dm)
	}
}

end