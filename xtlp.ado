*! version 1.1.0  xtlp
program define xtlp, eclass sortpreserve
    version 14.0
    
    * Extract vce(...) verbatim from the command line before syntax().
    * Stata's syntax varlist ... , VCE(...) chokes when vce() contains a comma
    * (e.g. vce(dkraay lag(2), ase)) because varlist tries to absorb
    * the inner token.
    * (reghdfe avoids this via [*] self-parsing; xtreg via _vce_parserun.)
    *
    * We locate "vce(" (case-insensitive) and scan for the matching ")" by
    * counting nested parentheses. A single regex cannot do this correctly:
    *   - greedy  vce\((.*)\)    over-matches: vce(dkraay) hor(0) -> "dkraay) hor(0"
    *   - lazy    vce\(.*?\)     under-matches: vce(dkraay lag(2), ase) -> "dkraay lag(2"
    * so we do a manual bracket-pairing scan instead.
    local cmdline : copy local 0
    local vce ""
    local vce_seen 0
    local vce0 `"`0'"'
    local vce0l = lower(`"`vce0'"')
    local vce_pos = strpos(`"`vce0l'"', "vce(")
    if `vce_pos' > 0 {
        local vce_seen 1
        local open_pos = `vce_pos' + 3            // index of "(" in "vce("
        local n = strlen(`"`vce0'"')
        local depth = 1
        local i = `open_pos' + 1
        local close_pos = 0
        while (`i' <= `n') & (`close_pos' == 0) {
            local ch = substr(`"`vce0'"', `i', 1)
            if "`ch'" == "(" local depth = `depth' + 1
            else if "`ch'" == ")" {
                local depth = `depth' - 1
                if `depth' == 0 local close_pos = `i'
            }
            local i = `i' + 1
        }
        if `close_pos' == 0 {
            di as error "Error: vce() has unbalanced parentheses."
            exit 198
        }
        * content between the parentheses
        local vce = substr(`"`vce0'"', `open_pos' + 1, `close_pos' - `open_pos' - 1)
        * reject a second vce() option
        local vce_rest = substr(`"`vce0l'"', `close_pos' + 1, `n' - `close_pos')
        if strpos(`"`vce_rest'"', "vce(") > 0 {
            di as error "Error: vce() specified more than once."
            exit 198
        }
        * remove "vce(...)" from `0' so syntax() does not see it
        local vce_full = substr(`"`vce0'"', `vce_pos', `close_pos' - `vce_pos' + 1)
        local 0 = subinword(`"`vce0'"', `"`vce_full'"', "", 1)
    }

    * vce() was present but empty (e.g. "vce()") -> clear error, not silent default
    if `vce_seen' & `"`vce'"' == "" {
        di as error "Error: vce() specified with no type; use one of unadjusted, robust, cluster, or dkraay."
        exit 198
    }

    syntax varlist(min=2 numeric fv ts) [if] [in], [ FE TFE KEEPSINgletons ] ///
		Method(string) [ Hor(numlist integer) YTRansf(string) SHock(numlist integer) ///
		Graph]

	**# 1. Setup and Checks
    * ----------------------
    qui xtset
    local idvar `r(panelvar)'
    local timevar `r(timevar)'
    
    if "`idvar'" == "" | "`timevar'" == "" {
		di as error "Error: You must xtset your data with both panel and time variables (e.g., xtset id time)."
        exit 459
    }
	
	if "`method'" == "" {
        di as error "Error: option method() is required."
        exit 198
    }

	**# 2. Pass data to Mata for processing
    * --------------------------------------
    gettoken depvar indepvars : varlist
	
	marksample touse, novarlist
	
	_fv_check_depvar `depvar'
	markout `touse' `idvar' `timevar'
	fvexpand `indepvars' if `touse'
	local indepvars_list "`r(varlist)'"
	
	markout `touse' `indepvars_list' `idvar' `timevar'
	
	**# fe_type
	* ----------
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

	**# singleton groups (Correia 2015; aligned to reghdfe default)
	* ---------------------------------------------------------------
	* By default (keepsingletons not specified) drop singleton groups, on the
	* SAME full regression sample reghdfe uses (depvar + indepvars marked out).
	* Singletons contribute 0 to the within-demeaned moment but inflate N and
	* the FE dof, biasing SEs downward; reghdfe drops them by default.
	* keepsingletons keeps them and warns (same text as reghdfe).
	*
	* NOTE: the drop is applied PER HORIZON (each horizon's lead differs, so
	* its sample -- and thus its singletons -- can differ). Each horizon that
	* drops >0 prints its own line tagged with h:
	* "(h=#: dropped # singleton observations)".
	local _sing_drop = ("`keepsingletons'" == "")
	if "`keepsingletons'" != "" {
		di as error `"WARNING: Singleton observations not dropped; statistical significance is biased"'
    }
	
	**# method
	* ---------
	* method() is case-sensitive (lowercase only), matching xtreg/reghdfe.
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

	**# vce
	* ------
	* vce_code: 0 = unadjusted (homoskedastic, default),
	*           1 = robust HC1 (heteroskedasticity-only),
	*           2 = one-way cluster,
	*           3 = two-way cluster (Cameron-Gelbach-Miller 2011),
	*           4 = Driscoll-Kraay
	* Syntax: the comma separates the vcetype (with its type-specific options)
	* from the generic SE suboptions. Before the comma, dkraay accepts only
	* lag(#). After the comma, "ase" and "nodfadj" are allowed (mutually
	* exclusive; nodfadj only with dkraay).
	*   vce(unadjusted [, ase])                unadjusted (homoskedastic); aligned to xtreg/reghdfe
	*   vce(robust [, ase])                    HC1 robust (heteroskedasticity-only); aligned to reghdfe
	*   vce(cluster clustervar [clustervar2] [, ase])   one- or two-way cluster;
	*           clustervar REQUIRED; clustervar2 REQUIRED for two-way (aligned to reghdfe)
	*   vce(dkraay [lag(#)] [, ase | nodfadj]) Driscoll-Kraay; lag omitted/-1 = automatic bandwidth
	*           (lag(#) is dkraay-specific and goes BEFORE the comma; ase and nodfadj
	*            are generic SE suboptions that go AFTER the comma, and are MUTUALLY
	*            EXCLUSIVE -- ase drops all small-sample brackets, nodfadj keeps the
	*            small-sample adj but does not subtract df_a)
	*   nodfadj = do NOT subtract df_a (absorbed FE dof) in the small-sample adjustment;
	*            i.e. use the xtscc-style factor (default subtracts df_a, aligned to reghdfe ivreg2-style).
	*            Only allowed with dkraay; not allowed together with ase.
	*   ase = no small-sample adjustment (asymptotic SE), aligned to xtscc `ase` / reghdfe `asymptotic`.
	*         Not allowed together with nodfadj.
	* vcetype abbreviations (each type accepts only its short form + full name):
	*   un/unadjusted -> unadjusted ; r/robust ; cl/cluster ; dk/dkraay
	local vce_code = 0
	local dk_lag = -1
	local ase_code = 0
	local nodfadj_code = 0

	* Parse the vce() content (already rejected empty vce() above).
	if `"`vce'"' != "" {
		* Split "main [, subopts]" (substr needs explicit length here)
		local vce_main  `"`vce'"'
		local vce_subopts ""
		local vce_comma = strpos(`"`vce'"', ",")
		local vce_len   = strlen(`"`vce'"')
		if `vce_comma' > 0 {
			local vce_main   = substr(`"`vce'"', 1, `vce_comma' - 1)
			local vce_subopts = substr(`"`vce'"', `vce_comma' + 1, `vce_len' - `vce_comma')
		}

		* Parse suboptions after the comma. Two generic SE suboptions are allowed
		* here: "ase" and "nodfadj". nodfadj is only meaningful with dkraay;
		* that restriction is enforced in each vcetype branch below (non-dkraay
		* types reject nodfadj first). The ase/nodfadj mutual-exclusion check
		* runs inside the dkraay branch (after the "only with dkraay" gate).
		* dkraay-specific lag(#) goes BEFORE the comma.
		foreach sub of local vce_subopts {
			if `"`sub'"' == "ase" {
				if `ase_code' == 1 {
					di as error "Error: vce() suboption ase specified more than once."
					exit 198
				}
				local ase_code = 1
			}
			else if `"`sub'"' == "nodfadj" {
				if `nodfadj_code' == 1 {
					di as error "Error: vce() suboption nodfadj specified more than once."
					exit 198
				}
				local nodfadj_code = 1
			}
			else {
				di as error "Error: vce() suboption `sub' not recognized; only ase and nodfadj are allowed after the comma."
				exit 198
			}
		}

		tokenize `"`vce_main'"'
		local vcetype `"`1'"'
		* vcetype is case-sensitive (lowercase only), matching xtreg/reghdfe.

		* vcetype abbreviations: each type accepts only its short form + full:
		*   un -> unadjusted ; r -> robust ; cl -> cluster ; dk -> dkraay
		if `"`vcetype'"' == "un" {
			local vcetype "unadjusted"
		}
		else if `"`vcetype'"' == "r" {
			local vcetype "robust"
		}
		else if `"`vcetype'"' == "cl" {
			local vcetype "cluster"
		}
		else if `"`vcetype'"' == "dk" {
			local vcetype "dkraay"
		}

		if `"`vcetype'"' == "unadjusted" {
			local vce_code = 0
			* unadjusted takes NO positional argument
			if `"`2'"' != "" {
				di as error "Error: vce(unadjusted) takes no argument; got: `2'"
				exit 198
			}
			* nodfadj only applies to dkraay; reject elsewhere to catch user mistakes
			if `nodfadj_code' == 1 {
				di as error "Error: vce() suboption nodfadj only allowed with dkraay."
				exit 198
			}
		}
		else if `"`vcetype'"' == "robust" {
			local vce_code = 1
			* robust takes NO positional argument
			if `"`2'"' != "" {
				di as error "Error: vce(robust) takes no argument; got: `2'"
				exit 198
			}
			if `nodfadj_code' == 1 {
				di as error "Error: vce() suboption nodfadj only allowed with dkraay."
				exit 198
			}
		}
		else if `"`vcetype'"' == "cluster" {
			* cluster takes up to TWO clustervars; either may be any existing
			* variable (aligned to reghdfe: the cluster dimensions are independent
			* of the absorbed FE). clustervar is REQUIRED; clustervar2 is REQUIRED 
			* for two-way.
			if `nodfadj_code' == 1 {
				di as error "Error: vce() suboption nodfadj only allowed with dkraay."
				exit 198
			}
			if `"`4'"' != "" {
				di as error "Error: vce(cluster) takes at most two clustervars; got: `2' `3' `4'"
				exit 198
			}
			if `"`2'"' == "" {
				di as error "Error: vce(cluster) requires at least one clustervar."
				exit 198
			}
			local clustervar `"`2'"'
			confirm variable `clustervar'
			if `"`3'"' != "" {
				* Two-way cluster (clustervar x clustervar2)
				local clustervar2 `"`3'"'
				confirm variable `clustervar2'
				local vce_code = 3
			}
			else {
				local vce_code = 2
			}
		}
		else if `"`vcetype'"' == "dkraay" {
			local vce_code = 4
			* dkraay's only pre-comma option is lag(#). nodfadj and ase go AFTER
			* the comma (parsed above into nodfadj_code/ase_code). So before the
			* comma, token 2 may be lag(#) and token 3 must be empty.
				if `"`3'"' != "" {
				di as error "Error: vce(dkraay) accepts only lag(#) before the comma; got: `2' `3'"
					exit 198
				}
			if `"`2'"' != "" {
				* token 2 must be lag(#) with # a non-negative integer (>=0).
				* lag(-1) is NOT accepted: -1 is only an internal sentinel for
				* "automatic bandwidth" (used when the user omits lag(#) entirely,
				* aligned to xtscc/reghdfe default). A negative user lag is rejected
				* (aligned to reghdfe, which rejects negative bandwidths).
				if !regexm(`"`2'"', "^lag\(([0-9]+)\)$") {
					di as error "Error: vce(dkraay) expects lag(#) with # a non-negative integer before the comma; got: `2'"
					exit 198
				}
				local dk_lag = regexs(1)
			}
			* ase and nodfadj are mutually exclusive: ase drops ALL small-sample
			* brackets (so the df_a convention is moot), nodfadj keeps them but
			* switches the df_a convention. Allowing both is contradictory.
			if `ase_code' == 1 & `nodfadj_code' == 1 {
				di as error "Error: vce(dkraay) suboptions ase and nodfadj are mutually exclusive."
					exit 198
			}
		}
		else {
			di as error "Error: vce() type must be unadjusted, robust, cluster or dkraay; got: `vcetype'"
			exit 198
		}
	}

	* Build canonical vce() string and display vcetype for e()
	if `vce_code' == 0 {
		local vce_string "unadjusted"
		local vcetype_string ""
	}
	else if `vce_code' == 1 {
		local vce_string "robust"
		local vcetype_string "Robust"
	}
	else if `vce_code' == 2 {
		local vce_string "cluster `clustervar'"
		local vcetype_string "Cluster"
	}
	else if `vce_code' == 3 {
		local vce_string "cluster `clustervar' `clustervar2'"
		local vcetype_string "Cluster"
	}
	else if `vce_code' == 4 {
		local vce_string "dkraay"
		if `dk_lag' >= 0 local vce_string "`vce_string' lag(`dk_lag')"
		if `nodfadj_code' == 1 local vce_string "`vce_string' nodfadj"
		local vcetype_string "Drisc/Kraay"
	}
	else {
		local vce_string ""
		local vcetype_string ""
	}

	**# horizon
	* ----------
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
	
	**# dep variables
	* ----------------
	if `hor' > 0 {
		* levels
		if "`ytransf'"=="" | "`ytransf'"=="level" {
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
	
	**# shock variables
	* ------------------
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
        local K : word count `indepvars_list'
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
			local thisvar : word `j' of `indepvars_list'
			local shockvarlist `shockvarlist' `thisvar'
		}
		local Nshock : word count `shocklist'
	}
	
	**# 3. Display
	* -------------
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
	
	**# 4. Estimate a single horizon or full horizon range
	* -------------------------------------------------------
	if `hor' == 0 {
		**## h=0

		* graph is only meaningful for multi-horizon IRFs (hor>0); at hor=0
		* there is a single point to plot, so warn and ignore the option.
		if "`graph'" != "" {
			di as txt "note: graph option ignored (only valid with hor>0, i.e. a multi-horizon IRF)."
		}

		markout `touse' `depvar' `indepvars_list' `idvar' `timevar'

		* Drop singletons on the FULL regression sample (depvar marked out),
		* matching reghdfe; otherwise id-singletons whose only non-missing row
		* is missing the depvar are miscounted.
		if `_sing_drop' {
			_xtlp_drop_singletons, touse(`touse') idvar(`idvar') timevar(`timevar') fetype(`fe_type')
			local _hdropped = r(dropped)
			if `_hdropped' > 0 {
				di as txt `"(dropped `_hdropped' singleton observations)"'
			}
		}
		
		mata: st_numscalar("XTLP_H", 0)
		mata: lp_work("`depvar'", "`indepvars_list'", "`touse'", "`idvar'", "`timevar'", `fe_type', `method_code', `vce_code', `dk_lag', `ase_code', `nodfadj_code', "`clustervar'", "`clustervar2'", 1)
		
		**### header
		local N_val   = scalar(N)
		local T_eff   = scalar(df_r) + 1
		local Ng      = scalar(N_g)
		// e(lag)/scalar(lag) exists only for dkraay (vce_code==4); read it
		// only then so the header does not error on other vce types.
		if `vce_code' == 4 {
			local lag_val = scalar(lag)
		}
		else {
			local lag_val = .
		}

		if `method_code' == 1 {
			di _n as txt "Method: " as result "Fixed effects" ///
					 _col(50) as txt "Number of obs" ///
					 _col(68) as txt "=" as res %10.0fc `N_val'
		}
		else {
			di _n as txt "Method: " as result "Split-panel jackknife" ///
					 _col(50) as txt "Number of obs" ///
					 _col(68) as txt "=" as res %10.0fc `N_val'
		}

		if `vce_code' == 0 {
			di as txt "Standard errors: " as result "Unadjusted"
		}
		else if `vce_code' == 1 {
			di as txt "Standard errors: " as result "Robust"
		}
		else if `vce_code' == 2 {
			di as txt "Standard errors: " as result "One-way clustered"
		}
		else if `vce_code' == 3 {
			di as txt "Standard errors: " as result "Two-way clustered"
		}
		else if `vce_code' == 4 {
			if `nodfadj_code' == 1 {
				di as txt "Standard errors: " as result "Driscoll-Kraay (nodfadj)"
			}
			else {
				di as txt "Standard errors: " as result "Driscoll-Kraay"
			}
			if `lag_val' < . {
				di as txt "Maximum lag: " as res `lag_val'
			}
		}

		if `ase_code' == 0 {
			di as txt "Small-sample adjustment: " as result "On" _n
		}
		else if `ase_code' == 1 {
			di as txt "Small-sample adjustment: " as result "Off" _n
		}

		if `vce_code' == 2 {
			local Nclust = scalar(N_clust_tw)
			local msg "(Std. err. adjusted for `Nclust' clusters in `clustervar')"
			local len = strlen("`msg'")
			local start_col = 79 - `len'
			di as txt _col(`start_col') "(Std. err. adjusted for " ///
			   as res `Nclust' as txt " clusters in " as res "`clustervar'" as txt ")"
		}
		else if `vce_code' == 3 {
			local Nclust = scalar(N_clust_tw)
			local msg "(Std. err. adjusted for `Nclust' clusters in `clustervar' `clustervar2')"
			local len = strlen("`msg'")
			local start_col = 79 - `len'
			di as txt _col(`start_col') "(Std. err. adjusted for " ///
			   as res `Nclust' as txt " clusters in " as res "`clustervar'" as txt " " ///
			   as res "`clustervar2'" as txt ")"
		}
		else if `vce_code' == 4 {
			local msg "(Std. err. adjusted for `Ng' clusters in `timevar')"
			local len = strlen("`msg'")
			local start_col = 79 - `len'
			di as txt _col(`start_col') "(Std. err. adjusted for " ///
			   as res `T_eff' as txt " clusters in " as res "`timevar'" as txt ")"
		}
		
		**### ereturn
		matrix colnames b = `indepvars_list'
		matrix colnames V = `indepvars_list'
		matrix rownames V = `indepvars_list'
	
		ereturn post b V, esample(`touse')
	
		ereturn scalar N      = scalar(N)
		ereturn scalar N_g    = scalar(N_g)
		ereturn scalar df_r   = scalar(df_r)
		if `vce_code' == 4 {
			ereturn scalar lag    = scalar(lag)   // e(lag) only exists for dkraay
		}
		if `vce_code' == 2 | `vce_code' == 3 {
			ereturn scalar N_clust = scalar(N_clust_tw)   // # clusters; aligned to reghdfe e(N_clust)
		}
	
		ereturn local depvar     "`depvar'"
		ereturn local indepvars  "`indepvars'"
		ereturn local vce        "`vce_string'"
		ereturn local vcetype    "`vcetype_string'"
		ereturn local cmd        "xtlp"
		ereturn local cmdline    `"xtlp `cmdline'"'
		ereturn local properties "b V"
		
		ereturn display
		
	}
	else if (`hor' != 0 & `nshock' == 1) {
		**## h>0 shock=1
		
		* plot data
		cap drop _birf _seirf _birf_lo _birf_up
		tempvar _t _zero birf seirf birf_up birf_lo 
		
        if `hs'<=0 loc h1 = `hor' + 1 - `hs'
        else       loc h1 = `hor'

        * Cap the stored IRF rows at the observation count: the IRF is laid out
        * one row per observation (_n==k maps horizon h to row k), so when the
        * requested horizon range exceeds _N there are no rows to hold the extra
        * horizons. Capping avoids an mkmat/rownames conformability error; the
        * over-range horizons are still estimated and saved to e(b`h')/e(V`h')
        * but cannot appear in the IRF table/graph.
        local nobs = _N
        if `h1' > `nobs' {
            di as txt "note: requested horizon range (`h1' rows) exceeds the number of observations (" `nobs' "); IRF table truncated to " `nobs' " rows."
            loc h1 = `nobs'
        }

        if `hs'<=0 qui gen `_t' = _n - 1 + `hs'
        else       qui gen `_t' = _n

		qui gen `_zero' = 0
		
		qui gen `birf'    = 0 if _n<=`h1'
		qui gen `seirf'   = 0 if _n<=`h1'
		qui gen `birf_up' = 0 if _n<=`h1'
		qui gen `birf_lo' = 0 if _n<=`h1'
		
		* estimation
		forval h=`hran' {
            if `hs'<=0 loc k = `h' + 1 - `hs'
            else       loc k = `h'

			loc hstr = `h' - `hs'
			
			local depvar_transf ``y'`hstr''
			tempvar touse_current
            qui gen byte `touse_current' = `touse'
            qui markout `touse_current' `depvar_transf' `indepvars_list' `idvar' `timevar'
            
			* Drop singletons on this horizon's full sample (depvar_transf marked
			* out), matching reghdfe. Each horizon's lead differs so its sample
			* and singletons can differ; report per-horizon with the h tag.
			if `_sing_drop' {
				_xtlp_drop_singletons, touse(`touse_current') idvar(`idvar') timevar(`timevar') fetype(`fe_type')
				local _hdropped = r(dropped)
				if `_hdropped' > 0 {
					di as txt `"(h=`h': dropped `_hdropped' singleton observations)"'
				}
			}

			mata: st_numscalar("XTLP_H", `h')
			mata: lp_work("`depvar_transf'", "`indepvars_list'", "`touse_current'", "`idvar'", "`timevar'", `fe_type', `method_code', `vce_code', `dk_lag', `ase_code', `nodfadj_code', "`clustervar'", "`clustervar2'", 0)
			
			matrix colnames b = `indepvars_list'
			matrix colnames V = `indepvars_list'
			matrix rownames V = `indepvars_list'
			
			matrix b`h' = b
			matrix V`h' = V
			ereturn matrix b`h' = b`h'
			ereturn matrix V`h' = V`h'
			ereturn scalar df_r`h' = scalar(df_r)
			if `vce_code' == 4 {
				ereturn scalar lag`h'  = scalar(lag)
			}

			* only write to the IRF plot variables when row k exists (k<=h1);
			* the estimate is already saved to e() above regardless.
			if `k' <= `h1' {
				local coef = b[1, 1]
				local var  = V[1, 1]
				local se   = sqrt(`var')
	            local ub   = `coef' + 1.96 * `se'
	            local lb   = `coef' - 1.96 * `se'
				quietly {
					replace `birf'    = `coef' if _n == `k'
	                replace `seirf'   = `se'   if _n == `k'
	                replace `birf_up' = `ub'   if _n == `k'
	                replace `birf_lo' = `lb'   if _n == `k'
				}
			}
		}
		
		**### header
		if `method_code' == 1 {
			di _n as txt "Method: " as result "Fixed effects"
		}
		else {
			di _n as txt "Method: " as result "Split-panel jackknife"
		}

		if `vce_code' == 0 {
			di as txt "Standard errors: " as result "Unadjusted"
		}
		else if `vce_code' == 1 {
			di as txt "Standard errors: " as result "Robust"
		}
		else if `vce_code' == 2 {
			di as txt "Standard errors: " as result "One-way clustered"
		}
		else if `vce_code' == 3 {
			di as txt "Standard errors: " as result "Two-way clustered"
		}
		else if `vce_code' == 4 {
			if `nodfadj_code' == 1 {
				di as txt "Standard errors: " as result "Driscoll-Kraay (nodfadj)"
			}
			else {
				di as txt "Standard errors: " as result "Driscoll-Kraay"
			}
		}

		if `ase_code' == 0 {
			di as txt "Small-sample adjustment: " as result "On"
		}
		else if `ase_code' == 1 {
			di as txt "Small-sample adjustment: " as result "Off"
		}

		**### ereturn
		mkmat `birf'    if _n<=`h1', mat(BIRF)
		mkmat `seirf'   if _n<=`h1', mat(SEIRF)
		mkmat `birf_lo' if _n<=`h1', mat(SEIRF_LO)
		mkmat `birf_up' if _n<=`h1', mat(SEIRF_UP)
		
		mat IRF = BIRF , SEIRF , SEIRF_LO , SEIRF_UP
		matrix colnames IRF = "IRF" "Std. err." "95% CI Lower" "95% CI Upper"

		* row labels: exactly h1 labels (hs, hs+1, ...) matching the h1 stored
		* rows (may be fewer than the requested range if capped at _N above).
		loc rows ""
		loc rstart = `hs'
		forval i=1/`h1' {
			loc rows `rows' `=`rstart' + `i' - 1'
		}
		matrix rownames IRF = `rows'
		matlist IRF, noheader format(%12.5f) title("Impulse Response Function") lines(oneline) rowtitle("Horizon")
		ereturn local vce     "`vce_string'"
		ereturn local vcetype "`vcetype_string'"
		ereturn local cmd      "xtlp"
		ereturn local cmdline  `"xtlp `cmdline'"'
		ereturn matrix irf = IRF
		
		**### graph
        loc mod = mod(`hor' - `hs', 2)
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
            qui twoway ///
                (rarea `birf_up' `birf_lo' `_t', fcolor(`lcolor'%15) lc(`lcolor'%7)) ///
				(line `_zero' `_t', lcolor(gs5) lpattern(dash)) ///
                (line `birf' `_t', lcolor(`lcolor') lpattern(solid) lwidth(medthick)) ///
                if _n<=`h1', ///
                graphregion(color(white)) ///
				plotregion(margin(zero)) ///
				legend(`off' order(3 "IRF of `y' (`depvar') to shock (`shockvar'), method(`method')") position(6)) ///
                tlabel(`hs'(`p')`hor') ///
                xtitle("Horizon") ///
				name("IRF_`method'", replace)
		}
	}
	else {
		**## h>0 shock=2
		
		* plot data
		cap drop _birf _seirf _birf_lo _birf_up
		tempvar _t _zero
		
        if `hs'<=0 loc h1 = `hor' + 1 - `hs'
        else       loc h1 = `hor'

        * Cap the stored IRF rows at the observation count (see single-shock
        * branch for rationale): when the requested horizon range exceeds _N
        * there are no rows to hold the extra horizons; cap to avoid an
        * mkmat/rownames conformability error.
        local nobs = _N
        if `h1' > `nobs' {
            di as txt "note: requested horizon range (`h1' rows) exceeds the number of observations (" `nobs' "); IRF table truncated to " `nobs' " rows."
            loc h1 = `nobs'
        }

        if `hs'<=0 qui gen `_t' = _n - 1 + `hs'
        else       qui gen `_t' = _n
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
            if `hs'<=0 loc k = `h' + 1 - `hs'
            else       loc k = `h'

			loc hstr = `h' - `hs'
			
			local depvar_transf ``y'`hstr''
			tempvar touse_current
            qui gen byte `touse_current' = `touse'
            qui markout `touse_current' `depvar_transf' `indepvars_list' `idvar' `timevar'
            
			* Drop singletons on this horizon's full sample (depvar_transf marked
			* out), matching reghdfe. Each horizon's lead differs so its sample
			* and singletons can differ; report per-horizon with the h tag.
			if `_sing_drop' {
				_xtlp_drop_singletons, touse(`touse_current') idvar(`idvar') timevar(`timevar') fetype(`fe_type')
				local _hdropped = r(dropped)
				if `_hdropped' > 0 {
					di as txt `"(h=`h': dropped `_hdropped' singleton observations)"'
				}
			}

			mata: lp_work("`depvar_transf'", "`indepvars_list'", "`touse_current'", "`idvar'", "`timevar'", `fe_type', `method_code', `vce_code', `dk_lag', `ase_code', `nodfadj_code', "`clustervar'", "`clustervar2'", 0)
			
			matrix colnames b = `indepvars_list'
			matrix colnames V = `indepvars_list'
			matrix rownames V = `indepvars_list'
			
			matrix b`h' = b
			matrix V`h' = V
			ereturn matrix b`h' = b`h'
			ereturn matrix V`h' = V`h'
			ereturn scalar df_r`h' = scalar(df_r)
			if `vce_code' == 4 {
				ereturn scalar lag`h'  = scalar(lag)
			}

			* row label only when this horizon has a storable row (k<=h1); the
			* estimate is saved to e() above regardless.
			if `k' <= `h1' {
				local row_names "`row_names' `h'"
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
		}
		
		**### header
		if `method_code' == 1 {
			di _n as txt "Method: " as result "Fixed effects"
		}
		else {
			di _n as txt "Method: " as result "Split-panel jackknife"
		}

		if `vce_code' == 0 {
			di as txt "Standard errors: " as result "Unadjusted"
		}
		else if `vce_code' == 1 {
			di as txt "Standard errors: " as result "Robust"
		}
		else if `vce_code' == 2 {
			di as txt "Standard errors: " as result "One-way clustered"
		}
		else if `vce_code' == 3 {
			di as txt "Standard errors: " as result "Two-way clustered"
		}
		else if `vce_code' == 4 {
			if `nodfadj_code' == 1 {
				di as txt "Standard errors: " as result "Driscoll-Kraay (nodfadj)"
			}
			else {
				di as txt "Standard errors: " as result "Driscoll-Kraay"
			}
		}

		if `ase_code' == 0 {
			di as txt "Small-sample adjustment: " as result "On"
		}
		else if `ase_code' == 1 {
			di as txt "Small-sample adjustment: " as result "Off"
		}

		**### ereturn
		local sidx = 0
		local first_matrix = 1
		foreach idx of local shocklist {
			local ++sidx
			local shockvar : word `idx' of `indepvars_list'
			
			tempname BIRF SEIRF SEIRF_LO SEIRF_UP IRF_s
			
			mkmat `birf`sidx''    if _n <= `h1', mat(`BIRF')
			mkmat `seirf`sidx''   if _n <= `h1', mat(`SEIRF')
			mkmat `birf_lo`sidx'' if _n <= `h1', mat(`SEIRF_LO')
			mkmat `birf_up`sidx'' if _n <= `h1', mat(`SEIRF_UP')
			
			matrix `IRF_s' = `BIRF', `SEIRF', `SEIRF_LO', `SEIRF_UP'
			local cname1 "IRF"
			local cname2 "Std. err."
			local cname3 "95% CI Lower"
			local cname4 "95% CI Upper"
			matrix colnames `IRF_s' = "`cname1'" "`cname2'" "`cname3'" "`cname4'"
			matrix rownames `IRF_s' = `row_names'
			
			if `first_matrix' == 1 {
                matrix `IRF_all' = `IRF_s'
                local first_matrix = 0
            }
            else {
                matrix `IRF_all' = `IRF_all', `IRF_s'
            }
	
			di _n as txt "Impulse Response Function for shock #" `sidx' ///
				  as txt " (" as result "`shockvar'" as txt ")"
			
			matlist `IRF_s', noheader format(%12.5f) ///
				lines(oneline) rowtitle("Horizon")
			
			**### graph
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
					(line `birf`sidx'' `_t', lcolor(`lcolor') lpattern(solid) lwidth(medthick)) ///
					if _n <= `h1', ///
					graphregion(color(white)) ///
					plotregion(margin(zero)) ///
					legend(`off' order(3 "IRF of `y' (`depvar') to shock (`shockvar'), method(`method')") position(6)) ///
					tlabel(`hs'(`p')`hor') ///
					xtitle("Horizon") ///
					name("IRF_`method'_`sidx'", replace)
			}
		}
		matrix rownames `IRF_all' = `row_names'
		ereturn local vce     "`vce_string'"
		ereturn local vcetype "`vcetype_string'"
		ereturn local cmd      "xtlp"
		ereturn local cmdline  `"xtlp `cmdline'"'
		ereturn matrix irf = `IRF_all'
	}
end

* -----------------------------------------------------------------------------
* _xtlp_drop_singletons
* Drop singleton groups from a 0/1 sample indicator, replicating reghdfe's
* iterative single-FE drop. Singletons are obs whose absorbed-FE group has
* exactly 1 obs in the current sample:
*   - fetype==1 (fe) : drop id-singletons only (time FE not absorbed)
*   - fetype==2 (tfe): drop id-singletons AND time-singletons, iteratively
*     (dropping one set can create new singletons in the other).
* Group size is counted with -egen total()- (broadcast to EVERY row of the
* group, including touse==0 rows) -- NOT gen sum() if touse, which misses
* singletons whose touse row is not last in the group.
* Modifies the sample var in place; returns the number dropped in r(dropped).
* -----------------------------------------------------------------------------
program define _xtlp_drop_singletons, rclass sortpreserve
	syntax , TOUSE(varname) IDvar(varname) TIMEvar(varname) FETYPE(integer)

	tempvar _gsize _tsize
	local _iter 0
	local _dropped 1
	local _total_dropped 0
	while `_dropped' > 0 {
		local _iter = `_iter' + 1
		local _dropped 0
		* id-singletons
		qui bysort `idvar': egen long `_gsize' = total(`touse')
		qui count if `touse' & `_gsize' == 1
		local _dropped = r(N)
		qui replace `touse' = 0 if `touse' & `_gsize' == 1
		drop `_gsize'
		* time-singletons (only under tfe; under single-fe the time FE is
		* not absorbed so time-singletons are not degenerate)
		if `fetype' == 2 {
			qui bysort `timevar': egen long `_tsize' = total(`touse')
			qui count if `touse' & `_tsize' == 1
			local _dropped = `_dropped' + r(N)
			qui replace `touse' = 0 if `touse' & `_tsize' == 1
			drop `_tsize'
		}
		local _total_dropped = `_total_dropped' + `_dropped'
	}
	return scalar dropped = `_total_dropped'
end

* -----------------------------------------------------------------------------
* MATA CODE BLOCK
* -----------------------------------------------------------------------------
version 14.0
set matalnum on
mata:

// Is the categorical C constant within each group defined by `info`
// (panelsetup on some FE variable)? If yes, that FE is nested within C
// (every FE level maps to one C value), so the FE does not reduce the
// cluster-robust DoF (aligned to reghdfe's dof_update_nested).
// NOTE: `info` must come from panelsetup() on a variable whose levels occupy
// CONTIGUOUS rows in the data (true for ID after the (ID,Time) sort, but NOT
// for Time, which repeats 1..T per individual). For Time use nested_in_byval().
real scalar nested_in(real colvector C, real matrix info)
{
	real scalar i, s, e
	for (i = 1; i <= rows(info); i++) {
		s = info[i,1]
		e = info[i,2]
		if (max(C[|s \ e|]) - min(C[|s \ e|]) != 0) return(0)
	}
	return(1)
}

// Sort-free version of nested_in: group rows by the unique values of `key`
// (which need NOT be contiguous) and test whether C is constant within each
// key group. Used for the time FE, because Time is not contiguous after the
// (ID,Time) sort. Returns 1 if the FE defined by `key` is nested within C.
real scalar nested_in_byval(real colvector C, real colvector key)
{
	real colvector ukey
	real scalar i
	real colvector idx
	ukey = uniqrows(key)
	for (i = 1; i <= rows(ukey); i++) {
		idx = selectindex(key :== ukey[i])
		if (max(C[idx]) - min(C[idx]) != 0) return(0)
	}
	return(1)
}

// Check the demeaned-regressor cross product G = Xd'Xd for collinearity /
// singularity. cholinv() returns a matrix of missings when G is not positive
// definite (perfectly collinear regressors, or a degenerate subsample such as
// the SPJ A/B halves when too few periods per individual, or a multi-horizon
// run where the leads have exhausted the data at this horizon).
// Returns 1 if X'X is non-singular (ok to proceed); 0 if singular.
//   fatal=1 (single-result h=0 call): on singular, give a readable error and
//     exit(198) instead of letting missing SEs propagate to rc=504 later.
//   fatal=0 (per-horizon h>0 call): on singular, post missing b/V and return
//     0 so lp_work returns early -- the IRF row shows "." for this horizon,
//     matching the established out-of-range-horizon behavior (do NOT abort the
//     whole multi-horizon run).
// `where` labels the offending sample (full / SPJ part A / SPJ part B).
real scalar check_nonsingular(real matrix XX_inv, real matrix X_dot, string scalar where, real scalar fatal)
{
	if (!hasmissing(XX_inv)) return(1)

	if (fatal) {
		errprintf("xtlp: regressors are collinear or the sample is degenerate")
		errprintf(" (" + where + "); cannot invert X'X.\n")
		if (strpos(where, "SPJ") > 0) {
			errprintf("         method(spj) needs at least 4 time periods per")
			errprintf(" individual (T>=4) so each split half has >=2 obs.\n")
		}
		else {
			errprintf("         drop the redundant regressor, or use fewer")
			errprintf(" covariates / a larger sample.\n")
		}
		exit(198)
	}
	else {
		// per-horizon: too little data at this horizon (e.g. 1 obs/individual
		// after FE demeaning). Post missing results and signal early return.
		printf("{txt}note: too few observations to estimate at this horizon")
		printf(" (" + where + "); SE set to missing.\n")
		st_matrix("b", J(1, cols(X_dot), .))
		st_matrix("V", J(cols(X_dot), cols(X_dot), .))
		st_numscalar("N",    0)
		st_numscalar("N_g",  0)
		st_numscalar("df_r", .)
		st_numscalar("lag",  .)
		return(0)
	}
}

void lp_work(string scalar depvar, ///
			 string scalar indepvars, ///
			 string scalar touse, ///
             string scalar idvar, ///
			 string scalar timevar, ///
             real   scalar fe_type, ///
             real   scalar method_code, ///
             real   scalar vce_code, ///
             real   scalar dk_lag, ///
             real   scalar ase_code, ///
             real   scalar nodfadj_code, ///
             string scalar clustervar, ///
             string scalar clustervar2, ///
             real   scalar fatal_collinear)
{
    real matrix Y, X, info, V_est, b_est
    real colvector ID, Time, C1, C2
    real scalar N, K, N_g
    string rowvector xnames
    
    // Load data
    // C1 = first cluster variable; C2 = second (twoway). Either may be any
    // variable (cluster dims are independent of absorbed FE, aligned to reghdfe).
    // Reuse ID/Time when the clustervar equals them to avoid a redundant load.
    Y    = st_data(., depvar, touse)
    X    = st_data(., indepvars, touse)
    ID   = st_data(., idvar, touse)
    Time = st_data(., timevar, touse)
    C1   = ((vce_code == 2 | vce_code == 3) & clustervar != idvar)   ? st_data(., clustervar, touse)   : ID
    C2   = (vce_code == 3 & clustervar2 != timevar)                  ? st_data(., clustervar2, touse)  : Time

    // Get variable names for labeling
    xnames = tokens(indepvars)
    
    // Panel info setup (sort by ID, Time)
    // We need to sort X, Y, ID, Time together to ensure structure
    real matrix ALL
    real scalar needC1, needC2
    needC1 = ((vce_code == 2 | vce_code == 3) & clustervar != idvar)
    needC2 = (vce_code == 3 & clustervar2 != timevar)
    if (needC1 & needC2) {
        ALL = Y, ID, Time, C1, C2, X
        _sort(ALL, (2,3)) // Sort by ID then Time
        Y    = ALL[., 1]
        ID   = ALL[., 2]
        Time = ALL[., 3]
        C1   = ALL[., 4]
        C2   = ALL[., 5]
        X    = ALL[., 6::cols(ALL)]
    }
    else if (needC1) {
        ALL = Y, ID, Time, C1, X
        _sort(ALL, (2,3)) // Sort by ID then Time
        Y    = ALL[., 1]
        ID   = ALL[., 2]
        Time = ALL[., 3]
        C1   = ALL[., 4]
        X    = ALL[., 5::cols(ALL)]
    }
    else if (needC2) {
        ALL = Y, ID, Time, C2, X
        _sort(ALL, (2,3)) // Sort by ID then Time
        Y    = ALL[., 1]
        ID   = ALL[., 2]
        Time = ALL[., 3]
        C2   = ALL[., 4]
        X    = ALL[., 5::cols(ALL)]
    }
    else {
        ALL = Y, ID, Time, X
        _sort(ALL, (2,3)) // Sort by ID then Time
        Y    = ALL[., 1]
        ID   = ALL[., 2]
        Time = ALL[., 3]
        X    = ALL[., 4::cols(ALL)]
    }

    N = rows(Y)
    K = cols(X)

    // Degenerate/empty sample (e.g. a multi-horizon run where the lead f`h'.y
    // is missing for every observation at this horizon, or an `if` that drops
    // all rows). Post missing results and return silently so the IRF row shows
    // "." rather than erroring mid-loop (preserves the out-of-range-horizon
    // behavior); genuine collinearity with a non-empty sample is caught later.
    if (N == 0) {
        st_matrix("b", J(1, K, .))
        st_matrix("V", J(K, K, .))
        st_numscalar("N",    0)
        st_numscalar("N_g",  0)
        st_numscalar("df_r", .)
        st_numscalar("lag",  .)
        return
    }

    // Panel setup info: [start_index, end_index] for each individual
    info = panelsetup(ID, 1)
    N_g  = rows(info) // Number of groups (individuals)

    // -------------------------------------------------------
	//# Step 1: Define Split Points  (SPJ only)
    // -------------------------------------------------------
    real colvector T_a_idx, T_b_idx, cut_i
    real colvector selA, selB
    real scalar i, start, end_t, Ti, cut_idx

    if (method_code == 2) {
        T_a_idx = J(N, 1, 0)
        T_b_idx = J(N, 1, 0)
        cut_i   = J(N_g, 1, .)

        for (i=1; i<=N_g; i++) {
            start = info[i, 1]
            end_t = info[i, 2]
            Ti    = end_t - start + 1

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

        selA = selectindex(T_a_idx :!= 0)
        selB = selectindex(T_b_idx :!= 0)
    }

    // -------------------------------------------------------
    //# Step 2: Within Transformation (Demeaning)
    // -------------------------------------------------------
	// Full sample  (always needed: FE uses it directly; SPJ uses it as X_dot)
    real matrix    YX, YX_dm
    real colvector Y_dot
    real matrix    X_dot
    YX    = Y, X
    YX_dm = twoway_demean(YX, ID, Time, info, fe_type)
    Y_dot = YX_dm[., 1]
    X_dot = YX_dm[., 2..(K+1)]

    // Split-sample demeaning (SPJ only)
    real colvector Y_A, ID_A, Time_A, Y_B, ID_B, Time_B, Y_a_dot, Y_b_dot
    real matrix    X_A, X_B, infoA, infoB, YX_a, YX_a_dm, X_a_dot, YX_b, YX_b_dm, X_b_dot
    real matrix    X_dot_a_full, X_dot_b_full
    real scalar NT
    NT = rows(X_dot)

    if (method_code == 2) {
        // A sample
        Y_A    = Y[selA, .]
        X_A    = X[selA, .]
        ID_A   = ID[selA, .]
        Time_A = Time[selA, .]
        infoA  = panelsetup(ID_A, 1)

        YX_a    = Y_A, X_A
        YX_a_dm = twoway_demean(YX_a, ID_A, Time_A, infoA, fe_type)
        Y_a_dot = YX_a_dm[., 1]
        X_a_dot = YX_a_dm[., 2..cols(YX_a_dm)]

        // B sample
        Y_B    = Y[selB, .]
        X_B    = X[selB, .]
        ID_B   = ID[selB, .]
        Time_B = Time[selB, .]
        infoB  = panelsetup(ID_B, 1)

        YX_b    = Y_B, X_B
        YX_b_dm = twoway_demean(YX_b, ID_B, Time_B, infoB, fe_type)
        Y_b_dot = YX_b_dm[., 1]
        X_b_dot = YX_b_dm[., 2..cols(YX_b_dm)]

        // Prepare full-length split demeaned X
        X_dot_a_full = J(NT, K, .)
        X_dot_b_full = J(NT, K, .)
        X_dot_a_full[selA, .] = X_a_dot
        X_dot_b_full[selB, .] = X_b_dot
    }

	// -------------------------------------------------------
    //# Step 3: Estimate Coefficients (OLS on demeaned data)
    // -------------------------------------------------------
    real colvector b_full, b_a, b_b
    real matrix XX_inv_full, XX_inv_a, XX_inv_b
  
    // Full sample  (always needed)
    XX_inv_full = cholinv(cross(X_dot, X_dot))
    if (!check_nonsingular(XX_inv_full, X_dot, "full sample", fatal_collinear)) return
    b_full      = XX_inv_full * cross(X_dot, Y_dot)

    // Part A / B  (SPJ only)
    if (method_code == 2) {
        XX_inv_a = cholinv(cross(X_a_dot, X_a_dot))
        if (!check_nonsingular(XX_inv_a, X_a_dot, "SPJ part A", fatal_collinear)) return
        b_a      = XX_inv_a * cross(X_a_dot, Y_a_dot)

        XX_inv_b = cholinv(cross(X_b_dot, X_b_dot))
        if (!check_nonsingular(XX_inv_b, X_b_dot, "SPJ part B", fatal_collinear)) return
        b_b      = XX_inv_b * cross(X_b_dot, Y_b_dot)
    }
  
    // -------------------------------------------------------
    //# Step 4: Estimator
    // -------------------------------------------------------
	if (method_code == 1) {
	    b_est = b_full
	}
	else if (method_code == 2) {
	    b_est = 2 * b_full - 0.5 * (b_a + b_b)
	}
  
    // -------------------------------------------------------
    //# Step 5: Variance Calculation
    // -------------------------------------------------------
	real matrix X_mat, X_a_mat, X_b_mat, X_sub, d_dot

    X_mat   = X_dot
    d_dot   = X_dot          // FE: d_dot = within-demeaned X (no split needed)

    if (method_code == 2) {
        // SPJ: d_dot = 2*X_dot - X_sub, where X_sub stitches A/B demeaned halves
        X_a_mat = X_dot_a_full
        X_b_mat = X_dot_b_full
        X_sub   = J(NT, K, .)

        for (i=1; i<=N_g; i++) {
            real scalar start1, end_t1, Ti1, ci

            start1 = info[i,1]
            end_t1 = info[i,2]
            Ti1    = end_t1 - start1 + 1
            ci     = cut_i[i]

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

        d_dot = 2:*X_mat :- X_sub
    }

	
	real colvector e
    e = Y_dot - X_dot * b_est      // NT x 1

	// W_N = first-cluster-dim meat (sum over clustervar groups of g g').
	// Grouped by C1 (the first cluster variable, any variable; aligned to reghdfe:
	// cluster dimensions are independent of the absorbed FE). Used by cluster
	// (vce_code==2) as M_c1 and by twoway (vce_code==3) as M_c1. Conventional/
	// robust/dkraay do not use it -> skip.
	real matrix W_N
	real colvector c1uniq
	real scalar G_c1
	if (vce_code == 2 | vce_code == 3) {
		c1uniq = uniqrows(sort(C1, 1))
		G_c1   = rows(c1uniq)
		W_N    = J(K, K, 0)

		real matrix d_i
		real colvector e_i, g_i, idx_c1
		real scalar cc1
		for (cc1 = 1; cc1 <= G_c1; cc1++) {
			idx_c1 = selectindex(C1 :== c1uniq[cc1])

			d_i = d_dot[idx_c1, .]                  // n_c1 x K
			e_i = e[idx_c1]                         // n_c1 x 1
			g_i = d_i' * e_i                        // K x 1

			W_N = W_N + g_i * g_i'                  // K x K
		}
	}

	// Sandwich: var0_mat = XX_inv * W * XX_inv
	real matrix XX_inv
	XX_inv = XX_inv_full            // = (X_dot'X_dot)^-1, computed in Step 3

	real scalar adj, rank
	rank = K + 1   // within absorbs FE, no constant kept => rank = K+1 (matches xtscc e(rank))
	real scalar df_r_val

	// df_r for conventional/robust: subtract absorbed FE dof and slopes (K).
	// This DOES deduct df_a, unlike cluster/dkraay/twoway which align to xtscc.
	//   fe  (fe_type==1): df_r = NT - N_g - K               (absorb individual FE only)
	//   tfe (fe_type==2): df_r = NT - N_g - G_t + 1 - K     (absorb individual + time FE;
	//                     +1 because the two-way absorption has one redundant global
	//                     constant, aligned to reghdfe absorb(id t) df_r).
	// Matches xtreg,fe (single FE) and reghdfe unadjusted/robust df_r.
	real scalar df_r_conv, rss, G_t_conv
	if (vce_code == 0 | vce_code == 1) {
		if (fe_type == 2) {
			G_t_conv   = rows(uniqrows(Time))
			df_r_conv  = NT - N_g - G_t_conv + 1 - K
		}
		else {
			df_r_conv  = NT - N_g - K
		}
		rss = quadcross(e, e)          // within residual sum of squares

		// residual dof must be positive, else the small-sample factor NT/df_r
		// divides by <=0 and SEs become missing (rc=504 later).
		if (df_r_conv <= 0) {
			if (fatal_collinear) {
				errprintf("xtlp: too few degrees of freedom for this vce.\n")
				errprintf("         N*T=" + strofreal(NT) + " - N_g=" + strofreal(N_g)
				          + " - K=" + strofreal(K) + " = " + strofreal(df_r_conv)
				          + " <= 0.\n")
				errprintf("         Reduce the number of regressors or use more")
				errprintf(" observations.\n")
				exit(198)
			}

			// per-horizon: too little data -> missing SE for this horizon
			st_matrix("b", J(1, K, .))
			st_matrix("V", J(K, K, .))
			st_numscalar("N", N); st_numscalar("N_g", N_g)
			st_numscalar("df_r", .); st_numscalar("lag", .)
			return
		}
	}

	if (vce_code == 0) {
		// ============================================================
		//## Conventional (homoskedastic) SE, aligned to xtreg,fe / reghdfe unadjusted
		//  M = sum_i d_i d_i' * e_i^2; V = XX_inv * M * XX_inv * (NT/df_r)
		// Homoskedastic replacement: e_i^2 -> SSR / NT
		// ase=1 => NT/(NT-1) instead of NT/df_r
		// ============================================================
		real matrix M_homo
		M_homo = (rss / NT) * cross(d_dot, d_dot)	
		if (ase_code == 1)  adj = NT / (NT - 1)
		else                adj = NT / df_r_conv
		V_est = XX_inv * (M_homo * adj) * XX_inv
		df_r_val = df_r_conv
	}
	else if (vce_code == 1) {
		// ============================================================
		//## Robust HC1 SE (heteroskedasticity-only, no clustering),
		// aligned to reghdfe vce(robust)
		//   M = sum_i d_i d_i' * e_i^2  ;  V = XX_inv * M * XX_inv * (NT/df_r)
		// d_dot = within-demeaned regressors (FE: X_dot; SPJ: 2*X_dot - X_sub)
		// ase=1 => NT/(NT-1) instead of NT/df_r
		// ============================================================
		real matrix M_hc1
		M_hc1 = quadcross(d_dot, (e :^ 2), d_dot)   // K x K, weighted by e_i^2
		if (ase_code == 1)  adj = NT / (NT - 1)
		else                adj = NT / df_r_conv
		V_est = XX_inv * (M_hc1 * adj) * XX_inv
		df_r_val = df_r_conv
	}
	else if (vce_code == 2) {
		// ============================================================
		//## One-way cluster SE by clustervar C1 (Arellano 1987), aligned
		// to reghdfe. Cluster dims are independent of the absorbed FE.
		//   dof_adj = (N-1)/(N - nested_adj - K - df_a_eff) * G_c1/(G_c1-1)
		//   df_a_eff = absorbed FE dof surviving nesting + redundant global-constant
		//   accounting (aligned to reghdfe dof_update_nested):
		//     fe : R_id = id_nested ? N_g : 0
		//          df_a_eff = N_g - R_id ;  nested_adj = id_nested ? 1 : 0
		//     tfe: R_t  = t_nested  ? G_t : 1   (time FE is the 2nd intercept)
		//          R_id = id_nested ? N_g : (t_nested ? 1 : 0)
		//          df_a_eff = (N_g - R_id) + (G_t - R_t)
		//          nested_adj = (id_nested | t_nested) ? 1 : 0
		//   ase => dof_adj = G_c1/(G_c1-1)
		// ============================================================
		real scalar id_nested_c1, t_nested_c1, G_t_c1, df_a_eff_c1, nested_adj_c1
		real scalar R_id_c1, R_t_c1
		// Need >=2 clusters, else G_c1/(G_c1-1) divides by zero -> missing SE.
		if (G_c1 <= 1) {
			if (fatal_collinear) {
				errprintf("xtlp: vce(cluster " + clustervar + ") has only "
				          + strofreal(G_c1) + " cluster group")
				if (G_c1 == 1) errprintf("; need at least 2.\n")
				else           errprintf("s; need at least 2.\n")
				exit(198)
			}
			st_matrix("b", J(1, K, .)); st_matrix("V", J(K, K, .))
			st_numscalar("N", N); st_numscalar("N_g", N_g)
			st_numscalar("df_r", .); st_numscalar("lag", .)
			return
		}

		id_nested_c1 = nested_in(C1, info)
		if (fe_type == 2) {
			// Time is NOT contiguous after the (ID,Time) sort (it repeats 1..T per
			// individual), so panelsetup(Time) is wrong here. Count unique time
			// values and test nesting sort-free.
			G_t_c1 = rows(uniqrows(Time))
			t_nested_c1 = nested_in_byval(C1, Time)
			// R_t: the absorbed time FE always gives up >=1 global-constant dof
			//      (it is the 2nd intercept); all of them if nested in the cluster.
			// R_id: all id dof if id nested; else if t is nested (so some FE is
			//      nested) id also surrenders 1 global-constant dof; else 0.
			R_t_c1  = t_nested_c1  ? G_t_c1 : 1
			R_id_c1 = id_nested_c1 ? N_g    : (t_nested_c1 ? 1 : 0)
			df_a_eff_c1  = (N_g - R_id_c1) + (G_t_c1 - R_t_c1)
			nested_adj_c1 = (id_nested_c1 | t_nested_c1) ? 1 : 0
		}
		else {
			df_a_eff_c1 = id_nested_c1 ? 0 : N_g
			nested_adj_c1 = id_nested_c1 ? 1 : 0
		}

		if (ase_code == 1) {
			adj = G_c1 / (G_c1 - 1.0)
		}
		else {
			// small-sample denominator must stay positive
			real scalar c1_denom
			c1_denom = NT - nested_adj_c1 - K - df_a_eff_c1
			if (c1_denom <= 0) {
				if (fatal_collinear) {
					errprintf("xtlp: vce(cluster " + clustervar + ") has too few degrees"
					          + " of freedom for the small-sample adjustment (denominator"
					          + " = " + strofreal(c1_denom) + " <= 0).\n")
					errprintf("         Use vce(cluster " + clustervar + ", ase) to drop"
					          + " the small-sample adjustment, or reduce regressors / add"
					          + " data.\n")
					exit(198)
				}
				st_matrix("b", J(1, K, .)); st_matrix("V", J(K, K, .))
				st_numscalar("N", N); st_numscalar("N_g", N_g)
				st_numscalar("df_r", .); st_numscalar("lag", .)
				return
			}
			adj = ((NT - 1.0) / c1_denom) * (G_c1 / (G_c1 - 1.0))
		}
		V_est = XX_inv * (W_N * adj) * XX_inv
		df_r_val = G_c1 - 1
		st_numscalar("N_clust_tw", G_c1)
	}
	else if (vce_code == 3) {
		// ============================================================
		//## Two-way cluster SE (Cameron-Gelbach-Miller 2011), c1 x c2
		//   M = M_c1 + M_c2 - M_(c1#c2)   (inclusion-exclusion)
		//   V = XX_inv * M * XX_inv * dof_adj
		//   dof_adj = (N-1)/(N - nested_adj - K - df_a_eff) * N_clust/(N_clust-1)
		//   ase => dof_adj = N_clust/(N_clust-1)
		//   PSD fix if V not positive semi-definite (reghdfe_fix_psd)
		// ============================================================
		real matrix M_time, M_idt, M_tw, g_idt, S_idt
		real colvector g_time, c2uniq, c1c2_key, ukey
		real scalar G_time, N_clust, tt2, g_idc2
		real matrix U_tw
		real colvector lam_tw

		// M_c1 = W_N (already summed over C1 groups above). Reuse.

		// M_time: group by the second cluster variable C2 (set in lp_work's
		// data-load section).
		c2uniq = uniqrows(sort(C2, 1))
		G_time = rows(c2uniq)
		M_time = J(K, K, 0)
		for (tt2 = 1; tt2 <= G_time; tt2++) {
			real colvector idx_c2
			idx_c2 = selectindex(C2 :== c2uniq[tt2])
			g_time = d_dot[idx_c2, .]' * e[idx_c2]   // K x 1
			M_time = M_time + g_time * g_time'
		}

		// M_(c1#c2): intersection of the two cluster dimensions. Group
		// observations by the (C1, C2) pair and sum the per-obs scores:
		//   M_idt = sum_g (sum_{i in g} d_i e_i)(sum_{i in g} d_i e_i)'
		// When C1=id & C2=time in a balanced panel, each (id,time) pair
		// has exactly one obs, so this collapses to HC0. For general
		// clustervars the grouping is essential.
		// composite key = C1 scaled up by (max C2 + 1) + C2; unique pairs
		// map to unique keys.
		c1c2_key = C1 :* (max(C2) + 1) :+ C2
		ukey     = uniqrows(c1c2_key)
		M_idt = J(K, K, 0)
		g_idt = d_dot :* e            // NT x K, per-obs score
		for (g_idc2 = 1; g_idc2 <= rows(ukey); g_idc2++) {
			idx_c2 = selectindex(c1c2_key :== ukey[g_idc2])
			S_idt = colsum(g_idt[idx_c2, .])'   // K x 1 group score
			M_idt = M_idt + S_idt * S_idt'
		}

		M_tw = W_N + M_time - M_idt

		N_clust = min((G_c1, G_time))
		// Need >=2 clusters in the smaller dimension, else N_clust/(N_clust-1)
		// divides by zero -> missing SE.
		if (N_clust <= 1) {
			if (fatal_collinear) {
				errprintf("xtlp: vce(cluster " + clustervar + " " + clustervar2
				          + ") has only " + strofreal(N_clust) + " cluster group")
				if (N_clust == 1) errprintf(" in the smaller dimension; need at least 2.\n")
				else           errprintf("s in the smaller dimension; need at least 2.\n")
				exit(198)
			}
			st_matrix("b", J(1, K, .)); st_matrix("V", J(K, K, .))
			st_numscalar("N", N); st_numscalar("N_g", N_g)
			st_numscalar("df_r", .); st_numscalar("lag", .)
			return
		}
		// Nested DoF (aligned to reghdfe): an absorbed FE nested within
		// EITHER cluster dim does not reduce the cluster-robust DoF.
		real scalar id_nested_tw, t_nested_tw, df_a_eff_tw, nested_adj_tw
		real scalar R_id_tw, R_t_tw, G_t_tw
		id_nested_tw = nested_in(C1, info) | nested_in(C2, info)
		if (fe_type == 2) {
			// Time is not contiguous after the (ID,Time) sort, so test nesting
			// of the absorbed time FE (Time) sort-free.
			// NOTE: G_t_tw is the # of absorbed TIME-FE levels (unique Time
			// values) -- NOT G_time, which is the # of levels of the SECOND
			// clustervar (C2, possibly != timevar, e.g. region).
			G_t_tw = rows(uniqrows(Time))
			t_nested_tw = nested_in_byval(C1, Time) | nested_in_byval(C2, Time)
			// Same nesting/redundant-constant accounting as the one-way case:
			//   R_t  = G_t_tw if t nested in either cluster dim, else 1
			//   R_id = N_g    if id nested in either dim, else (t_nested ? 1 : 0)
			R_t_tw  = t_nested_tw  ? G_t_tw : 1
			R_id_tw = id_nested_tw ? N_g    : (t_nested_tw ? 1 : 0)
			df_a_eff_tw  = (N_g - R_id_tw) + (G_t_tw - R_t_tw)
			nested_adj_tw = (id_nested_tw | t_nested_tw) ? 1 : 0
		}
		else {
			df_a_eff_tw = id_nested_tw ? 0 : N_g
			nested_adj_tw = id_nested_tw ? 1 : 0
		}
        
		if (ase_code == 1) {
			adj = N_clust / (N_clust - 1)
		}
		else {
			// small-sample denominator must stay positive
			real scalar tw_denom
			tw_denom = NT - nested_adj_tw - K - df_a_eff_tw
			if (tw_denom <= 0) {
				if (fatal_collinear) {
					errprintf("xtlp: vce(cluster " + clustervar + " " + clustervar2
					          + ") has too few degrees of freedom for the small-sample"
					          + " adjustment (denominator = " + strofreal(tw_denom)
					          + " <= 0).\n")
					errprintf("         Use vce(cluster " + clustervar + " "
					          + clustervar2 + ", ase) to drop the small-sample"
					          + " adjustment, or reduce regressors / add data.\n")
					exit(198)
				}
				st_matrix("b", J(1, K, .)); st_matrix("V", J(K, K, .))
				st_numscalar("N", N); st_numscalar("N_g", N_g)
				st_numscalar("df_r", .); st_numscalar("lag", .)
				return
			}
			adj = ((NT - 1) / tw_denom) * (N_clust / (N_clust - 1))
		}
		V_est = XX_inv * (M_tw * adj) * XX_inv

		// PSD fix (eigenvalue method, aligned to reghdfe_fix_psd).
		//
		// IMPORTANT: reghdfe standardizes the regressors/residual by their
		// column standard deviations before partialling out, builds the VCE
		// in that *standardized* space, applies the Cameron-Gelbach-Miller PSD
		// fix there, and only then undoes the standardization
		// (V = V :/ (stdev_x' * stdev_x)). The element-wise scaling
		// V[i,j] -> V[i,j]/(stx_i*stx_j) does NOT commute with the eigenvalue
		// zeroing, so doing PSD in the unscaled space (as a naive CGM impl
		// would) gives a different diagonal in the degenerate cases where V is
		// not PSD (twoway cluster with an absorbed FE nested in / cross-cutting
		// a cluster dim: A8/B8/E8/H2/H5). To align with reghdfe bit-for-bit we
		// therefore carry V into the standardized space, PSD it there, and
		// bring it back. The standardized meat/bread is algebraically identical
		// to the unscaled one (verified: max diff ~1e-14), so we only need to
		// change *where* PSD runs, not how the meat is built.
		//   stx = sx / sy,  sx = sd(raw X cols), sy = sd(raw Y)
		//   V_std = V_raw :* (stx' * stx)   (element-wise, since
		//          V_raw[i,j] = V_std[i,j] / (stx_i*stx_j))
		real rowvector sx_tw, stx_tw
		real scalar sy_tw
		sx_tw = sqrt((diagonal(cross(X, X))' :-
		             (cross(1, X):^2) / NT) :/ (NT - 1))
		sy_tw = sqrt((cross(Y, Y) :-
		             (cross(1, Y):^2) / NT) :/ (NT - 1))
		sx_tw = colmax((sx_tw \ J(1, K, 1e-3)))
		stx_tw = sx_tw :/ sy_tw

		_makesymmetric(V_est)
		// PSD in the standardized space, then undo the scaling.
		V_est = V_est :* (stx_tw' * stx_tw)
		symeigensystem(V_est, U_tw=., lam_tw=.)
		if (min(lam_tw) < 0) {
			printf("{txt}Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied. (h=%g)\n", st_numscalar("XTLP_H"))
			lam_tw = lam_tw :* (lam_tw :>= 0)
			V_est = quadcross(U_tw', lam_tw, U_tw')
		}
		V_est = V_est :/ (stx_tw' * stx_tw)
		df_r_val = N_clust - 1
		st_numscalar("N_clust_tw", N_clust)
	}
	else if (vce_code == 4) {
		// ============================================================
		//## Driscoll-Kraay (1998) SE, aligned to xtscc
		// ============================================================
		// Cross-sectional sums of moment conditions:
		//   h_t = sum_i d_it * e_it      (K x 1, not divided by sqrt(N))
		// Data is sorted by (ID, Time); group periods by Time value.
		real colvector tuniq
		real scalar T_eff, nObs, j, m_lag
		tuniq = uniqrows(sort(Time, 1))
		T_eff = rows(tuniq)
		nObs  = NT

		// Need >=2 time periods, else T_eff/(T_eff-1) divides by zero and the
		// HAC meat is undefined (a single cross-section has no autocovariance).
		if (T_eff <= 1) {
			if (fatal_collinear) {
				errprintf("xtlp: vce(dkraay) requires at least 2 time periods (T>=2);"
				          + " found " + strofreal(T_eff) + ".\n")
				exit(198)
			}
			st_matrix("b", J(1, K, .)); st_matrix("V", J(K, K, .))
			st_numscalar("N", N); st_numscalar("N_g", N_g)
			st_numscalar("df_r", .); st_numscalar("lag", .)
			return
		}

		// Time-continuity check (compromise: warn, do not exit)
		if (T_eff < (tuniq[T_eff] - tuniq[1] + 1)) {
			printf("{txt}warning: time variable is not regularly spaced (gaps detected).\n")
			printf("{txt}         DK lags are computed by row position, so results may be unreliable.\n")
		}

		// Bandwidth: -1 => Newey-West automatic, else user-specified lags
		if (dk_lag < 0) {
			m_lag = floor(4 * (T_eff / 100)^(2/9))
		}
		else {
			m_lag = dk_lag
		}
		if (m_lag < 0)  m_lag = 0
		if (m_lag >= T_eff) {
			// user-requested lag exceeds the available periods; cap at T_eff-1
			// (the maximum that fits the Bartlett kernel) and warn so the silent
			// truncation is visible.
			if (dk_lag >= 0) {
				printf("{txt}note: requested lag(" + strofreal(dk_lag)
				       + ") is too large for T=" + strofreal(T_eff)
				       + " periods; capped at " + strofreal(T_eff - 1) + ".\n")
			}
			m_lag = T_eff - 1
		}

		real matrix H
		H = J(T_eff, K, 0)
		real scalar tt
		real colvector idx_t
		real matrix d_t
		real colvector e_t
		for (tt = 1; tt <= T_eff; tt++) {
			idx_t = selectindex(Time :== tuniq[tt])
			d_t = d_dot[idx_t, .]          // N_t x K
			e_t = e[idx_t]                 // N_t x 1
			H[tt, .] = (d_t' * e_t)'       // 1 x K
		}

		// Meat with Bartlett (Newey-West) kernel, xtscc normalization
		//   Shat = (T/n^2) * [ sum_t h_t h_t' + sum_j w_j (Omegaj + Omegaj') ]
		real matrix Shat, Omegaj
		real scalar w_j
		Shat = cross(H, H) :/ ((nObs^2) / T_eff)        // Omega_0
		for (j = 1; j <= m_lag; j++) {
			w_j    = 1 - j / (m_lag + 1)
			Omegaj = cross(H[|j+1, 1 \ T_eff, .|], H[|1, 1 \ T_eff-j, .|]) :/ ((nObs^2) / T_eff)
			Shat   = Shat + w_j * (Omegaj + Omegaj')
		}

		// Sandwich + small-sample adjustment.
		// Two DoF conventions (see SE_summary.md / xtlp-dk-implementation memory):
		//   * default (nodfadj_code==0): reghdfe ivreg2-style, SUBTRACTS df_a
		//       V = (X'X)^-1 Shat (X'X)^-1 * (n^2/T) * (n-1)/(n-K-df_a) * T/(T-1)
		//     where df_a = absorbed-FE dof surviving nesting in the time cluster
		//     (dkraay clusters on Time, so the absorbed time FE is nested and fully
		//      relieved; only the individual FE counts, minus one redundant global
		//      constant when tfe). Matches reghdfe vce(dkraay) bit-for-bit.
		//   * nodfadj_code==1: xtscc-style, does NOT subtract df_a
		//       V = (X'X)^-1 Shat (X'X)^-1 * (n^2/T) * T/(T-1) * (n-1)/(n-rank)
		//     rank = K+1 (slope + one global constant). Matches xtscc bit-for-bit.
		// ase=1 drops ALL small-sample brackets -> only the n^2/T normalization
		//   remains (aligned to xtscc `ase`; reghdfe v6 has no asymptotic option).
		//   nodfadj is immaterial under ase (both df_a factors -> 1 asymptotically).
		real scalar dk_adj, df_a_dk
		real scalar id_nested_dk, t_nested_dk, G_t_dk, R_id_dk, R_t_dk
		// dkraay clusters on Time: the absorbed time FE is always nested in it
		// (t_nested_dk=1), and the individual FE is never nested in Time for a
		// panel (id is constant within individual, varies across time).
		id_nested_dk = 0
		if (fe_type == 2) {
			G_t_dk = rows(uniqrows(Time))
			t_nested_dk = 1
			R_t_dk  = t_nested_dk  ? G_t_dk : 1
			R_id_dk = id_nested_dk ? N_g    : (t_nested_dk ? 1 : 0)
			df_a_dk = (N_g - R_id_dk) + (G_t_dk - R_t_dk)
		}
		else {
			df_a_dk = id_nested_dk ? 0 : N_g
		}

		if (ase_code == 1) {
			// asymptotic: drop ALL small-sample brackets and keep only the
			// n^2/T normalization (aligned to xtscc `ase`; reghdfe v6 no longer
			// exposes an asymptotic option). nodfadj is immaterial here because
			// both df_a factors -> 1 asymptotically.
			dk_adj = (nObs^2) / T_eff
		}
		else {
			// small-sample factor divides by (n-K-df_a) or (n-rank); guard the
			// denominator so it cannot be <=0 (would yield missing SE).
			real scalar dk_denom
			dk_denom = (nodfadj_code == 1) ? (nObs - rank) : (nObs - K - df_a_dk)
			if (dk_denom <= 0) {
				if (fatal_collinear) {
					errprintf("xtlp: vce(dkraay) has too few degrees of freedom for the"
					          + " small-sample adjustment (denominator = "
					          + strofreal(dk_denom) + " <= 0).\n")
					errprintf("         Use vce(dkraay, ase) to drop the small-sample")
					errprintf(" adjustment, or reduce regressors / add data.\n")
					exit(198)
				}
				st_matrix("b", J(1, K, .)); st_matrix("V", J(K, K, .))
				st_numscalar("N", N); st_numscalar("N_g", N_g)
				st_numscalar("df_r", .); st_numscalar("lag", .)
				return
			}
			if (nodfadj_code == 1) {
				// xtscc-style: do NOT subtract df_a
				dk_adj = (nObs^2) / T_eff * (T_eff / (T_eff - 1)) * ((nObs - 1) / (nObs - rank))
			}
			else {
				// reghdfe ivreg2-style: subtract df_a
				dk_adj = (nObs^2) / T_eff * ((nObs - 1) / (nObs - K - df_a_dk)) * (T_eff / (T_eff - 1))
			}
		}
		V_est = XX_inv * Shat * XX_inv * dk_adj
		df_r_val = T_eff - 1
		st_numscalar("lag", m_lag)
	}

	real colvector se
	se = sqrt(diagonal(V_est))
	
    // -------------------------------------------------------
    // Step 6: Post Results to Stata
    // -------------------------------------------------------
	st_matrix("b", b_est')    // 1 x K
	st_matrix("V", V_est)     // K x K
	
	st_numscalar("N",    N)
	st_numscalar("N_g",  N_g)
	st_numscalar("df_r", df_r_val)
}

//# function: twoway demean
real matrix twoway_demean(real matrix Z, ///
						  real colvector ID, ///
                          real colvector Time, ///
						  real matrix info, ///
                          real scalar fe_type)
{
    real scalar i, start, end_t
    real matrix Z_dm, Z_prev

	real scalar iter, max_iter, tol, diff
    real colvector tuniq
    real scalar Nt, j, tval

    Z_dm = Z

	tol = 1e-9
    iter = 0
    max_iter = (fe_type == 1 ? 1 : 1000) 

    if (fe_type != 1) {
        tuniq = uniqrows(sort(Time, 1))
        Nt = rows(tuniq)
    }

    while (iter < max_iter) {

        Z_prev = Z_dm 

    // -------- 1. individual FE (within-id) --------
    	for (i = 1; i <= rows(info); i++) {
    	    start = info[i,1]
    	    end_t = info[i,2]

    	    real matrix Zi
    	    real rowvector mean_i

    	    Zi     = Z_dm[|start,1 \ end_t,.|]
    	    mean_i = mean(Zi)

    	    Z_dm[|start,1 \ end_t,.|] = Zi :- mean_i
    	}

    	// -------- 2. time FE (within-time) --------
        if (fe_type != 1) {
			for (j = 1; j <= Nt; j++) {
				tval = tuniq[j]

				real colvector idx_t
				idx_t = selectindex(Time :== tval)

				if (rows(idx_t) == 0) continue

				real matrix Zt
				real rowvector mean_t

                Zt     = Z_dm[idx_t, .]
				mean_t = mean(Zt)

				Z_dm[idx_t, .] = Zt :- mean_t
			}
		}

		iter++
        
        if (fe_type == 1) break 

        diff = max(abs(Z_dm - Z_prev))

        if (diff < tol) break
	}

	return(Z_dm)
}

end