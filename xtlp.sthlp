{smcl}
{* *! version 1.0.0  10dec2025}{...}
{viewerjumpto "Syntax" "xtspj##syntax"}{...}
{viewerjumpto "Description" "xtspj##description"}{...}
{viewerjumpto "Options" "xtspj##options"}{...}
{viewerjumpto "Examples" "xtspj##examples"}{...}
{viewerjumpto "Stored results" "xtspj##results"}{...}
{p2colset 1 15 17 2}{...}
{p2col:{bf:[XT] xtlp} {hline 2}}Panel local projections with fixed-effect (FE) estimator and split-panel jackknife (SPJ) estimator{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:xtlp} {depvar} {indepvars} {ifin}
{cmd:,} {opt m:ethod(method_name)} [{opt fe} {opt tfe} {opt h:or(numlist)} {opt ytr:ansf(transf_name)} {opt sh:ock(integer)} {opt g:raph}]

{synoptset 23 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Estimation}
{synopt :{opt m:ethod(method_name)}}{cmd:fe} or {cmd:spj} (required){p_end}

{synopt :{opt fe}|{opt tfe}}include individual fixed effects (default) or two-way fixed effects (individual and time){p_end}

{syntab:Multiple Horizons}
{synopt :{opt h:or(numlist)}}horizon(s) for impulse response functions: specify {it:#} for max horizon or {it:#_start #_end}{p_end}
{synopt :{opt ytr:ansf(transf_name)}}transform dependent variable: {cmd:level} (default), {cmd:diff}, {cmd:cmltdiff}, or {cmd:cmltsum}{p_end}
{synopt :{opt sh:ock(integer)}}number of leading variables in {it:indepvars} to treat as shocks; default is {cmd:shock(1)}{p_end}
{synopt :{opt g:raph}}graph the impulse response functions{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
A panel variable and a time variable must be specified using {helpb xtset}.{p_end}

{pstd}
{it:depvar} and {it:indepvars} may not contain factor variables and time-series operators; see {help fvvarlist} and {help tsvarlist}.{p_end}

{pstd}
The command requires exactly one dependent variable and at least one independent variable.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:xtlp} estimates the dynamic impulse response functions (IRFs) in panel data 
using the Local Projection (LP) method. It offers two estimators via {opt method()}:
the standard fixed-effect estimator ({cmd:m(fe)}) and the split-panel 
jackknife estimator ({cmd:m(spj)}). The SPJ estimator addresses the intrinsic Nickell bias in dynamic settings ({help xtlp##MSS2025:Mei, Sheng, and Shi, 2025}).
{p_end}

{pstd}
When LPs are estimated with fixed effects in short panels, the dynamic structure
of the predictive equation induces the Nickell bias in the FE estimator, even if no lagged dependent variable appears explicitly in {it:indepvars}.
This bias invalidates standard inference based on the FE t-statistics. The SPJ estimator implemented here in this command provides a simple and effective bias-correction. It restores valid statistical inference in panel LPs, following
{help xtlp##MSS2025:Mei, Sheng, and Shi (2025)}.{p_end}

{pstd}
The command performs a single-equation estimation under the specified fixed-effect 
structure ({opt fe} or {opt tfe}). Given {it:depvar} and {it:indepvars},
{cmd:xtlp} applies the chosen estimator ({cmd:m(fe)} or {cmd:m(spj)}) to produce coefficient 
estimates.
{p_end}

{pstd}
For multiple horizons, {cmd:xtlp} automates the IRF construction over the range
specified in {opt hor()}. It generates horizon-specific transformed dependent
variables via {opt ytransf()}, runs a regression for each
horizon, and compiles the results. The option {opt shock()} allows users to
treat several leading regressors as shocks; {cmd:xtlp} then reports the IRFs
and, if requested, produces IRF plots via {opt graph}.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimation}

{phang}
{opt method(method_name)} is required and specifies the estimator.

{p2colset 9 21 23 2}{...}
{p2col:{cmd:fe}}requests the standard fixed-effects (within) estimator.{p_end}

{p2col:{cmd:spj}}requests the split-panel jackknife (SPJ) estimator. This method
splits each individual time series into two subpanels and combines the full-sample
and subsample FE estimates to deliver a bias-corrected estimator for dynamic
panel LPs with fixed effects; see {help xtlp##MSS2025:Mei, Sheng, and Shi (2025)}.{p_end}
{p2colreset}{...}

{phang}
{opt fe} includes individual fixed effects in the model. This is the default
if {opt tfe} is not specified. It cannot be combined with {opt tfe}.

{phang}
{opt tfe} includes two-way fixed effects (both individual and time fixed effects) in the model. It cannot be combined with {opt fe}.

{dlgtab:Multiple horizons}

{phang}
{opt hor(numlist)} specifies the horizons for the LPs. This option accepts either one or two integers. The default is {cmd:hor(0)}. 

{phang2}
If {cmd:hor(0)} is specified (or implied by default), only a single estimation
is performed, and horizon-specific IRF options (i.e., {opt ytransf()}, {opt shock()}, {opt graph}) do not apply.{p_end}

{phang2}
If one integer {it:H} is specified (e.g., {cmd:hor(5)}), LPs are estimated for horizons 0 to {it:H}.

{phang2}
If two integers {it:S} and {it:H} are specified (e.g., {cmd:hor(1 5)}), LPs are estimated for horizons {it:S} to {it:H}. 
The start horizon {it:S} must be 0 or 1.

{phang}
{opt ytransf(transf_name)} specifies the transformation applied to the dependent variable {it:depvar} for the LP at each horizon {it:h}.

{p2colset 9 21 23 2}{...}
{p2col:{cmd:level}}(default) uses the level of {it:depvar}, {cmd:{it:y_{i,t+h}}}, as the dependent variable.{p_end}

{p2col:{cmd:diff}}uses the first difference of {it:depvar}, {cmd:{it:y_{i,t+h} - y_{i,t+h-1}}}, as the dependent variable.{p_end}

{p2col:{cmd:cmltdiff}}uses the cumulative difference of {it:depvar} relative to period t-1, {cmd:{it:y_{i,t+h} - y_{i,t-1}}}, as the dependent variable, which captures the cumulative response of {it:depvar}.{p_end}

{p2col:{cmd:cmltsum}}uses the cumulative sum of {it:depvar}, {cmd:{it:Σ_{k=0}^h y_{i,t+k}}}, as the dependent variable.{p_end}
{p 20 22 2}Note: This option is typically useful when {it:depvar} is already a
first-differenced variable (e.g., growth rate), so that the cumulative
sum recovers the level impact over the horizon.{p_end}
{p2colreset}{...}

{phang}
{opt shock(integer)} specifies that the first {it:#} variables in {it:indepvars}
are treated as shocks when constructing IRFs. The default is {cmd:shock(1)}.
For example, {cmd:shock(2)} means the first two variables in {it:indepvars}
are treated as separate shocks, and the command reports an IRF for each of them.{p_end}

{phang}
{opt graph} requests that IRFs be graphed after estimation. For each shock,
the graph plots the point estimates together with 95% confidence intervals
over the specified horizons.{p_end}

{dlgtab:Other}

{pstd}
The current version of {cmd:xtlp} does not accept a user-specified
{cmd:vce()} option. The variance–covariance matrix of the coefficients is
computed using a panel-robust sandwich estimator with clustering at the
individual level.


{marker examples}{...}
{title:Examples}

{pstd}Download four {it:.dta} files from the {cmd:applications/data_preparation} folder in the {browse "https://github.com/metricshilab/panel-lp-replication":replication package} of {help xtlp##MSS2025:Mei, Sheng, and Shi (2025)}{p_end}

{phang2}{it:./applications/data_preparation/RR_f4data.dta }{p_end}
{phang2}{it:./applications/data_preparation/BVX_t1data.dta}{p_end}
{phang2}{it:./applications/data_preparation/MSV_f2data.dta}{p_end}
{phang2}{it:./applications/data_preparation/CS_f3data.dta }{p_end}


    {title:Example 1: FE vs. SPJ with {opt fe} (single estimation)}

{phang2}{stata "use BVX_t1data, clear"}{p_end}
{phang2}{stata "keep if smp==1"}{p_end}

{pstd}Estimate using FE estimator ({cmd:m(fe)}) with individual fixed effects ({opt fe}){p_end}
{phang2}
{stata "xtlp Fd6y R_B L1R_B L2R_B L3R_B R_N L1R_N L2R_N L3R_N D1y L1D1y L2D1y L3D1y D1d_y L1D1d_y L2D1d_y L3D1d_y, fe m(fe)"}
{p_end}

{pstd}Estimate using SPJ estimator ({cmd:m(spj)}) with individual fixed effects ({opt fe}){p_end}
{phang2}
{stata "xtlp Fd6y R_B L1R_B L2R_B L3R_B R_N L1R_N L2R_N L3R_N D1y L1D1y L2D1y L3D1y D1d_y L1D1d_y L2D1d_y L3D1d_y, fe m(spj)"}
{p_end}

    {title:Example 2: FE vs. SPJ with ({opt tfe}) (single estimation)}

{phang2}{stata "use RR_f4data, replace"}{p_end}

{pstd}Estimate using FE estimator ({cmd:m(fe)}) with two-way fixed effects ({opt tfe}){p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe)"}{p_end}

{pstd}Estimate using SPJ estimator ({cmd:m(spj)}) with two-way fixed effects ({opt tfe}){p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(spj)"}{p_end}

    {title:Example 3: Estimating IRFs (multiple horizons)}

{phang2}{stata "use RR_f4data, replace"}{p_end}

{pstd}Estimate IRF from horizon 0 to 10 ({cmd:h(0 10)}) and plot graph ({cmd:g}){p_end}
{phang2}{stata "xtlp f0LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) h(0 10) g"}{p_end}
{phang2}{stata "xtlp f0LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(spj) h(0 10) g"}{p_end}

{pstd}Estimate IRF from horizon 1 to 10 ({cmd:h(1 10)}) and plot graph ({cmd:g}){p_end}
{phang2}{stata "xtlp f0LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) h(1 10) g"}{p_end}
{phang2}{stata "xtlp f0LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(spj) h(1 10) g"}{p_end}

    {title:Example 4: Dependent variable transformation (multiple horizons)}

{phang2}{stata "use CS_f3data, clear"}{p_end}

{pstd}Using cumulative sum transformation ({cmd:cmltsum}) for growth rates{p_end}
{phang2}{stata "xtlp GRRT_WB CRISIS l1CRISIS l2CRISIS l3CRISIS l4CRISIS l1GRRT_WB l2GRRT_WB l3GRRT_WB l4GRRT_WB, fe m(fe) h(0 10) ytr(cmltsum) g"}{p_end}
{phang2}{stata "xtlp GRRT_WB CRISIS l1CRISIS l2CRISIS l3CRISIS l4CRISIS l1GRRT_WB l2GRRT_WB l3GRRT_WB l4GRRT_WB, fe m(spj) h(0 10) ytr(cmltsum) g"}{p_end}

    {title:Example 5: Multiple shocks (multiple horizons)}

{phang2}{stata "use MSV_f2data, clear"}{p_end}
{phang2}{stata "keep CountryCode year F1y F2y F3y F4y F5y F6y F7y F8y F9y F10y L0HHD_L1GDP L1HHD_L1GDP L2HHD_L1GDP L3HHD_L1GDP L4HHD_L1GDP L0NFD_L1GDP L1NFD_L1GDP L2NFD_L1GDP L3NFD_L1GDP L4NFD_L1GDP L0y L1y L2y L3y L4y"}{p_end}

{pstd}Specifying two shock variables using {cmd:sh(2)}{p_end}
{phang2}{stata "xtlp F1y L0HHD_L1GDP L0NFD_L1GDP L1HHD_L1GDP L2HHD_L1GDP L3HHD_L1GDP L4HHD_L1GDP L1NFD_L1GDP L2NFD_L1GDP L3NFD_L1GDP L4NFD_L1GDP L0y L1y L2y L3y L4y, fe m(fe) h(0 9) sh(2) g"}{p_end}
{phang2}{stata "xtlp F1y L0HHD_L1GDP L0NFD_L1GDP L1HHD_L1GDP L2HHD_L1GDP L3HHD_L1GDP L4HHD_L1GDP L1NFD_L1GDP L2NFD_L1GDP L3NFD_L1GDP L4NFD_L1GDP L0y L1y L2y L3y L4y, fe m(spj) h(0 9) sh(2) g"}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:xtlp} is an {cmd:eclass} command. The contents of {cmd:e()} depend on
the {cmd:hor()} setting and the {cmd:shock()} setting.

{dlgtab:Case 1: Single estimation}

{pstd}
When {cmd:hor(0)} is specified (or implied by default), {cmd:xtlp} runs one FE or SPJ estimation.{p_end}

{pstd}
{cmd:xtlp} posts results in a way similar to other linear regression commands.{p_end}

{pstd}
{cmd:xtlp} stores the following in {cmd:e()}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2:Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of panels (individuals){p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}

{p2col 5 20 24 2:Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtlp}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 20 24 2:Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector (1 × K){p_end}
{synopt:{cmd:e(V)}}variance–covariance matrix (K × K){p_end}

{p2col 5 20 24 2:Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{dlgtab:Case 2: Multiple horizons}

{pstd}
When {cmd:hor()} specifies multiple horizons, {cmd:xtlp} estimates the
model separately for each horizon {it:h}.{p_end}

{pstd}
{cmd:xtlp} posts a consolidated matrix of IRFs, including point estimates,
standard errors, and the lower and upper bounds of 95% confidence intervals.{p_end}

{pstd}
{cmd:xtlp} stores the following in {cmd:e()}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2:Matrices}{p_end}
{synopt:{cmd:e(b}{it:h}{cmd:)}}coefficient vector for horizon {it:h} (1 × K){p_end}
{synopt:{cmd:e(V}{it:h}{cmd:)}}variance–covariance matrix for horizon {it:h} (K × K){p_end}
{synopt:{cmd:e(irf)}}IRF results: estimate, standard error, lower and upper bounds (H × 4 or H × (4×#shocks)){p_end}
{p2colreset}{...}

{pstd}
When {cmd:shock(1)} is specified (or implied by default), the columns of
{cmd:e(irf)} are named{p_end}

{p2colset 7 27 29 2}{...}
{p2col:{cmd:"IRF"}}IRF point estimate{p_end}
{p2col:{cmd:"Std.Err."}}standard error of IRF{p_end}
{p2col:{cmd:"IRF LOW"}}lower 95% confidence interval{p_end}
{p2col:{cmd:"IRF UP"}}upper 95% confidence interval{p_end}
{p2colreset}{...}

{pstd}
When {cmd:shock(#)} specifies more than one shock (i.e., {it:#} > 1), {cmd:e(irf)}
is organized as blocks of four columns for each shock with names{p_end}

{p2colset 7 27 29 2}{...}
{p2col:{cmd:"IRF_#"}}IRF point estimate for shock {it:#}{p_end}
{p2col:{cmd:"Std.Err._#"}}standard error of IRF for shock {it:#}{p_end}
{p2col:{cmd:"IRF LOW_#"}}lower 95% confidence interval for shock {it:#}{p_end}
{p2col:{cmd:"IRF UP_#"}}upper 95% confidence interval for shock {it:#}{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker MSS2025}{...}
{phang}
Ziwei Mei, Liugang Sheng, Zhentao Shi, 2025,
{browse "https://arxiv.org/abs/2302.13455":"{it:Nickell bias in panel local projection: Financial crises are worse than you think}"}, Journal of International Economics.{p_end}

{phang}
{browse "https://github.com/metricshilab/panel-lp-replication":Replication package} for Ziwei Mei, Liugang Sheng, Zhentao Shi, 2025, 
"{it:Nickell bias in panel local projection: Financial crises are worse than you think}", Journal of International Economics.{p_end}

{marker author}{...}
{title:Author}

{pstd}
Shu SHEN{break}
shushen@link.cuhk.edu.hk
{p_end}