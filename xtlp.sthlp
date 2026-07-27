{smcl}
{* *! version 1.2.0  26jul2026}{...}
{viewerjumpto "Syntax" "xtlp##syntax"}{...}
{viewerjumpto "Description" "xtlp##description"}{...}
{viewerjumpto "Options" "xtlp##options"}{...}
{viewerjumpto "Notes" "xtlp##notes"}{...}
{viewerjumpto "Examples" "xtlp##examples"}{...}
{viewerjumpto "Stored results" "xtlp##results"}{...}
{p2colset 1 15 17 2}{...}
{p2col:{bf:[XT] xtlp} {hline 2}}Panel local projections with fixed-effect (FE) estimator and split-panel jackknife (SPJ) estimator{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:xtlp} {depvar} {indepvars} {ifin} {cmd:,} {opt m:ethod(method_name)}
[{opt fe} {opt tfe} {cmd:vce(}{it:vcetype} [{cmd:,} {opt ase}]{cmd:)} {opt keepsin:gletons}{break}
{opt h:or(numlist)} {opt ytr:ansf(transf_name)}{break}
{opt sh:ock(integer)} {opt g:raph}]
{p_end}

{synoptset 23 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Estimation}
{synopt :{opt m:ethod(method_name)}}{cmd:fe} or {cmd:spj} (required){p_end}

{synopt :{opt fe}|{opt tfe}}include individual fixed effects (default) or two-way fixed effects (individual and time){p_end}

{synopt :{cmd:vce(}{it:vcetype} [{cmd:,} {opt ase}]{cmd:)}}variance estimator for the coefficients; default is {cmd:vce(unadjusted)}{p_end}

{synopt :{opt keepsin:gletons}}do not drop singleton groups before estimation{p_end}

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
{it:indepvars} may contain factor variables; see {help fvvarlist}.{p_end}
{pstd}
{it:depvar} and {it:indepvars} may contain time-series operators; see {help tsvarlist}.{p_end}

{pstd}
{cmd:xtlp} requires exactly one dependent variable and at least one independent
variable; if more than two variables are supplied, the first is taken as
{it:depvar} and the rest as {it:indepvars}.{p_end}

{pstd}
{cmd:xtlp} drops
observations with missing values in {it:depvar}, {it:indepvars}, the panel id,
or the time variable.{p_end}

{pstd}
{cmd:xtlp} requires Stata 14 or newer and does {it:not} provide
{cmd:e(predict)}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:xtlp} estimates the dynamic impulse response functions (IRFs) in panel data
using the Local Projection (LP) method. It offers two estimators via {opt method()}:
the standard fixed-effect estimator {cmd:m(fe)} and the split-panel
jackknife estimator {cmd:m(spj)}. The SPJ estimator addresses the intrinsic Nickell bias in dynamic settings ({help xtlp##MSS2026:Mei, Sheng, and Shi, 2026}).
{p_end}

{pstd}
When LPs are estimated with fixed effects in short panels, the dynamic structure
of the predictive equation induces the Nickell bias in the FE estimator, even if no lagged dependent variable appears explicitly in {it:indepvars}.
This bias invalidates standard inference based on the FE t-statistics. The SPJ estimator implemented here in this command provides a simple and effective bias-correction. It restores valid statistical inference in panel LPs, following
{help xtlp##MSS2026:Mei, Sheng, and Shi (2026)}.{p_end}

{pstd}
The command performs a single-equation estimation under the specified fixed-effect
structure ({opt fe} or {opt tfe}). Given {it:depvar} and {it:indepvars},
{cmd:xtlp} applies the chosen estimator ({cmd:m(fe)} or {cmd:m(spj)}) to produce coefficient estimates. The variance–covariance matrix (VCE) of the coefficients
is controlled by {opt vce()}. Five variance estimators are supported:
homoskedastic ({opt un:adjusted}, the default), heteroskedasticity-robust 
({opt r:obust}), one- and two-way cluster-robust ({opt cl:uster}), and 
Driscoll–Kraay ({opt dk:raay}).

{pstd}
For multiple horizons, {cmd:xtlp} automates the IRF construction over the range
specified in {opt h:or()}. It generates horizon-specific transformed dependent
variables via {opt ytr:ansf()}, runs a regression for each
horizon, and compiles the results. The option {opt sh:ock()} allows users to
treat several leading regressors as shocks; {cmd:xtlp} then reports the IRFs
and, if requested, produces IRF plots via {opt g:raph}.{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Estimation}

{phang}
{opt method(method_name)} is required and specifies the estimator. {it:method_name} may be:

{p2colset 9 21 23 2}{...}
{p2col:{cmd:fe}}requests the standard fixed-effects (within) estimator.{p_end}

{p2col:{cmd:spj}}requests the split-panel jackknife (SPJ) estimator. This method
splits each individual time series into two subpanels and combines the full-sample
and two half-sample FE estimates as {it:b_spj = 2*b_full - 0.5*(b_a + b_b)} to
deliver a bias-corrected estimator for dynamic panel LPs with fixed effects; see
{help xtlp##MSS2026:Mei, Sheng, and Shi (2026)}.{p_end}
{p2colreset}{...}

{phang}
{opt fe} includes individual fixed effects in the model. This is the default
if {opt tfe} is not specified. It cannot be combined with {opt tfe}.

{phang}
{opt tfe} includes two-way fixed effects (both individual and time fixed effects)
in the model, removed by iterative demeaning. It cannot be combined with {opt fe}.

{phang}
{cmd:vce(}{it:vcetype} [{cmd:,} {opt ase}]{cmd:)} specifies the standard-error estimator. {it:vcetype} may be:

{p2colset 9 21 23 2}{...}
{p2col:{opt un:adjusted}}conventional (homoskedastic) standard errors; the default. {p_end}
{p 20 22 2}Aligned with {cmd:xtreg ..., fe} and {cmd:reghdfe} without {cmd:vce()}.{p_end}

{p2col:{opt r:obust}}heteroskedasticity-robust (HC1) standard errors.{p_end}
{p 20 22 2}Aligned with {cmd:reghdfe ..., vce(robust)}.{p_end}

{p2col:{opt cl:uster} {it:clustervar1} [{it:clustervar2}]}{p_end}
{p 20 22 2}one- or two-way
cluster-robust standard errors; {it:clustervar1} is required and
{it:clustervar2} is required for two-way clustering (at most two clustervars). Two-way clustering follows {help xtlp##CGM2011:Cameron, Gelbach, and Miller (2011)}.{p_end}
{p 20 22 2}Aligned with {cmd:reghdfe ..., vce(cluster ...)}. {p_end}

{p2col:{opt dk:raay} [{cmd:lag(}{it:#}{cmd:)}] [, {opt nodfadj}]}{p_end}
{p 20 22 2}Driscoll–Kraay ({help xtlp##DK1998:1998})
standard errors, robust to cross-sectional and serial correlation, clustered on
the time variable. If gaps are detected in the time variable, {cmd:xtlp} issues a warning but does {it:not} abort execution. {cmd:lag(#)} sets the maximum lag order; the default is
{it:floor(4(T/100)^(2/9))}. {cmd:nodfadj} switches the small-sample adjustment
from the {helpb reghdfe} style (default, subtracts the absorbed-FE
degrees of freedom) to the {helpb xtscc} style. {cmd:nodfadj} and {cmd:ase} are mutually exclusive.{p_end}
{p 20 22 2}{cmd:xtlp ..., vce(dkraay lag(2))} is aligned with {cmd:reghdfe ..., vce(dkrray 3)}, since {helpb reghdfe}'s vce(dkraay #) uses bandwidth = lags + 1; {cmd:xtlp ..., vce(dkraay lag(2), nodfadj)} is aligned with {cmd:xtscc ..., lag(2)}. {p_end}

{p2col:{cmd:[,} {opt ase}{cmd:]}}returns asymptotic standard
errors with {it:no} small-sample adjustment. Allowed with all {it:vcetype}s. {cmd:ase} and {cmd:nodfadj} are mutually exclusive.{p_end}
{p 20 22 2}{cmd:xtlp ..., vce(dkraay lag(2), ase)} is aligned with {cmd:xtscc ..., lag(2) ase}.{p_end}
{p2colreset}{...}

{phang}
{opt keepsingletons} requests that singleton groups {it:not} be dropped before
estimation. The default is to drop them, aligned with {help xtlp##COR2015:Correia (2015)} and {helpb reghdfe}.

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
The start horizon {it:S} must be 0 or 1, and {it:H} must be
greater than {it:S}; at most two integers are accepted.{p_end}

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
are treated as shocks when constructing IRFs. The argument must be a single
integer between 1 and {it:K} (the number of covariates in {it:indepvars}); other
values are an error. The default is {cmd:shock(1)}.
For example, {cmd:shock(2)} means the first two variables in {it:indepvars}
are treated as separate shocks, and the command reports an IRF for each of them.{p_end}

{phang}
{opt graph} requests that IRFs be graphed after estimation. For each shock,
the graph plots the point estimates together with 95% confidence intervals over the specified horizons. {p_end}


{marker notes}{...}
{title:Notes}

{pstd} If the regressors are collinear or the sample is degenerate, a single estimation exits with error 198. In a multi-horizon run, a degenerate horizon does {it:not} abort execution: {cmd:xtlp} posts missing coefficients and standard errors for that horizon and continues with the remaining horizons. {p_end}


{marker examples}{...}
{title:Examples}

{pstd}Download four {it:.dta} files from the {cmd:applications/data_preparation} folder in the {browse "https://github.com/metricshilab/panel-lp-replication":replication package} of {help xtlp##MSS2026:Mei, Sheng, and Shi (2026)}{p_end}

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

{pstd}Specify two shock variables using {cmd:sh(2)}{p_end}
{phang2}{stata "xtlp F1y L0HHD_L1GDP L0NFD_L1GDP L1HHD_L1GDP L2HHD_L1GDP L3HHD_L1GDP L4HHD_L1GDP L1NFD_L1GDP L2NFD_L1GDP L3NFD_L1GDP L4NFD_L1GDP L0y L1y L2y L3y L4y, fe m(fe) h(0 9) sh(2) g"}{p_end}
{phang2}{stata "xtlp F1y L0HHD_L1GDP L0NFD_L1GDP L1HHD_L1GDP L2HHD_L1GDP L3HHD_L1GDP L4HHD_L1GDP L1NFD_L1GDP L2NFD_L1GDP L3NFD_L1GDP L4NFD_L1GDP L0y L1y L2y L3y L4y, fe m(spj) h(0 9) sh(2) g"}{p_end}

    {title:Example 6: Variance estimators ({opt vce()})}

{phang2}{stata "use RR_f4data, replace"}{p_end}

{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, fe m(fe) vce(un)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(r)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(cl COUNTDUMS)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(cl COUNTDUMS halfyear)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(dk)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(dk, nodfadj)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(dk, ase)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(dk lag(2))"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(dk lag(2), nodfadj)"}{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(dk lag(2), ase)"}{p_end}

    {title:Example 7: Keeping singleton observations ({opt keepsingletons})}

{phang2}{stata "use RR_f4data, replace"}{p_end}

{pstd}By default singleton groups are dropped;
{opt keepsingletons} retains them and warns{p_end}
{phang2}{stata "xtlp f10LNGDP CRISIS l1LNGDP l2LNGDP l3LNGDP l4LNGDP l1CRISIS l2CRISIS l3CRISIS l4CRISIS, tfe m(fe) vce(robust) keepsin"}{p_end}


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
{synopt:{cmd:e(N_clust)}}number of clusters; only with {cmd:vce(cluster ...)}{p_end}
{synopt:{cmd:e(lag)}}Driscoll–Kraay maximum lag order; only with {cmd:vce(dkraay)}{p_end}

{p2col 5 20 24 2:Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtlp}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(vce)}}canonical {cmd:vce()} string, e.g. {cmd:unadjusted}, {cmd:robust}, {cmd:cluster clustervar1}, {cmd:cluster clustervar1 clustervar2}, {cmd:dkraay lag(2)}, or {cmd:dkraay lag(2) nodfadj}. The {cmd:ase} suboption is {it:not} recorded in {cmd:e(vce)}{p_end}
{synopt:{cmd:e(vcetype)}}title displayed above the standard errors: {cmd:Robust}, {cmd:Cluster}, or {cmd:Drisc/Kraay}{p_end}
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
{p2col 5 20 24 2:Scalars}{p_end}
{synopt:{cmd:e(df_r}{it:h}{cmd:)}}residual degrees of freedom for horizon {it:h}{p_end}
{synopt:{cmd:e(lag}{it:h}{cmd:)}}Driscoll–Kraay maximum lag order for horizon {it:h}; only with {cmd:vce(dkraay)}{p_end}

{p2col 5 20 24 2:Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtlp}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(vce)}}canonical {cmd:vce()} string (see Case 1){p_end}
{synopt:{cmd:e(vcetype)}}title displayed above the standard errors (see Case 1){p_end}

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
{p2col:{cmd:"Std. err."}}standard error of IRF{p_end}
{p2col:{cmd:"95% CI Lower"}}lower 95% confidence interval{p_end}
{p2col:{cmd:"95% CI Upper"}}upper 95% confidence interval{p_end}
{p2colreset}{...}

{pstd}
When {cmd:shock(#)} specifies more than one shock (i.e., {it:#} > 1), {cmd:e(irf)}
is built by horizontally concatenating one 4-column block per shock. Each block
reuses the same four column names ({cmd:"IRF"}, {cmd:"Std. err."},
{cmd:"95% CI Lower"}, {cmd:"95% CI Upper"}), so {cmd:e(irf)} has {it:#} sets of
repeatedly-named columns; the
shock order follows the first {it:#} variables of {it:indepvars}.{p_end}


{marker references}{...}
{title:References}

{marker CGM2011}{...}
{phang}
Cameron, A. C., Gelbach, J. B., and Miller, D. L. (2011). Robust inference with multiway clustering. {it:Journal of Business & Economic Statistics}, 29(2), 238–249.{p_end}

{marker COR2015}{...}
{phang}
Correia, S. (2015). Singletons, cluster-robust standard errors and fixed effects: A bad mix.{p_end}

{marker DK1998}{...}
{phang}
Driscoll, J. C., and Kraay, A. C. (1998). Consistent covariance matrix estimation with spatially dependent panel data. {it:Review of Economics and Statistics}, 80(4), 549–560.{p_end}

{marker MSS2026}{...}
{phang}
Ziwei Mei, Liugang Sheng, Zhentao Shi (2026). {browse "https://doi.org/10.1016/j.jinteco.2025.104210":Nickell bias in panel local projection: Financial crises are worse than you think}. {it:Journal of International Economics}, 104210.{p_end}
{phang}
{browse "https://github.com/metricshilab/panel-lp-replication":Replication package} for {help xtlp##MSS2026:Mei, Sheng, and Shi (2026)}. {p_end}


{marker author}{...}
{title:Author}

{pstd}
Shu SHEN{break}
shushen@link.cuhk.edu.hk
{p_end}