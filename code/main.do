/*
================================================================================
DO FILE: Financial Inclusion and Consumption Smoothing under External Income Shocks
AUTHOR:  Malik Jibran Tariq
DATE:    May 2026
CONTEXT: Columbia SIPA, MPA in Economic Policy Management

USAGE:   Run from the repository root. Reads data/panel_FINAL.csv and
         writes tables, figures and the log to output/.

DATA:    panel_FINAL.csv
         Source: World Bank WDI + IMF Financial Access Survey + WB Pink Sheet
         9 countries · 2004–2024 · Annual · Balanced raw panel (189 obs)

VARIABLE DEFINITIONS (as in CSV):
  country_id   Numeric country identifier (for xtset)
  iso3         ISO3 country code
  country      Country name
  year         Calendar year
  cons_growth  Household final consumption per capita growth, % real (WDI)
  tot_index    Terms of trade index, 2015=100 (WDI)
  inflation    CPI inflation, annual % (WDI)
  trade_openness (Exports+Imports)/GDP, % (WDI)
  exch_rate    Official exchange rate, LCU per USD (WDI)
  gdp_pc       GDP per capita, constant USD (WDI)
  oil_price    Brent crude, USD/bbl annual avg (WB Pink Sheet)
  FA30         Deposit accounts, commercial banks per 1,000 adults (IMF FAS)
  FDEPOID      Depositors, credit unions & cooperatives per 1,000 adults (FAS)
  FMOBAB       Registered mobile money accounts per 1,000 adults (IMF FAS)
  fi_composite FI composite = FA30 + FDEPOID + FMOBAB
 */

clear all
capture log close
capture mkdir "output"
log using "output/main_log.smcl", replace


set more off


* Required packages (uncomment on first run)
* ssc install reghdfe,  replace
* ssc install ftools,   replace
* ssc install estout,   replace
* ssc install winsor2,  replace
* ssc install coefplot, replace
* ssc install boottest, replace


/*==============================================================================
   SECTION 1: LOAD DATA
==============================================================================*/

import delimited "data/panel_FINAL.csv", clear varnames(2) case(lower) bindquote(nobind)

summarize
describe
list in 1/5

list country year cons_growth tot_index fi_composite in 1/10

label var cons_growth    "HH consumption per capita growth, % real (WDI NE.CON.PRVT.PC.KD.ZG)"
label var tot_index      "Terms of trade index, 2015=100 (WDI TT.PRI.MRCH.XD.WD)"
label var inflation      "CPI inflation, annual % (WDI FP.CPI.TOTL.ZG)"
label var trade_openness "Trade openness, % of GDP (WDI NE.TRD.GNFS.ZS)"
label var exch_rate      "Exchange rate, official LCU/USD (WDI PA.NUS.FCRF)"
label var gdp_pc         "GDP per capita, constant USD (WDI NY.GDP.PCAP.KD)"
label var oil_price      "Brent crude, USD/bbl annual avg (WB Pink Sheet)"
label var fa30           "Deposit accounts, commercial banks per 1,000 adults (IMF FAS)"
label var fdepoid        "Depositors, credit unions & cooperatives per 1,000 adults (IMF FAS)"
label var fmobab         "Registered mobile money accounts per 1,000 adults (IMF FAS)"
label var fi_composite   "FI composite = FA30 + FDEPOID + FMOBAB"


/*==============================================================================
   SECTION 2: DECLARE PANEL STRUCTURE
==============================================================================*/

xtset country_id year


xtdescribe

/*==============================================================================
   SECTION 3: SAMPLE AUDIT
  
==============================================================================*/

di "=== SAMPLE AUDIT ==="
misstable summarize cons_growth tot_index fi_composite inflation trade_openness


gsort -fi_composite
list country year fi_composite in 1/15

tab iso3

xtsum cons_growth tot_index fi_composite

/*==============================================================================
   SECTION 4: VARIABLE CONSTRUCTION
==============================================================================*/

/* ----- 4.1 TERMS OF TRADE SHOCK ----- */

gen ln_tot = ln(tot_index)


gen shock = 100*D.ln_tot
summ shock
   

label var ln_tot "ln(ToT index)"
label var shock  "ToT shock: 100 × Δln(ToT index), pp"

/* ----- 4.2 OIL SHOCK ----- */

gen ln_oil = ln(oil_price)
gen shock_oil = 100*(ln_oil - L.ln_oil)

label var ln_oil    "ln(Brent crude price)"
label var shock_oil "Oil shock: 100 × Δln(oil price), global, absorbed by year FE"

/* ----- 4.3 WINSORIZED SHOCK ----- */

winsor2 shock, suffix(_w) cuts(1 99)

label var shock_w "ToT shock, winsorized at p1/p99"

/* ----- 4.4 FI COMPOSITE AND TRANSFORMATIONS ----- */


sum fi_composite, detail

gen ln_fi = ln(fi_composite)

summ ln_fi, detail

gen ln_fi_lag1 = L.ln_fi

label var ln_fi      "ln(FI composite)"
label var ln_fi_lag1 "ln(FI composite), t-1 — pre-existing conditions"

/* ----- 4.5 STANDARDIZE FI ON ESTIMATION SAMPLE ----- */


sum fi_composite if !missing(cons_growth) & !missing(shock) & ///
    !missing(fi_composite) & !missing(inflation) & !missing(trade_openness)
local fi_mean = r(mean)
local fi_sd   = r(sd)
gen fi_std = (fi_composite - `fi_mean') / `fi_sd'

label var fi_std "FI composite, standardized (mean=0, SD=1) on estimation sample"

/* Store percentiles of fi_std for marginal effects */
sum fi_std if !missing(cons_growth) & !missing(shock) & ///
    !missing(fi_composite) & !missing(inflation) & !missing(trade_openness), detail
local fi_p10 = r(p10)
local fi_p25 = r(p25)
local fi_p50 = r(p50)
local fi_p75 = r(p75)
local fi_p90 = r(p90)
di "FI std percentiles: p10=`fi_p10' p25=`fi_p25' p50=`fi_p50' p75=`fi_p75' p90=`fi_p90'"

/* Standardize bank-only FI (FA30 alone) for robustness Model M8 */
sum fa30 if !missing(cons_growth) & !missing(shock) & ///
    !missing(fa30) & !missing(inflation) & !missing(trade_openness)
local fa30_mean = r(mean)
local fa30_sd   = r(sd)
gen fi_bank_std = (fa30 - `fa30_mean') / `fa30_sd'

label var fi_bank_std "FA30 only, standardized — bank-only robustness"

/* ----- 4.6 INTERACTION TERMS ----- */

gen shock_x_fi_std    = shock   * fi_std
gen shock_x_lnfi      = shock   * ln_fi
gen shock_x_lnfi_lag  = shock   * ln_fi_lag1
gen shock_w_x_fi_std  = shock_w * fi_std
gen oil_x_fi_std      = shock_oil * fi_std
gen shock_x_fi_bank   = shock   * fi_bank_std

label var shock_x_fi_std   "SHOCK × fi_std          ← MAIN INTERACTION"
label var shock_x_lnfi     "SHOCK × ln(FI)"
label var shock_x_lnfi_lag "SHOCK × ln(FI) at t-1"
label var shock_w_x_fi_std "SHOCK_w × fi_std"
label var oil_x_fi_std     "SHOCK_oil × fi_std"
label var shock_x_fi_bank  "SHOCK × fi_bank_std"

/* ----- 4.7 CONTROL VARIABLES ----- */



gen crisis_year = (year==2008 | year==2009 | year==2020)
label var crisis_year "=1 for GFC (2008,2009) and COVID (2020)"


/*==============================================================================
   SECTION 5: PRE-REGRESSION DIAGNOSTICS
==============================================================================*/

/* ----- 5.1 WITHIN-COUNTRY VARIATION (identification prerequisite) ----- */

di "=== WITHIN-COUNTRY VARIATION ==="
xtsum shock fi_composite fi_std cons_growth


/* ----- 5.2 MULTICOLLINEARITY CHECK ----- */

di "=== MULTICOLLINEARITY CHECK ==="
reg cons_growth shock fi_std shock_x_fi_std inflation trade_openness ///
    if !missing(cons_growth) & !missing(shock) & !missing(fi_std)
vif

pwcorr shock fi_std shock_x_fi_std inflation trade_openness, star(0.05)


/*==============================================================================
   SECTION 6: MAIN REGRESSIONS
   All main models use xtreg, fe with country fixed effects and year dummies.
   Standard errors clustered at country level throughout.
==============================================================================*/



eststo clear


/*------------------------------------------------------------------------------
   MODEL M1: BASELINE — No interaction
   Purpose: Establishes the direct effects of the shock and FI before
   adding the interaction. Provides a reference point and tests whether
   FI has a direct effect on consumption growth (β₂).
------------------------------------------------------------------------------*/


*========================================================
* MODEL M1
*========================================================

* M1: Baseline, no interaction
eststo M1: xtreg cons_growth ///
    shock fi_std inflation trade_openness ///
    i.year, fe cluster(country_id)

estadd local cfe "Yes"
estadd local yfe "Yes"
estadd local clust "Country (N=9)"


/*------------------------------------------------------------------------------
   MODEL M2: MAIN SPECIFICATION — SHOCK × FI (standardized)
   This is the primary model. β₃ on shock_x_fi_std is the main coefficient
   of interest.
*/

*========================================================
* MODEL M2: MAIN SPECIFICATION
*========================================================

* M2: Main interaction
eststo M2: xtreg cons_growth ///
    c.shock##c.fi_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

estadd local cfe "Yes"
estadd local yfe "Yes"
estadd local clust "Country (N=9)"

/*------------------------------------------------------------------------------
   MODEL M3: LOG FI INTERACTION
   Tests whether the functional form of FI matters.
*/

*========================================================
* MODEL M3: LOG FI
*========================================================

* M3: Log FI interaction
eststo M3: xtreg cons_growth ///
    c.shock##c.ln_fi ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

estadd local cfe "Yes"
estadd local yfe "Yes"
estadd local clust "Country (N=9)"

/*==============================================================================
   SECTION 7: ROBUSTNESS CHECKS
==============================================================================*/


/*------------------------------------------------------------------------------
   MODEL M4: BANK-ONLY FI (FA30 alone, no FDEPOID or FMOBAB)
   */
*========================================================
* MODEL M4: BANK-ONLY FI
*========================================================
* M4: Banks-only FI interaction
eststo M4: xtreg cons_growth ///
    c.shock##c.fi_bank_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

estadd local cfe "Yes"
estadd local yfe "Yes"
estadd local clust "Country (N=9)"



/*------------------------------------------------------------------------------
   MODEL M5: EXCLUDE CRISIS YEARS (2008, 2009, 2020)
  */
*========================================================
* MODEL M5: EXCLUDE CRISIS YEARS
*========================================================

* M5: Exclude crisis years
eststo M5: xtreg cons_growth ///
    c.shock##c.fi_std ///
    inflation trade_openness ///
    i.year if crisis_year==0, fe cluster(country_id)

estadd local cfe "Yes"
estadd local yfe "Yes"
estadd local clust "Country (N=9)"

/*------------------------------------------------------------------------------
   MODEL M6: LAGGED FI — Pre-existing financial access
  */
*========================================================
* MODEL M6: LAGGED LOG FI
*========================================================

* M6: Lagged log FI interaction
eststo M6: xtreg cons_growth ///
    c.shock##c.ln_fi_lag1 ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

estadd local cfe "Yes"
estadd local yfe "Yes"
estadd local clust "Country (N=9)"

/*------------------------------------------------------------------------------
   MODEL M7: WINSORIZED SHOCK
   */
*========================================================
* MODEL M7: WINSORIZED SHOCK
*========================================================

eststo M7: xtreg cons_growth ///
    c.shock_w##c.fi_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

estadd local cfe   "Yes"
estadd local yfe   "Yes"
estadd local clust "Country (N=9)"


*========================================================
* MODEL M8: OIL SHOCK
*========================================================

eststo M8: xtreg cons_growth ///
    c.shock_oil##c.fi_std ///
    inflation trade_openness, ///
    fe cluster(country_id)

estadd local cfe   "Yes"
estadd local yfe   "No — oil collinear with year FE"
estadd local clust "Country (N=9)"
/*------------------------------------------------------------------------------*/


/*==============================================================================*/
   /*SECTION 8: MARGINAL EFFECTS
   The interaction β₃ alone is not the full story. The total marginal effect
   of the shock on consumption depends on the level of FI:
   ME = ∂(cons_growth)/∂(SHOCK) = β₁ + β₃ × fi_std
   We compute this at five points in the FI distribution.
==============================================================================*/



xtreg cons_growth ///
    c.shock##c.fi_bank_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)
	
margins, dydx(shock) at(fi_bank_std=(-1 0 1 2))
	


marginsplot, ///
    recast(line) recastci(rarea) ciopt(color(%20)) ///
    title("Marginal Effect of ToT Shock on Consumption Growth") ///
    subtitle("At different levels of Financial Inclusion (standardized)") ///
    xtitle("Financial Inclusion (SD units)") ///
    ytitle("Marginal Effect (pp)") ///
    yline(0, lpattern(dash) lcolor(gray)) ///
    note("Two-way FE. SE clustered by country. 95% CI shaded.")
graph export "output/Figure_MarginalEffects.png", replace width(1400)


*Bank Only:

xtreg cons_growth ///
    c.shock##c.fi_bank_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

margins, dydx(shock) at(fi_bank_std=(-2 -1 0 1 2))


marginsplot, ///
    recast(line) ///
    recastci(rarea) ///
    ciopt(color(%20)) ///
    title("Marginal Effect of ToT Shock on Consumption Growth") ///
    subtitle("At different levels of Bank-Based Financial Inclusion") ///
    xtitle("Bank-based Financial Inclusion (SD units)") ///
    ytitle("Marginal Effect (pp)") ///
    yline(0, lpattern(dash) lcolor(gray)) ///
    note("Two-way FE. SE clustered by country. 95% CI shaded.")

graph export "output/Figure_MarginalEffects_BankOnly.png", ///
    replace width(1400)

di "=== MARGINAL EFFECTS: β₁ + β₃ × fi_std ==="
local b1 = _b[shock]
local b3 = _b[c.shock#c.fi_bank_std]
di "β₁ (shock) = `b1'"
di "β₃ (shock_x_fi_std) = `b3'"
di " "
local b1 = _b[shock]
local b3 = _b[c.shock#c.fi_bank_std]

display "ME at -1 SD = " `b1' + `b3'*(-1)
display "ME at mean  = " `b1' + `b3'*(0)
display "ME at +1 SD = " `b1' + `b3'*(1)
display "ME at +2 SD = " `b1' + `b3'*(2)


/*==============================================================================
   SECTION 9: OUTPUT TABLES
==============================================================================*/

/* TABLE 1: Main results (M1–M6) */

esttab M1 M2 M3 M4 M5 M6 ///
    using "output/Table1_Main.rtf", replace ///
    b(3) se(3) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    label ///
    mtitles("Baseline" "Main" "Log FI" "Banks Only" "No Crisis" "Lagged FI") ///
    title("Financial Inclusion and Consumption Smoothing: Panel Fixed Effects") ///
    keep(shock ///
         fi_std c.shock#c.fi_std ///
         ln_fi c.shock#c.ln_fi ///
         fi_bank_std c.shock#c.fi_bank_std ///
         ln_fi_lag1 c.shock#c.ln_fi_lag1 ///
         inflation trade_openness) ///
    order(shock ///
          fi_std c.shock#c.fi_std ///
          ln_fi c.shock#c.ln_fi ///
          fi_bank_std c.shock#c.fi_bank_std ///
          ln_fi_lag1 c.shock#c.ln_fi_lag1 ///
          inflation trade_openness) ///
    stats(N r2_w cfe yfe clust, ///
          labels("Observations" "Within R-squared" "Country FE" "Year FE" "SE clustered")) ///
    addnotes("Dependent variable: household final consumption per capita growth (%)." ///
             "shock = 100 × Δln(terms-of-trade index)." ///
             "fi_std = standardized financial inclusion composite." ///
             "Banks Only model uses FA30: deposit accounts at commercial banks per 1,000 adults." ///
             "Standard errors clustered by country in parentheses." ///
             "* p<0.10, ** p<0.05, *** p<0.01")


/* TABLE 2: Robustness (M7–M8) */
esttab M7 M4 using "output/Table2_Robustness.rtf", replace ///
    label b(3) se(3) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    title("Table 2: Robustness Checks") ///
    mtitles("Winsorized Shock" "Banks Only") ///
    keep(shock shock_w ///
         fi_std c.shock#c.fi_std ///
         fi_bank_std c.shock#c.fi_bank_std ///
         inflation trade_openness) ///
    order(shock shock_w ///
          fi_std c.shock#c.fi_std ///
          fi_bank_std c.shock#c.fi_bank_std ///
          inflation trade_openness) ///
    stats(N r2_w cfe yfe clust, fmt(0 3) ///
          labels("Observations" ///
                 "Within R-squared" ///
                 "Country FE" ///
                 "Year FE" ///
                 "SE clustered")) ///
    addnotes("Dependent variable: household final consumption per capita growth (%)." ///
             "shock = 100 × Δln(terms-of-trade index)." ///
             "M7: shock winsorized at the 1st and 99th percentiles to test sensitivity to extreme shocks." ///
             "M8: bank-based financial inclusion uses FA30 only (commercial bank deposit accounts per 1,000 adults)." ///
             "Standard errors clustered by country in parentheses." ///
             "* p<0.10, ** p<0.05, *** p<0.01") ///
    nogaps compress

/* COEFFICIENT PLOT */
coefplot ///
    (M2, label("Main specification") ///
         msymbol(O) mcolor(navy) ciopts(color(navy%25))) ///
    (M5, label("Excluding crisis years") ///
         msymbol(S) mcolor(maroon) ciopts(color(maroon%25))) ///
    (M7, label("Winsorized ToT shock") ///
         msymbol(D) mcolor(dkgreen) ciopts(color(dkgreen%25))) ///
    (M4, label("Bank-based FI only") ///
         msymbol(T) mcolor(orange) ciopts(color(orange%25))), ///
    keep(c.shock#c.fi_std ///
         c.shock_w#c.fi_std ///
         c.shock#c.fi_bank_std) ///
    rename(c.shock#c.fi_std      = "Shock × FI composite" ///
           c.shock_w#c.fi_std    = "Winsorized shock × FI composite" ///
           c.shock#c.fi_bank_std = "Shock × bank-based FI") ///
    xline(0, lpattern(dash) lcolor(gray)) ///
    title("Interaction Coefficients Across Specifications") ///
    subtitle("Negative coefficients imply stronger consumption smoothing") ///
    note("90% confidence intervals shown. Two-way fixed effects. Standard errors clustered by country.") ///
    legend(pos(6) rows(2)) ///
    levels(90)

graph export "output/Figure_CoefficientPlot.png", ///
    replace width(1400)

/*==============================================================================
  */

*========================================================
* Wild cluster bootstrap: main interaction
*========================================================

* With 9 clusters, boottest enumerates all 512 Rademacher sign combinations.

xtreg cons_growth ///
    c.shock##c.fi_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

boottest c.shock#c.fi_std, ///
    cluster(country_id) reps(9999) ///
    boottype(wild) seed(42) noci

di "Wild bootstrap p-value, main interaction: " r(p)




*========================================================
* Wild cluster bootstrap: bank-only interaction
*========================================================

xtreg cons_growth ///
    c.shock##c.fi_bank_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

boottest c.shock#c.fi_bank_std, ///
    cluster(country_id) reps(9999) ///
    boottype(wild) seed(42) noci

di "Wild bootstrap p-value, bank-only interaction: " r(p)


*========================================================
* Lead-shock placebo test
*========================================================

gen shock_f1 = F.shock
gen shock_f2 = F2.shock

xtreg cons_growth ///
    c.shock##c.fi_std ///
    c.shock_f1##c.fi_std ///
    c.shock_f2##c.fi_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

test c.shock_f1#c.fi_std c.shock_f2#c.fi_std

di "Lead-shock placebo test p-value: " r(p)
di "Interpretation: p > 0.10 suggests no evidence that future shocks predict current consumption through FI."


*========================================================
* Final summary to log
*========================================================

di "========================================================"
di "FINAL ESTIMATION SUMMARY"
di "========================================================"

xtreg cons_growth ///
    c.shock##c.fi_std ///
    inflation trade_openness ///
    i.year, fe cluster(country_id)

di "N observations = " e(N)
di "N countries    = " e(N_g)
di "Within R2      = " e(r2_w)
di "Beta shock     = " _b[shock] "   SE = " _se[shock]
di "Beta FI        = " _b[fi_std] "   SE = " _se[fi_std]
di "Beta interaction = " _b[c.shock#c.fi_std] "   SE = " _se[c.shock#c.fi_std]

di " "
di "KEY INTERPRETIVE NOTE:"
di "The interaction coefficient tests whether financial inclusion changes the pass-through of ToT shocks to consumption growth."
di "Negative interaction = higher FI is associated with weaker shock transmission."

log close

