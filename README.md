# Financial Inclusion and Consumption Smoothing under External Income Shocks

Evidence from a nine-country panel, 2004-2024. Malik Jibran Tariq, Columbia SIPA, May 2026.

## Question
Does formal financial inclusion reduce the pass-through of terms-of-trade shocks to household consumption growth?

## Data
- World Bank WDI: consumption per capita growth, terms of trade, inflation, trade openness
- IMF Financial Access Survey: deposit accounts (FA30), credit union depositors (FDEPOID), mobile money accounts (FMOBAB)
- Countries: Bangladesh, Egypt, Ghana, India, Kenya, Pakistan, Philippines, Tanzania, Uganda

## Method
- Two-way fixed effects with a Shock x FI interaction
- Standard errors clustered by country
- Wild bootstrap inference (boottest, 512 Rademacher permutations) for nine clusters
- Robustness checks: log FI, bank-only FI, crisis-year exclusion, lagged FI, winsorized shocks, placebo leads

## Key findings
- The interaction coefficient is negative in all specifications
- The lagged FI interaction is significant at 10%
- The bank-only effect is larger than the composite, consistent with a savings plus credit channel

## Repository
- `data/`   analysis panel (panel_FINAL.csv)
- `code/`   annotated Stata do-file
- `output/` created when the code runs (tables, figures, log)
- `paper/`  PDF

## Replication
From the repository root, run `do code/main.do` in Stata. Uncomment the `ssc install` lines on first run

## Contact
mmt2197@columbia.edu | linkedin.com/in/malik-jibran-tariq
