# Peru Inflation and Macroeconomic Stabilization (1980–2022)

A regime-aware macro analytics case study using World Bank data.

## Project question
How sensitive are simple inflation models to crisis years, regime shifts, and specification choices?

## Main findings
- Peru’s inflation series is dominated by the 1988–1991 hyperinflation episode.
- Full-sample OLS is highly sensitive to a few influential crisis years, especially 1990.
- Broad money growth and FX depreciation remain positively associated with inflation, but coefficient magnitudes shrink materially outside the crisis regime.

## Methods
- World Bank API via `wbstats`
- Feature engineering for FX depreciation, terms-of-trade growth, rolling volatility, and regime indicators
- OLS, robust SEs, Newey-West SEs
- Influence diagnostics
- Specification comparisons across full, restricted, and split samples

## Files
- `run.R`: full analysis pipeline
- `peru_macro_stability_report.Rmd`: report source
- `peru_macro_stability_report.pdf`: final report
