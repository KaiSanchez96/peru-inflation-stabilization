# Peru Inflation Model Audit

This project audits a simple inflation model for Peru using publicly available World Bank data. Rather than treating a high in-sample fit as proof of model quality, the analysis shows how hyperinflation years, influential observations, and regime shifts can distort interpretation.

## Why this project matters

Analysts often trust models with high R² without checking whether a small number of extreme observations drive the result. This project shows why that can be misleading, especially in historical macroeconomic data with structural breaks and crisis episodes.

## Questions

- Does a simple inflation model look strong only because crisis years dominate the sample?
- How sensitive are coefficients to influential observations and regime shifts?
- What changes when we compare alternative specifications instead of relying on one baseline model?

## Dataset

The dataset was built from publicly available World Bank indicators using the `wbstats` package in R. Annual observations were collected for Peru from 1980 to 2022.

Variables include:
- Inflation
- Broad money growth
- Exchange-rate depreciation
- Government consumption
- Terms-of-trade growth

## Methods

- World Bank API data collection
- Data cleaning and feature engineering in R
- Exploratory visual analysis
- Baseline linear model
- Influence diagnostics
- Specification comparisons
- Reproducible reporting in R Markdown

## Main findings

- The full-sample model appears strong, but much of its fit is driven by hyperinflation years.
- A small number of influential observations dominate the regression.
- Estimated relationships shift across periods and specifications.
- High fit alone is not enough to trust a model.

## Skills demonstrated

- Reproducible data collection
- Data wrangling
- Feature engineering
- Diagnostic analysis
- Robustness checks
- Analytical communication

## Repository structure
- `run.R`: full analysis pipeline
- `peru_macro_stability_report.Rmd`: report source
- `peru_macro_stability_report.pdf`: final report
