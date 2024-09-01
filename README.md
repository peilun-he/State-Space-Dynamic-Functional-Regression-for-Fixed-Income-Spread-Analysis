# State-Space Dynamic Functional Regression for Multicurve

This repository contains the code for the paper titled "State-Space Dynamic
Functional Regression for Multicurve Fixed Income Spread Analysis and Stress Testing".

For test purpose, one can run the *test.m* file in Matlab. This file will apply the DNS
model and the DNS-FR model on UK bond yields and calculate the RMSE. The following table
gives the expected results after running this file. 

| Maturity | DNS model | DNS-FR model |
|----------|-----------|--------------|
| 1 year   | 0.0379    | 0.0802       |
| 2 years  | 0.0236    | 0.0121       |
| 3 years  | 0.0182    | 0.0083       |
| 5 years  | 0.0251    | 0.0110       |
| 7 years  | 0.0220    | 0.0109       |
| 10 years | 0.0860    | 0.0361       |
| 20 years | 0.2780    | 0.2272       |

The entire code used in the paper are stored in the *main.m* and *NSFR.R* files. 

## Data

1. US_yields.csv: US Treasury yields from 2005-01-31 to 2010-12-31. Sourced from
   the "YieldCurve" package in R.
2. UK_yields.csv: UK gilt yields from 2005-01-31 to 2010-12-31. The monthly government
   liability nominal yield curves are used. Downloaded from
   [Bank of England](https://www.bankofengland.co.uk/statistics/yield-curves)
   official website.
3. US_factors.csv: Extracted factors using US yields. This is the exactly same data
   after running *test.R* file. 

