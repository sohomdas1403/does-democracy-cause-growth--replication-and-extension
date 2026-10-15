# Democracy Does Cause Growth Extension & Replication Project


An extension of Acemoglu, Naidu, Restrepo & Robinson (2019), *"Democracy Does Cause Growth"* (Journal of Political Economy). The project tests whether the paper's instrument, regional waves of democratisation, satisfies the exclusion restriction. It does this by adding lagged regional GDP as a control in both stages of the 2SLS estimation. Includes global and regional estimates (LAC, EU, Africa, Asia, MENA, etc.) with first-stage diagnostics and long-run effects. 

## Motivation

The original paper instruments a country's democracy status with the share of democracies among other countries in the same region and with the same initial regime. The exclusion restriction requires that, conditional on country fixed effects, year fixed effects and lags of GDP, the instrument affects a country's GDP only through its own democracy status. This can fail if regional economic conditions drive both regional democratisation and a country's own growth.

## Model

We add $p$ lags of regional GDP, $y^R$, to both stages.

**First stage**

$$D_{ct} = \sum_{j=1}^{q} \pi_j Z_{c,t-j} + \sum_{j=1}^{p} \phi_j y_{c,t-j} + \sum_{j=1}^{p} \kappa^d_j y^R_{c,t-j} + \theta_c + \eta_t + v_{ct}$$

**Second stage**

$$y_{ct} = \beta D_{ct} + \sum_{j=1}^{p} \gamma_j y_{c,t-j} + \sum_{j=1}^{p} \kappa^y_j y^R_{c,t-j} + \alpha_c + \delta_t + \varepsilon_{ct}$$

Specifications use 4 lags of GDP, 4 lags of regional GDP and 4 lags of the democracy wave instrument. The **Extension_Regression_Results.pdf** file contains all produced results and regression tables, including global and regional estimates.

## Data & How to run

- Country-year panel, 1960–2010 (174 countries in the main sample), from the original paper's replication files
- The original paper's replication files are titled "replication_files_ddcg". Set base directory of the .do file to this folder.
- Run the stata code
- Tables are saved as .csv files to 'results/'

## Reference

Acemoglu, D., Naidu, S., Restrepo, P., & Robinson, J. A. (2019). Democracy Does Cause Growth. *Journal of Political Economy*, 127(1), 47–100.

This extension was part of a larger group project for *Foundations of Econometrics and Machine Learning*, with Raffier J., Sudan M., Sindler D., Mena G., Das S. and Ducotterd M.
