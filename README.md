# High-Dimensional Granger Causality for Climatic Attribution

This GitHub Repo contains all the Rdata (.rds) and the R codes to replicate the results in Friedrich, M., Margaritella, L., Smeekes, S., "High-Dimensional Granger Causality for Climatic Attribution" 2024.

Arxiv version: <https://arxiv.org/abs/2302.03996>

# Data

Two .rds datasets are provided:

I)  AggregateGHG.rds, containing yearly data from 1871 to 2018 on: Temperature Anomaly, Greenhouse gases (aggregate), Solar Activity, Stratospheric Aerosols from Volcanic Activity, Tropospheric Aerosols and Surface Albedo, GDP, El Niño–Southern Oscillation index (ENSO), Ocean Heat Content

II) DisaggregateGHG.rds, containing yearly data from 1871 to 2014 on: Temperature Anomaly, CO2, N2O, CH4, Solar Activity, Stratospheric Aerosols from Volcanic Activity, Tropospheric Aerosols and Surface Albedo, GDP, El Niño–Southern Oscillation index (ENSO), Ocean Heat Content

Data sources and the data pre-processing are described in the main reference paper.

# R Code

Two .R codes are provided: AggregateGHG_CODE.R and DisaggregateGHG_CODE.R based on the two datasets above; both are based on the "HDGCvar" R package, created by two of the authors: <https://github.com/Marga8/HDGCvar> and based on two papers:

[1] A. Hecq, L. Margaritella, S.Smeekes, “Granger Causality Testing in High Dimensional VARs: a Post Double Selection Procedure”, Journal of Financial Econometrics 21 (3), 915-958

[2] A. Hecq, L. Margaritella, S.Smeekes, “Inference in Non-stationary High-Dimensional VARs" arXiv preprint arXiv:2302.01434 (2023).
