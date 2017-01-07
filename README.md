# SNsims
[![Build Status](https://travis-ci.org/rbiswas4/SNsims.svg?branch=master)](https://travis-ci.org/rbiswas4/SNsims.svg?branch=master) [![Coverage Status](https://coveralls.io/repos/github/rbiswas4/SNsims/badge.svg)](https://coveralls.io/github/rbiswas4/SNsims)

A set of classes to provide structure to SN simulations using SNCosmo or LSST sims packages.

# Requirements 
- `lsst_sims`: package repository: This is a set of packages that provide the bulk of the calculations necessary for computing fluxes and uncertainties of a SN given the parameters of an allowed model. This can sometimes be achieved by using a number of third party packages. Some of these examples will be listed below.
- `sncosmo` : This provides methods to calculate the spectra for models with a spectral time series. These spectra are then used to calculate band fluxes using `lsst.sims.photUtils`.
- `Healpy` : This is made available through the `lsst_sims` package but can be installed externally as well. 

# Data Products
- While not necessary in principle, the use of precomputed sets of pointings per tile are useful to save time. Such precomputed pointings are made available as downloadable sqlite databases. For the current baseline `minion_1016`, the recommended data product is currently [available](http://astro.washington.edu/sim_data/hpOpSim/healpixelized_MINION_1016_256_64_indexed_ddfnodith.db.gz') as a link.

# Contributors
- R. Biswas
- D. Rubin
