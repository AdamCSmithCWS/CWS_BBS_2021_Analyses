# bbsStanBayes

early attempt to translate bbsBayes models from JAGS to Stan (<https://github.com/BrandonEdwards/bbsBayes>)

Initially working through the GAMYE and Slope models from bbsBayes package.

Currently the plotting script "plotting_w\_alternate_bbsBayes.R" relies on a development version of bbsBayes at this repo <https://github.com/AdamCSmithCWS/bbsBayes/tree/testing_Stan>, which includes modifications to the `generate_indices()` function and a function on which it depends `extract_index_data()`
