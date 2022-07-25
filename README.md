# bbsStanBayes

early attempt to translate bbsBayes models from JAGS to Stan (<https://github.com/BrandonEdwards/bbsBayes>)

Initially working through the GAMYE model from bbsBayes package.

Eventually planning to add a spatial prior to the information on trajectory (GAM coefficients) and abundance (intercepts), using a conditional autoregressive, neighbourhood structure. Also improving the GAM basis function, using the mgcv package to generate the basis. Mixing is much improved.
