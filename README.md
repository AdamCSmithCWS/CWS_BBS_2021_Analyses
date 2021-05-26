# bbsStanBayes
early attempt to translate bbsBayes models to Stan

Initially working through the GAMYE model from bbsBayes.
Adding a spatial prior to the information on trajectory (GAM coefficients) and abundance (intercepts), using a conditional autoregressive, neighbourhood structure. 
Also improving the GAM basis function, using the mgcv package to generate the basis. Mixing is much improved.

