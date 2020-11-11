// This is a Stan implementation of the bbsBayes slope model

data {
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
 
  int<lower=0> firstyr[ncounts]; // first year index
  
  int<lower=1> obser[ncounts];              // observer indicators
  int<lower=1> nobservers[nstrata];
  
  real nonzeroweight[nstrata]; //proportion of the routes included - scaling factor
 
 // new data not in jags
 int<lower=1> max_nobservers; //dimension for obs parameter
 
  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basispred; // basis function matrix
 

}

parameters {
  real noise[ncounts];             // over-dispersion
  real lambda[ncounts];             // Poisson means
  
 
  real strata_p[nstrata];
  real STRATA; 

  real eta; //first-year intercept
  
  real yeareffect[nstrata,nyears];

  real obs[nstrata,max_nobservers]; //observer effects

  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta[nknots_year];    // sd of GAM coefficients among strata 
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdyear_gam;    // sd of GAM coefficients
 
  real<lower=0> sdyear[nstrata];    // sd of year effects

  
  vector[nknots_year] BETA; 
  matrix[nstrata,nknots_year] beta_p;         // GAM strata level parameters

}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
 vector[nstrata] strata;
  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
    matrix[nyears,nstrata] year_pred;
  vector[nyears] Y_pred;  


  for(k in 1:nknots_year){
    beta[,k] = (sdbeta[k] * beta_p[,k]) + BETA[k];
  }
  Y_pred = year_basispred * BETA; 
  
      for(s in 1:nstrata){
     year_pred[,s] = year_basispred * transpose(beta[s,]);
}



// intercepts and slopes
for(s in 1:nstrata){
  strata[s] = strata_p[s] + STRATA;
}
  

  for(i in 1:ncounts){
    E[i] =  year_pred[year[i],strat[i]] + strata[strat[i]] + yeareffect[strat[i],year[i]] + obs[strat[i],obser[i]] + eta*firstyr[i] + noise[i];
  }
  
  }
  
  
  
model {

  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  noise ~ normal(0,sdnoise); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ normal(0,1); //prior on sd of gam hyperparameters
  sdyear ~ normal(0,1); // prior on sd of yeareffects - stratum specific
  sdyear_gam ~ normal(0,1); // prior on sd of GAM parameters
  
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed route-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
 for(s in 1:nstrata){
  obs[s,] ~ normal(0,sdobs);//observer effects
    sum(obs[s,]) ~ normal(0,0.001*nobservers[s]);
 
  yeareffect[s,] ~ normal(0,sdyear[s]);
   sum(yeareffect[s,]) ~ normal(0,0.001*nyears);
  
 }
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,sdyear_gam);// prior on fixed effect mean GAM parameters
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdstrata ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of GAM parameter variation

for(k in 1:nknots_year){
  beta_p[,k] ~ normal(0,sdbeta[k]); //prior on stratum-level slopes
    sum(beta_p[,k]) ~ normal(0,0.001*nstrata);
}
  strata_p ~ normal(0,sdstrata); //prior on stratum-level slopes
  
  //sum to zero constraints
  sum(strata_p) ~ normal(0,0.001*nstrata);
  
}

//  generated quantities {
// 
//   real<lower=0> n[nstrata,nyears];
//   real<lower=0> nsmooth[nstrata,nyears];
//   real<lower=0> retrans_noise;
//   
//   retrans_noise = 0.5*(sdnoise^2);
// 
// for(y in 1:nyears){
//   
//       for(s in 1:nstrata){
// 
//   real n_o[nobservers[s]];
//   real nsmooth_o[nobservers[s]];
// 
//         for(o in 1:nobservers[s]){
//       n_o[o] = exp(strata[s]+ year_pred[y,s] + obs[s,o] + yeareffect[s,y] + retrans_noise );
//       nsmooth_o[o] = exp(strata[s] + year_pred[y,s] + obs[s,o] + retrans_noise );
//         }
//         n[s,y] = nonzeroweight[s] * mean(n_o);
//         nsmooth[s,y] = nonzeroweight[s] * mean(nsmooth_o);
//         
//         
//     }
//   }
// 
// 
// 
//  }

