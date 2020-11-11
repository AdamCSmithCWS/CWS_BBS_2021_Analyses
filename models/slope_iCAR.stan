// This is a Stan implementation of the bbsBayes slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends


//iCAR function
functions {
  real icar_normal_lpdf(vector bb, int nstrata, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nstrata); //soft sum to zero constraint on phi
 }
}


data {
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> fixedyear; // centering value for years
  
  int<lower=0> firstyr[ncounts]; // first year index
  
  int<lower=1> obser[ncounts];              // observer indicators
  int<lower=1> nobservers[nstrata];
  
  real nonzeroweight[nstrata]; //proportion of the routes included - scaling factor
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]


 // new data not in jags
 int<lower=1> max_nobservers; //dimension for obs parameter
}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
  real lambda[ncounts];             // Poisson means
  
  vector[nstrata] beta_p;
  real BETA; 

  vector[nstrata] strata_p;
  real STRATA; 

  real eta; //first-year intercept
  
  matrix[nstrata,nyears] yeareffect_raw;

  matrix[nstrata,max_nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta;    // sd of slopes 
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdyear[nstrata];    // sd of year effects

  
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[nstrata] beta;
  vector[nstrata] strata;
  matrix[nstrata,nyears] yeareffect;
  matrix[nstrata,max_nobservers] obs;
  vector[ncounts] noise;

// covariate effect on intercepts and slopes
   beta = (sdbeta*beta_p) + BETA;
   strata = (sdstrata*strata_p) + STRATA;
   noise = sdnoise*noise_raw;
   
for(s in 1:nstrata){
    yeareffect[s,] = sdyear[s]*yeareffect_raw[s,];
    obs[s,] = sdobs*obs_raw[s,];
}

  for(i in 1:ncounts){
    E[i] =  beta[strat[i]] * (year[i]-fixedyear) + yeareffect[strat[i],year[i]] + strata[strat[i]] + obs[strat[i],obser[i]] + eta*firstyr[i] + noise[i];
  }
  
  }
  
model {

  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ normal(0,1); //prior on sd of gam hyperparameters
  sdyear ~ normal(0,1); // prior on sd of yeareffects - stratum specific
  
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed route-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
 for(s in 1:nstrata){
  obs_raw[s,] ~ normal(0,1);//observer effects
  sum(obs_raw[s,]) ~ normal(0,0.001*nobservers[s]);
  
  yeareffect_raw[s,] ~ normal(0,1);
  sum(yeareffect_raw[s,]) ~ normal(0,0.001*nyears);
  
 }
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdstrata ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of slope variation

  beta_p ~ icar_normal_lpdf(nstrata, node1, node2);
  strata_p ~ icar_normal_lpdf(nstrata, node1, node2);

  //sum to zero constraints
  sum(strata_p) ~ normal(0,0.001*nstrata);
  sum(beta_p) ~ normal(0,0.001*nstrata);
  
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
//       n_o[o] = exp(strata[s]+ beta[s]*(y-fixedyear) + obs[s,o] + yeareffect[s,y] + retrans_noise );
//       nsmooth_o[o] = exp(strata[s]+ beta[s]*(y-fixedyear) + obs[s,o] + retrans_noise );
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

