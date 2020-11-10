// This is a model for the slopes nad intercepts of route-specific BBS trends based on covariates
// plus a full intrinsic CAR component that estimates mean slopes per stratum 
// 
// while building this, consider running a CAR route-specific slope model to estimate route-specific trends of species
// then run subsequent model that explains trends with route-specific covariates - that is separate this model into two steps


functions {
  real icar_normal_lpdf(vector bb, int nstrata, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nstrata); //soft sum to zero constraint on phi
 }
}


data {
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nroutes;
  int<lower=1> nyears;
  int<lower=1> nspecies;
  int<lower=1> nobservers;
  // data for spline s(year)

 
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];              // strata indicators
  int<lower=1> species[ncounts];              // species indicators
 int<lower=1> observer[ncounts];              // observer indicators
  int<lower=1> route[ncounts];              // route indicators
  real year[ncounts];              // centered years
  int<lower=1> year_raw[ncounts]; // year index
  int<lower=1> date[ncounts];  // day indicator in the season
  
  int<lower=1> strat_rt[nroutes]; //vector of strata assignements for each route
  
  // Predictors
  real cov_i[nroutes]; //route-level intercept predictor values (recent human footprint values)
  real cov_s[nroutes]; //route-level slope predictor values (human footprint change values)
  real cov_sw[nspecies]; // species wintering region slope predictor (human footprint change values)
  
  //indexes for re-scaling predicted counts within strata based on route-level intercepts
  int<lower=1> max_routes; //dimension 1 of routes matrix
  int<lower=0> routes[max_routes,nstrata]; //matrix of which routes are in each stratum
  int<lower=1> nroutes_strat[nstrata]; //number of unique routes in each stratum
  
  int species_strat[nstrata,nspecies]; //matrix of 1s and 0s - species by strata inclusion matrix
}

parameters {
  real noise[nspecies,ncounts];             // over-dispersion
 

  real<lower=0> sigma_alpha[nspecies];    // spatial standard deviation of intercepts
  real<lower=0> sigma_beta[nspecies];    // spatial standard deviation of slopes
  real<lower=0> sdnoise[nspecies];    // sd of over-dispersion
 //real<lower=1> nu; 
  real<lower=0> sdobs;    // sd of observer effects
  
  real ALPHA[nspecies]; // overall intercept for each species - fixed effect
  matrix[nstrata,nspecies] alpha_raw;             // spatial intercepts - non-covariate component at a stratum level
  
  real BETA_f[nspecies]; // overall slope for each species - fixed effect
  matrix[nstrata,nspecies] beta_raw;             // spatial intercepts
  
  real cov_i_eff_raw[nspecies]; //coefficients for intercept effects
  real cov_I_HYPER; //hyperparameter for the effect of covariate on the intercept
  real cov_s_eff_raw[nspecies]; //coefficients for slope effects
  real cov_S_HYPER; //hyperparameter for the effect of covariate on the slope
  
  real cov_sw_eff; //coefficient for slope effect of wintering region change - only one parameter, each species provides a single data-point

  real obs[nobservers]; //observer effects
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  real intercept[nroutes,nspecies]; 
    real slope[nroutes,nspecies]; 
    real cov_i_eff[nspecies];
    real cov_s_eff[nspecies];
    matrix[nstrata,nspecies] alpha;
    matrix[nstrata,nspecies] beta;
    real BETA[nspecies]; //combination of fixed species trend and wintering region footprint change effects
    
    
// covariate effect on intercepts and slopes
for(sp in 1:nspecies){
cov_i_eff[sp] = cov_I_HYPER + cov_i_eff_raw[sp];

cov_s_eff[sp] = cov_S_HYPER + cov_s_eff_raw[sp];

BETA[sp] = BETA_f[sp] + cov_sw[sp]*cov_sw_eff; // species overall mean trend 


 // spatial components of trends and intercepts   
  
     alpha[,sp] = (sigma_alpha[sp] * alpha_raw[,sp]) + ALPHA[sp];
     
     beta[,sp] = (sigma_beta[sp] * beta_raw[,sp]) + BETA[sp];
    
  }
  
  
//intercepts
// strata intercept by sp, covariate effect, species mean effect, 
// OBSERVER-ROUTE[obs_rt,sp] effect is also needed?
// observer effect would be better. assume it's the same across species...
// x = a ? b : c; conditional operator - if a != 0 then do b, if a == 0 then do c.
for(sp in 1:nspecies){
  for(rt in 1:nroutes){
  intercept[rt,sp] = species_strat[strat_rt[rt],sp] ? alpha[strat_rt[rt],sp] + cov_i[rt]*cov_i_eff[sp] : 0 ; //conditional 0-intercept if species is not present

}
}
//slopes model
for(sp in 1:nspecies){
  for(rt in 1:nroutes){
  slope[rt,sp] = species_strat[strat_rt[rt],sp] ? beta[strat_rt[rt],sp] + cov_s[rt]*cov_s_eff[sp] : 0 ; //conditional 0-intercept if species is not present

}
} 

  for(i in 1:ncounts){
    E[i] =  slope[route[i],species[i]] * year[i] + intercept[route[i],species[i]] + obs[observer[i]] + noise[species[i],i];
  }
  
  }
model {
  for(sp in 1:nspecies){
  sdnoise[sp] ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  noise[sp,] ~ normal(0,sdnoise[sp]); //heavy tailed extra Poisson log-normal variance
  }
  sdobs ~ normal(0,1); //prior on sd of gam hyperparameters
  
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed route-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  obs ~ normal(0,sdobs);//observer effects - same across all species - no route random effect beyond the covariates
  
  count ~ poisson_log(E); //vectorized count likelihood
  
  BETA_f ~ normal(0,0.1);// prior on species fixed effect mean slopes
  ALPHA ~ normal(0,1);// prior on species fixed effect mean intercepts
  
  
  //spatial iCAR intercepts and slopes by strata
  sigma_alpha ~ normal(0,1); //prior on scale of spatial variation
  sigma_beta ~ normal(0,1); //prior on scale of spatial variation

  for(sp in 1:nspecies){
  alpha_raw[,sp] ~ icar_normal_lpdf(nstrata, node1, node2);
  beta_raw[,sp] ~ icar_normal_lpdf(nstrata, node1, node2);
  }
}

// generated quantities {
//   real<lower=0> N[nyears];
//   real<lower=0> NSmooth[nyears];
//   real<lower=0> n[nstrata,nyears];
//   real<lower=0> nsmooth[nstrata,nyears];
//   
//   
//       for(s in 1:nstrata){
// 
//   for(y in 1:nyears){
// 
//       n[s,y] = exp(ALPHA1 + year_pred[y,s] + year_effect[y] + season_pred[50] ) + a[s];
//       nsmooth[s,y] = exp(ALPHA1 + year_pred[y,s] + season_pred[50] ) + a[s];
//     }
//   }
//   
//     for(y in 1:nyears){
// 
//       N[y] = exp(ALPHA1 + Y_pred[y] + year_effect[y] + season_pred[50] );
//       NSmooth[y] = exp(ALPHA1 + Y_pred[y] + season_pred[50] );
//       
//     }
//     
// }

