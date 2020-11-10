// This is a Stan implementation of the bbsBayes slope model

data {
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nroutes;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> fixedyear; // centering value for years
  
  int<lower=1> firstyr[ncounts]; // first year index
  
  int<lower=1> obser[ncounts];              // observer indicators
  int<lower=1> nobservers[nstrata];
  
  real nonzeroweight[nstrata]; //proportion of the routes included - scaling factor
 
 // new data not in jags
 int<lower=1> max_nobservers; //dimension for obs parameter
}

parameters {
  real noise[ncounts];             // over-dispersion
  real lambda[ncounts];             // Poisson means
  
  real beta_p[nstrata];
  real BETA; 

  real strata_p[nstrata];
  real STRATA; 

  real eta; //first-year intercept
  
  real yeareffect[nstrata,nyears];

  real obs[nstrata,max_nobservers]; //observer effects

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


// covariate effect on intercepts and slopes
for(s in 1:nstrata){
  beta[s] = beta_p[s] + BETA;
  strata[s] = strata_p[s] + STRATA;
}
  

  for(i in 1:ncounts){
    E[i] =  beta[strat[i]] * (year[i]-fixedyear) + strata[strat[i]] + obs[strat[i],obser[i]] + eta*firstyr[i] + noise[i];
  }
  
  }
  
model {

  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  noise ~ normal(0,sdnoise); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ normal(0,1); //prior on sd of gam hyperparameters
  sdyear ~ normal(0,1); // prior on sd of yeareffects - stratum specific
  
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed route-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
 for(s in 1:nstrata){
  obs[s,] ~ normal(0,sdobs);//observer effects
  yeareffect[s,] ~ normal(0,sdyear[s]);
 }
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdstrata ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of slope variation

  beta_p ~ normal(0,sdbeta); //prior on stratum-level slopes
  strata_p ~ normal(0,sdstrata); //prior on stratum-level slopes

}

 generated quantities {

  real<lower=0> n[nstrata,nyears];
  real<lower=0> nsmooth[nstrata,nyears];
  real<lower=0> n_o[nstrata,nyears,max_nobservers];
  real<lower=0> nsmooth_o[nstrata,nyears,max_nobservers];
  real<lower=0> retrans_noise;
  
  retrans_noise = 0.5*(sdnoise^2);

for(y in 1:nyears){
  
      for(s in 1:nstrata){

        for(o in 1:nobservers[s]){
      n_o[s,y,o] = exp(strata[s]+ beta[s]*(y-fixedyear) + obs[s,o] + yeareffect[s,y] + retrans_noise );
      nsmooth_o[s,y,o] = exp(strata[s]+ beta[s]*(y-fixedyear) + obs[s,o] + retrans_noise );
        }
        n[s,y] = nonzeroweight[s] * mean(n_o[s,y,1:nobservers[s]]);
        nsmooth[s,y] = nonzeroweight[s] * mean(nsmooth_o[s,y,1:nobservers[s]]);
        
        
    }
  }



 }

