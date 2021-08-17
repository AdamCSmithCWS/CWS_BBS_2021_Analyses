// This is a Stan implementation of the bbsBayes gamye model
// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues
// with iCAR component for the stratum smooths
// iCAR function

functions {
  real icar_normal_lpdf(vector bb, int ns, int[] n1, int[] n2) {
    return -0.5 * dot_self(bb[n1] - bb[n2])
      + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
 }
}


data {
  int<lower=1> nroutes;
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> route[ncounts]; // route index
 
  int<lower=0> firstyr[ncounts]; // first year index
  
  int<lower=1> observer[ncounts];              // observer indicators
  int<lower=1> nobservers;
//  int<lower=1> sum_observers; //dimension for obs_raw_v vector
//  int<lower=1> max_nobservers; //dimension for obs_raw matrix
//  int<lower=1> obs_mat[nstrata,max_nobservers] ;
  
 int<lower=0> nroutes_strata[nstrata]; // number of routes in each stratum
  int<lower=0> maxnroutes_strata; //largest value of nroutes_strata

  int rte_mat[nstrata,maxnroutes_strata]; //matrix identifying which routes are in each stratum
  // above is actually a ragged array, but filled with 0 values so that it works
  // but throws an error if an incorrect strata-route combination is called
  real nonzeroweight[nstrata]; //proportion of the routes included - scaling factor
 
  // spatial neighbourhood information
  int<lower=1> N_edges;
  int<lower=1, upper=nstrata> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nstrata> node2[N_edges];  // and node1[i] < node2[i]



  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basispred; // basis function matrix
 

}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
 
 vector[nstrata] strata_raw;
   real STRATA; 

  real eta; //first-year intercept
  
  matrix[nstrata,nyears] yeareffect_raw;

  vector[nobservers] obs_raw;    // sd of year effects
  vector[nroutes] rte_raw;   // 
  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdrte;    // sd of route effects
  real<lower=0> sdbeta[nknots_year];    // sd of GAM coefficients among strata 
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdBETA;    // sd of GAM coefficients
  real<lower=0> sdyear[nstrata];    // sd of year effects
 // real<lower=4,upper=500> nu; // df of t-distribution > 4 so that it has a finite mean, variance, kurtosis
  
  vector[nknots_year] BETA_raw;//_raw; 
  matrix[nstrata,nknots_year] beta_raw;         // GAM strata level parameters

}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[nstrata] strata;
  vector[nroutes] rte; // route intercepts
  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  matrix[nyears,nstrata] year_pred;
  vector[nyears] Y_pred;  

  //vector[nobservers] obs; //observer effects
  matrix[nstrata,nyears] yeareffect;
  //vector[ncounts] noise;             // over-dispersion
  vector[nknots_year] BETA;
  
  
  BETA = sdBETA*BETA_raw;
  
  for(k in 1:nknots_year){
    beta[,k] = (sdbeta[k] * beta_raw[,k]) + BETA[k];
  }
  Y_pred = year_basispred * BETA; 
  
      for(s in 1:nstrata){
     year_pred[,s] = year_basispred * transpose(beta[s,]);
}

for(s in 1:nstrata){
    yeareffect[s,] = sdyear[s]*yeareffect_raw[s,];

}

     rte = sdrte*rte_raw;


// intercepts and slopes

  strata = (sdstrata*strata_raw) + STRATA;
  
  

  for(i in 1:ncounts){
    real noise = sdnoise*noise_raw[i];
    real obs = sdobs*obs_raw[observer[i]];

    E[i] =  year_pred[year[i],strat[i]] + strata[strat[i]] + yeareffect[strat[i],year[i]] + rte[route[i]] + obs + eta*firstyr[i] + noise;
  }
  
  }
  
  
  
model {
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ student_t(4,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
   
  sdobs ~ normal(0,0.5); //prior on sd of observer effects
  sdrte ~ std_normal(); //prior on sd of route effects
  sdyear ~ gamma(2,2); // prior on sd of yeareffects - stratum specific, and boundary-avoiding with a prior mode at 0.5 (1/2) - recommended by https://doi.org/10.1007/s11336-013-9328-2 
  sdBETA ~ std_normal(); // prior on sd of GAM parameters
  
  //nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed route-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution

  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  rte_raw ~ std_normal();//route effects
  sum(rte_raw) ~ normal(0,0.001*nroutes);
 
 for(s in 1:nstrata){

  yeareffect_raw[s,] ~ std_normal();
  //soft sum to zero constraint on year effects within a stratum
  sum(yeareffect_raw[s,]) ~ normal(0,0.001*nyears);
  
 }
  
  BETA_raw ~ std_normal();// prior on fixed effect mean GAM parameters
  //sum to zero constraint
  //sum(BETA_raw) ~ normal(0,0.001*nknots_year);
  
  STRATA ~ std_normal();// prior on fixed effect mean intercept
  eta ~ normal(0,0.2);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and gam parameters by strata
  sdstrata ~ std_normal(); //prior on sd of intercept variation
  sdbeta ~ normal(0,1); //prior on sd of GAM parameter variation

for(k in 1:nknots_year){
    beta_raw[,k] ~ icar_normal(nstrata, node1, node2);
}
   strata_raw ~ icar_normal(nstrata, node1, node2);


  count ~ poisson_log(E); //vectorized count likelihood with log-transformation


}

 generated quantities {

   real<lower=0> n[nstrata,nyears];
   real<lower=0> nsmooth[nstrata,nyears];
   real<lower=0> retrans_noise;
   real<lower=0> retrans_obs;
  vector[ncounts] log_lik;
  
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  
  
retrans_noise = 0.5*(sdnoise^2);
retrans_obs = 0.5*(sdobs^2);

for(y in 1:nyears){

      for(s in 1:nstrata){

  real n_t[nroutes_strata[s]];
  real nsmooth_t[nroutes_strata[s]];

        for(t in 1:nroutes_strata[s]){
            real retrans_yr = 0.5*(sdyear[s]^2);

      n_t[t] = exp(strata[s]+ year_pred[y,s] + rte[rte_mat[s,t]] + yeareffect[s,y] + retrans_noise + retrans_obs);
      nsmooth_t[t] = exp(strata[s] + year_pred[y,s] + rte[rte_mat[s,t]] + retrans_yr + retrans_noise + retrans_obs);
        }
        n[s,y] = nonzeroweight[s] * mean(n_t);
        nsmooth[s,y] = nonzeroweight[s] * mean(nsmooth_t);



    }
  }



 }

