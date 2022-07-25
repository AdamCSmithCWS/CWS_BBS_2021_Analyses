
data {
int<lower=1> N;
  vector[N] x;             
  vector[N] y;               
}

parameters {
  real a; 
  real b;
  real<lower=0> sigma;    
}

model {
  sigma ~ student_t(3,0,1); 
  b ~ std_normal();
  a ~ std_normal();

 y ~ normal(a+b*x,sigma); 

}

generated quantities {

   vector[N] log_lik; 

  for(i in 1:N){
   log_lik[i] = normal_lpdf(y[i] | a+b*x[i], sigma);
   }
 
  }

