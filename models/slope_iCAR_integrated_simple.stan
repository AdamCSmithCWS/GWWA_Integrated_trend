// This is a Stan implementation of a site-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only


// iCAR function, from Morris et al. 2019
// Morris, M., K. Wheeler-Martin, D. Simpson, S. J. Mooney, A. Gelman, and C. DiMaggio (2019). 
// Bayesian hierarchical spatial models: Implementing the Besag York Mollié model in stan. 
// Spatial and Spatio-temporal Epidemiology 31:100301.

 functions {
   real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
     return -0.5 * dot_self(bb[n1] - bb[n2])
       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
  }
 }



data {
  // dimension definitions
  int<lower=1> nsites;
  int<lower=1> ncounts;
  int<lower=1> ncounts_bbs;//only for setting length of noise_raw_bbs vector
  int<lower=1> ncounts_gwwa;//only for setting length of noise_raw_bbs vector
  int<lower=1> nyears;
  int<lower=1> nobservers_bbs;
  int<lower=1> nobservers_gwwa;
  
 
  // count-level data
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=0, upper=ncounts_bbs> inds_bbs; // separate bbs count-level indicators to  track which gwwa observations link to which count. inds_bbs[i] == 0 if survey[i] == 0 (gwwa)
  array [ncounts] int<lower=0, upper=ncounts_gwwa> inds_gwwa;  // separate gwwa count-level indicators to track which gwwa observations link to which count. inds_gwwa[i] == 0 if survey[i] == 1 (bbs)
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> site; // site index
  array [ncounts] int<lower=0, upper=1> survey; //survey index 1 == BBS 0 == gwwa
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators - survey specific and so must be combined with survey[ncounts] to make sense
  array [ncounts] real off_set; // log(ncounts) - only applies to gwwa survey, off_set[i] == 0 if survey[i] == 1

  // site-level data
  array [nsites] int<lower=0, upper=1> survey_sites; //survey index for sites 1 == BBS 0 == gwwa
  array [nsites] int<lower=1> site_bbs; // vector of BBS site indices for each integrated site
  array [nsites] int<lower=1> site_gwwa; // vector of GWWA site indices for each integrated site
  
  // constant
  int<lower=1> fixedyear; // centering value for years
 
  // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nsites> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nsites> node2;  // and node1[i] < node2[i]


}

parameters {
  vector[ncounts_bbs] noise_raw_bbs;             // over-dispersion for BBS
  vector[ncounts_gwwa] noise_raw_gwwa;             // over-dispersion for GWWA

  vector[nsites] beta_raw_space;
  //vector[nsites] beta_raw_rand;
  real BETA; 

  vector[nsites] alpha_raw;  //route/quad-level abundance 
  real ALPHA_bbs; 
  real ALPHA_gwwa; 

  real eta; //first-year intercept
  
  vector[nobservers_bbs] obs_raw_bbs; //observer effects
  vector[nobservers_gwwa] obs_raw_gwwa; //observer effects

  real<lower=0> sdnoise_bbs;    // scale of over-dispersion for BBS counts
  real<lower=0> sdnoise_gwwa;    // scale of over-dispersion for GWWA counts
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs_bbs;    // sd of observer effects BBS
  real<lower=0> sdobs_gwwa;    // sd of observer effects GWWA
  real<lower=0> sdbeta_space;    // sd of slopes in space 
  //real<lower=0> sdbeta_rand;    // sd of slopes random
  real<lower=0> sdalpha;    // sd of intercepts

  
}




transformed parameters{

  vector[ncounts] E;           // log_scale additive likelihood
   //vector[nsites] beta_rand;
  vector[nsites] beta_space;
 vector[nsites] beta;
  vector[nsites] alpha;

// covariate effect on intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   //beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_space + BETA;
   alpha = (sdalpha*alpha_raw);// + ALPHA;
   //obs = sdobs*obs_raw;

// if statement to allow each count to be modeled by either BBS parameters or GWWA parameters
// while integrating information on beta 
// survey == 1 if BBS, survey == 0 if gwwa
  for(i in 1:ncounts){
       if(survey[i]){
   real obs = sdobs_bbs * obs_raw_bbs[observer[i]];
   real noise = sdnoise_bbs*noise_raw_bbs[inds_bbs[i]];
   E[i] =  beta[site[i]] * (year[i]-fixedyear) + ALPHA_bbs +  alpha[site[i]] + off_set[i] + obs + eta*firstyr[i] + noise;
    }
else    {
   real obs = sdobs_gwwa * obs_raw_gwwa[observer[i]];
   real noise = sdnoise_gwwa*noise_raw_gwwa[inds_gwwa[i]];
   E[i] =  beta[site[i]] * (year[i]-fixedyear) + ALPHA_gwwa + alpha[site[i]] + obs + off_set[i] + noise;
    }
  }
  
  
}


model {


  
  // beta_raw_rand ~ normal(0,1);//observer effects
  // sum(beta_raw_rand) ~ normal(0,0.001*nsites);

  
  sdnoise_bbs ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw_bbs ~ normal(0,1); //~ student_t(4,0,1); //normal tailed extra Poisson log-normal variance
  sdnoise_gwwa ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw_gwwa ~ normal(0,1); //~ student_t(4,0,1); //normal tailed extra Poisson log-normal variance
  
  sdobs_bbs ~ std_normal(); //prior on sd of gam hyperparameters
  obs_raw_bbs ~ normal(0,1);//observer effects
  sum(obs_raw_bbs) ~ normal(0,0.001*nobservers_bbs);


  sdobs_gwwa ~ std_normal(); //prior on sd of gam hyperparameters
  obs_raw_gwwa ~ normal(0,1);//observer effects
  sum(obs_raw_gwwa) ~ normal(0,0.001*nobservers_gwwa);

  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA_bbs ~ normal(0,1);// prior on fixed effect mean intercept
  ALPHA_gwwa ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ student_t(4,0,1); //prior on sd of intercept variation
  sdbeta_space ~ gamma(2,20);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.1 and 99% < 0.33
  //sdbeta_rand  ~ gamma(2,20);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation

  beta_raw_space ~ icar_normal(nsites, node1, node2);
  alpha_raw ~ icar_normal(nsites, node1, node2);


}

 generated quantities {

     vector[ncounts] log_lik;
     matrix[nsites,nyears] indices_plot; //site level annual index scaled to a common survey
     matrix[nsites,nyears] indices; //site level annual index for comparison to observed means
     
   
     
for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  
  //calculating site-level annual indices, scaled to GWWA counts to generate a composite trend
for(y in 1:nyears){
 
 for(s in 1:nsites){
    
     indices_plot[s,y] = exp(beta[s] * (y-fixedyear) + ALPHA_gwwa + alpha[s] + sdnoise_gwwa*sdnoise_gwwa*0.5 + sdobs_gwwa*sdobs_gwwa*0.5); 


    if(survey_sites[s]){
      indices[s,y] = exp(beta[s] * (y-fixedyear) + ALPHA_bbs + alpha[s] + sdnoise_bbs*sdnoise_bbs*0.5 + sdobs_bbs*sdobs_bbs*0.5 + log(50)); 
    }else{
     indices[s,y] = exp(beta[s] * (y-fixedyear) + ALPHA_gwwa + alpha[s] + sdnoise_gwwa*sdnoise_gwwa*0.5 + sdobs_gwwa*sdobs_gwwa*0.5 + log(5)); 
    }
 }
 

  }
 
 }
 
