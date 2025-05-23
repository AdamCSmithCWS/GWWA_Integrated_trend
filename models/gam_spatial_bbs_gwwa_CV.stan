// This is a Stan implementation of the gamye model

// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues and to reduce the model output size (greatly)
// althought nice to have them here where Rhat and ess_ are calculated

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
  
  ///BBS only data
  int<lower=1> nsites;
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  
  array[ncounts] int<lower=0> count;              // count observations
  array[ncounts] int<lower=1> strat;               // strata indicators
  array[ncounts] int<lower=1> year; // year index
  array[ncounts] int<lower=1> site; // site index
  array[ncounts] int<lower=0> firstyr; // first year index
  array[ncounts] int<lower=1> observer;              // observer indicators
  
  int<lower=1> nobservers;// number of observers

  // array data to estimate annual indices using only observer-site combinations that are in each stratum
  array[nstrata] int<lower=0> nobs_sites_strata; // number of observer-site combinations in each stratum
  int<lower=0> maxnobs_sites_strata; //largest value of nobs_sites_strata 
  array[nstrata,maxnobs_sites_strata] int ste_mat; //matrix identifying which sites are in each stratum
  array[nstrata,maxnobs_sites_strata] int obs_mat; //matrix identifying which sites are in each stratum
  // above are effectively ragged arrays, but filled with 0 values so that Stan will accept it as data
  // but throws an error if an incorrect strata-site combination is called
  
  array[nstrata] real nonzeroweight; //proportion of the sites in the stratum included - scaling factor
 

  // data for spline s(year)
  int<lower=1> nknots_year;  // number of knots in the basis function for year
  matrix[nyears, nknots_year] year_basis; // basis function matrix

  //data for spatial iCAR among strata
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nstrata> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nstrata> node2;  // and node1[i] < node2[i]


  // Extra Poisson variance options
  int<lower=0,upper=1> heavy_tailed; //indicator if extra poisson variance should be t-distributed or normal (yes = 1, no = 0 and therefore normal)
  int<lower=0,upper=1> calc_nu; //indicator if nu should be calculated (yes = 1, no = 0)
  int<lower=0,upper=1> use_pois; //indicator if count variation should be based on over-dispersed Poisson (if ==1) or Negative binomial (if == 0)
  
  // loo or CV calculations
  int<lower=0,upper=1> calc_log_lik; //indicator if log_lik should be calculated (log_lik for all obs to support loo = 1, no log-lik = 0)
  int<lower=0,upper=1> calc_CV; //indicator if CV should be calculated (CrossValidation = 1, no CV = 0)
  // CV folds - if calc_CV == 1 then the following values define the training and test sets
  int<lower=1, upper=ncounts> ntrain; //
  int<lower=1, upper=ncounts> ntest; //
  array[ntrain] int<lower=1, upper=ncounts> train; // indices of counts to include in train data
  array[ntest] int<lower=1, upper=ncounts> test; // indices of counts to include in test data
  
  // This approach to CV only works if all observers and routes are included in each training-fold
  // So, CV folds must be nested within observers and routes
  // Could implement leave future observation style CV within observers and routes if indexing was done carefully


// GWWA only data
  int<lower=1> nsites_gwwa;
  int<lower=1> ncounts_gwwa;
  int<lower=1> nobservers_gwwa;// number of observers
 
  array[ncounts_gwwa] int<lower=0> count_gwwa;              // count observations
  array[ncounts_gwwa] int<lower=1> strat_gwwa;               // strata indicators
  array[ncounts_gwwa] int<lower=1> year_gwwa; // year index
  array[ncounts_gwwa] int<lower=1> site_gwwa; // site index
  array[ncounts_gwwa] int<lower=1> observer_gwwa;              // observer indicators
  array[ncounts_gwwa] real off_set; // log(ncounts) - only applies to gwwa survey, off_set[i] == 0 if survey[i] == 1
  




}

transformed data {
   
     array[ntrain] int count_tr = count[train];
     array[ntrain] int strat_tr = strat[train];
     array[ntrain] int year_tr = year[train];
     array[ntrain] int site_tr = site[train];
     array[ntrain] int firstyr_tr = firstyr[train];
     array[ntrain] int observer_tr = observer[train];
     
     array[ntest] int count_te = count[test];
     array[ntest] int strat_te = strat[test];
     array[ntest] int year_te = year[test];
     array[ntest] int site_te = site[test];
     array[ntest] int firstyr_te = firstyr[test];
     array[ntest] int observer_te = observer[test];
     
  
  
  
}


parameters {
  vector[ntrain*use_pois] noise_raw;             // over-dispersion if use_pois == 1
 
 vector[nstrata] strata_raw;
   real STRATA; 
   real STE_gwwa; 

  real eta; //first-year intercept
  

  vector[nobservers] obs_raw;    // sd of year effects
  vector[nsites] ste_raw;   // 
  real<lower=0> sdnoise;    // sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdste;    // sd of site effects
  //array[nknots_year] real<lower=0> sdbeta;    // sd of GAM coefficients among strata 
  real<lower=0> sdbeta;    // sd of GAM coefficients among strata - one value across all k and strata
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdBETA;    // sd of GAM coefficients
  real<lower=3> nu; // df of t-distribution > 3 so that it has a finite mean, variance, kurtosis
  
  vector[nknots_year] BETA_raw;//_raw; 
  matrix[nstrata,nknots_year] beta_raw;         // GAM strata level parameters

//GWWA parameters
  vector[ncounts_gwwa*use_pois] noise_gwwa_raw;             // over-dispersion if use_pois == 1
  real<lower=3> nu_gwwa; // df of t-distribution > 3 so that it has a finite mean, variance, kurtosis
  vector[nobservers_gwwa] obs_gwwa_raw;    // sd of year effects
  vector[nsites_gwwa] ste_gwwa_raw;   // 
  real<lower=0> sdnoise_gwwa;    // sd of over-dispersion
  real<lower=0> sdobs_gwwa;    // sd of observer effects
  real<lower=0> sdste_gwwa;    // sd of site effects

}

transformed parameters { 
  vector[ntrain] E;           // log_scale additive likelihood
  matrix[nstrata,nknots_year] beta;         // spatial effect slopes (0-centered deviation from continental mean slope B)
  matrix[nyears,nstrata] smooth_pred;
  vector[nyears] SMOOTH_pred; 
  vector[nknots_year] BETA;
  real<lower=0> phi; //transformed sdnoise

  vector[ncounts_gwwa] E_gwwa;           // log_scale additive likelihood
  real<lower=0> phi_gwwa; //transformed sdnoise
  
  if(use_pois){
    phi = 0;
    phi_gwwa = 0;
  }else{
    phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    phi_gwwa = 1/sqrt(sdnoise_gwwa); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  }
  
  
  BETA = sdBETA*BETA_raw;
  
  for(k in 1:nknots_year){
    beta[,k] = (sdbeta * beta_raw[,k]) + BETA[k];
  }
  SMOOTH_pred = year_basis * BETA; 
  
  // for(s in 1:nstrata){
  //   beta[s,] = (sdbeta[s] * beta_raw[s,]) + transpose(BETA);
  // } 
  
  for(s in 1:nstrata){
     smooth_pred[,s] = year_basis * transpose(beta[s,]);
}

// intercepts and slopes

  


  for(i in 1:ntrain){
    real noise;
    real obs = sdobs*obs_raw[observer_tr[i]];
    real strata = (sdstrata*strata_raw[strat_tr[i]]) + STRATA;
    real ste = sdste*ste_raw[site_tr[i]]; // site intercepts
    if(use_pois){
    noise = sdnoise*noise_raw[i];
    }else{
    noise = 0;
    }
    
    E[i] =  smooth_pred[year_tr[i],strat_tr[i]] + strata + eta*firstyr_tr[i] + ste + obs + noise;
  }
  
  for(i in 1:ncounts_gwwa){
    real noise;
    real obs = sdobs_gwwa*obs_gwwa_raw[observer_gwwa[i]];
    real strata = (sdstrata*strata_raw[strat_gwwa[i]]) + STRATA;
    real ste = sdste_gwwa*ste_gwwa_raw[site_gwwa[i]] + STE_gwwa; // site intercepts
    if(use_pois){
    noise = sdnoise_gwwa*noise_gwwa_raw[i];
    }else{
    noise = 0;
    }
    
    E_gwwa[i] =  smooth_pred[year_gwwa[i],strat_gwwa[i]] + off_set[i] + ste + obs + noise;
  }
  
  
  }
  

model {
  nu ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
  nu_gwwa ~ gamma(2,0.1); // prior on df for t-distribution of heavy tailed site-effects from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#prior-for-degrees-of-freedom-in-students-t-distribution
 
  if(use_pois){
    if(heavy_tailed){
    if(calc_nu){
      noise_raw ~ student_t(nu,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
      noise_gwwa_raw ~ student_t(nu_gwwa,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
   }else{
      noise_raw ~ student_t(3,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
      noise_gwwa_raw ~ student_t(3,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
    }
   }else{
    noise_raw ~ std_normal();//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
    noise_gwwa_raw ~ std_normal();//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
    }
  }
  if(use_pois){
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  sdnoise_gwwa ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  }else{
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  sdnoise_gwwa ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance or inverse sqrt(phi) for negative binomial
  }  
  sdobs ~ normal(0,0.3); // informative prior on scale of observer effects - suggests observer variation larger than 3-4-fold differences is unlikely
  sdste ~ student_t(3,0,1); //prior on sd of site effects
  sdBETA ~ student_t(3,0,1); // prior on sd of GAM parameters
  sdbeta ~ student_t(3,0,1); // prior on sd of GAM parameters
  sdstrata ~ student_t(3,0,1); //prior on sd of intercept variation
 
 
  obs_raw ~ std_normal(); // ~ student_t(3,0,1);//observer effects
  //sum(obs_raw) ~ normal(0,0.001*nobservers); // constraint isn't useful here

  ste_raw ~ std_normal();//site effects
  //sum(ste_raw) ~ normal(0,0.001*nsites); //constraint isn't useful here

//GWWA 
  sdobs_gwwa ~ normal(0,0.3); // informative prior on scale of observer effects - suggests observer variation larger than 3-4-fold differences is unlikely
  sdste_gwwa ~ student_t(3,0,1); //prior on sd of site effects
  obs_gwwa_raw ~ std_normal(); // ~ student_t(3,0,1);//observer effects
  ste_gwwa_raw ~ std_normal();//site effects
  sum(ste_gwwa_raw) ~ normal(0,0.001*nsites_gwwa); //constraint isn't useful here
  STE_gwwa ~ normal(0,1);// prior on fixed effect mean intercept




  BETA_raw ~ std_normal();// prior on fixed effect mean GAM parameters
  //sum to zero constraint built into the basis function

  
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
 
for(k in 1:nknots_year){
    beta_raw[,k] ~ icar_normal(nstrata, node1, node2);;
}
   strata_raw ~ icar_normal(nstrata, node1, node2);
    //sum(strata_raw) ~ normal(0,0.001*nstrata);
 
if(use_pois){
  count_tr ~ poisson_log(E); //vectorized count likelihood with log-transformation
  count_gwwa ~ poisson_log(E_gwwa); //vectorized count likelihood with log-transformation
}else{
   count_tr ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
   count_gwwa ~ neg_binomial_2_log(E_gwwa,phi_gwwa); //vectorized count likelihood with log-transformation
 
}

}

 generated quantities {
// only apply to BBS
   array[nstrata,nyears] real<lower=0> n;
   real<lower=0> retrans_noise;
   real<lower=0> retrans_obs;
   real<lower=0> retrans_ste;
   vector[ncounts*calc_log_lik] log_lik; // alternative value to track the observervation level log-likelihood
   vector[ntest*calc_CV] log_lik_cv; // alternative value to track the log-likelihood of the coutns in the test dataset
   real adj;
 
  if(calc_log_lik){
  // potentially useful for estimating loo-diagnostics, such as looic
  if(use_pois){
  for(i in 1:ncounts){
   log_lik[i] = poisson_log_lpmf(count_tr[i] | E[i]);
   }
  }else{
   for(i in 1:ncounts){
   log_lik[i] = neg_binomial_2_log_lpmf(count_tr[i] | E[i] , phi);
   } 
  }
  }
  
  if(calc_CV){
    for(i in 1:ntest){
      
    real noise;
    real obs = sdobs*obs_raw[observer_te[i]];
    real strata = (sdstrata*strata_raw[strat_te[i]]) + STRATA;
    real ste = sdste*ste_raw[site_te[i]]; // site intercepts
   
   if(use_pois){
      if(heavy_tailed){
        if(calc_nu){
    noise = student_t_rng(nu,0,sdnoise);
        }else{
    noise = student_t_rng(3,0,sdnoise);
        }
      }else{
    noise = normal_rng(0,sdnoise);
      }
   
      
   log_lik_cv[i] = poisson_log_lpmf(count_te[i] | smooth_pred[year_te[i],strat_te[i]] + strata + eta*firstyr_te[i] + ste + obs + noise);
  
   }else{
     noise = 0;
  log_lik_cv[i] = neg_binomial_2_log_lpmf(count_te[i] | smooth_pred[year_te[i],strat_te[i]] + strata + eta*firstyr_te[i] + ste + obs + noise, phi);
  
   }
  
  }
  
  }
  
  if(use_pois){
  if(heavy_tailed){
      if(calc_nu){
         adj = (1.422*(nu^0.906))/(1+(1.422*(nu^0.906)));
        }else{
         adj = (1.422*(3^0.906))/(1+(1.422*(3^0.906)));
        }
    }else{
      adj = 1;
    }
     retrans_noise = 0.5*(sdnoise/adj)^2;
}else{
  adj = 1;
  retrans_noise = 0;
}
     
retrans_obs = 0.5*(sdobs^2);
retrans_ste = 0.5*(sdste^2);

// Annual indices of abundance - strata-level annual predicted counts


for(y in 1:nyears){

      for(s in 1:nstrata){

  array[nobs_sites_strata[s]] real n_t;
  real strata = (sdstrata*strata_raw[s]) + STRATA;
  
        for(t in 1:nobs_sites_strata[s]){  //nobs_sites_strata maxnobs_sites_strata

  real ste = sdste*ste_raw[ste_mat[s,t]]; // site intercepts
  real obs = sdobs*obs_raw[obs_mat[s,t]]; // site intercepts


      n_t[t] = exp(strata+ smooth_pred[y,s] + retrans_noise + obs + ste);// + retrans_obs);
        }
        n[s,y] = nonzeroweight[s] * mean(n_t);//mean of exponentiated predictions across sites in a stratum
        //using the mean of hte exponentiated values, instead of including the log-normal
        // retransformation factor (0.5*sdste^2), because this retransformation makes 2 questionable assumptions:
          // 1 - assumes that sdste is equal among all strata
          // 2 - assumes that the distribution of site-effects is normal
        // As a result, these annual indices reflect predictions of mean annual abundance within strata of the routes that are included in the stratum


    }
  }



 }

