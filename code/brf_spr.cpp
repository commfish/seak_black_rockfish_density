// determine selectivity parameters conditioned on catch at age
// determine current SPR level

#include<TMB.hpp>
template<class Type>

Type objective_function<Type>::operator()()
{
  // data -----
  
  // von b
  DATA_VECTOR(age);
  DATA_VECTOR(length);
  int von_n = age.size();
  
  //selectivity
  DATA_VECTOR(X);
  DATA_VECTOR(prop);
  int n = X.size();
  
  // weight & maturity

  DATA_INTEGER(lw_int);
  DATA_INTEGER(lw_slope);
  DATA_INTEGER(mat_b0);
  DATA_INTEGER(mat_b1);
  
  // natural mortality prior (assume normal prior with mu and sigma)
  DATA_SCALAR(mu_M)           
  DATA_SCALAR(sigma_M)        

  // parameters ----
  
  //von b
  PARAMETER(Linf);
  PARAMETER(kappa);
  PARAMETER(t0)
  PARAMETER(log_von_sigma);
  
  PARAMETER(Fcur);
  PARAMETER(M);
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER(log_sel_sigma);
  
  // procedures (transform parameters) -----
  Type von_sigma = exp(log_von_sigma);
  Type sel_sigma = exp(log_sel_sigma);
  
  // age-based natural mortality M
  vector<Type> L(n);
  vector<Type> log_L(n);
  //vector<Type> M(n);
  //Type log_Linf = log(Linf);
  //Type log_kappa = log(kappa);

  // vectors for predictions
  // von b
  vector<Type> yfit(von_n); 
  
  // selectivity
  vector<Type> Na(n);
  vector<Type> UnNa(n); //unfished abundance
  vector<Type> Sa(n);
  vector<Type> Ca(n);
  vector<Type> fit_prop(n);
  
  // fishing mortality at age
  vector<Type> Fa(n);
  vector<Type> expN(n);
  vector<Type> totBio(n);
  vector<Type> expBio(n);

  // weight 
  vector<Type> Wa(n);
  vector<Type> Mat(n);
  vector<Type> CWa(n);
  
  vector<Type> unfished(n);
  vector<Type> fished(n); 
  
  // negative log likelihood - set to 0 ----
  // von b
  Type von_nll = 0.0;
  
  // selectivity
  Type nll = 0.0;
  Type nll1 = 0.0;
 
  // total nll
  Type tot_nll = 0.0;
  
  // von b -------------------------------------------------------------
  for(int i=0; i<von_n; i++){
    yfit(i) = Linf * (Type(1.0) - exp(-kappa * (age(i) - t0)));
  }
  
  von_nll = -sum(dnorm(length, yfit, von_sigma, true));
  
  // M --------------------------------------------------------------------

  for(int i=0; i<n; i++){
    L(i) = Linf * (Type(1.0) - exp(-kappa * (X(i) - t0)));
  }

  for(int i=0; i<n; i++){
  log_L(i) = log(L(i));
	}
  
// for(int i=0; i<n; i++){
//  M(i) = exp(Type(0.55) - Type(1.61) * log_L(i) + Type(1.44) * log_Linf + log_kappa);
//  }
  
  // weight at age
  for(int i=0; i<n; i++){
  	Wa(i) = exp(lw_int + lw_slope * log_L(i));
  }

  
  //selectivity ----------------------------------------------------------------
  
  for(int i=0; i<n; i++){
    Sa(i) = Type(1.0) / (Type(1.0) + exp(alpha + beta * X(i) ));
  } 
  
  nll = -sum(dnorm(prop, Sa, sel_sigma, true));
  
//  for(int i=0; i<n; i++){
//  	Fa(i) = Sa(i) * Fcur;
//  }
 
  for(int j=0; j<n; j++){
    
    if(j == 0) Na(j) = 1000;
    if(j > 0) 
      Na(j) = Na(j-1) * exp(-(M + Sa(j-1) * Fcur));
  }

  for(int j=0; j<n; j++){
    
    if(j == 0) UnNa(j) = 1000;
    if(j > 0) 
      UnNa(j) = UnNa(j-1) * exp(-M);
  }

  // maturity at age -----------------------------------------------------------

  for(int i=0; i<n; i++){
  	Mat(i) = exp(mat_b0 + mat_b1 * X(i)) / (Type(1.0) + exp(mat_b0 + mat_b1 * X(i)));
  }

  // catch at age --------------------------------------------------------------

  for(int i=0; i<n; i++){
    Ca(i) = Na(i) * (1 - exp(-(M + Sa(i) * Fcur))) * Fcur * Sa(i) / (M + Sa(i) * Fcur);
   }
  
  Type maxCa = max(Ca);
  
  for(int i = 0; i<n; i++){
    fit_prop(i) = Ca(i) / maxCa; 
  }
  
  nll1 = -sum(dnorm(prop, fit_prop, sel_sigma, true));

  // Prior for natural mortality
  Type prior_M = 0;
  prior_M += Type(0.5) * pow(M - mu_M, 2) / pow(sigma_M, 2);
  
  // biomass estimates --------------------
	for(int i = 0; i<n; i++){
		CWa(i) = Ca(i) * Wa(i);
		unfished(i) = UnNa(i) * Wa(i) * Mat(i);
		fished(i) = Na(i) * Wa(i) * Mat(i);
		Fa(i) = Sa(i) * Fcur;
		expN(i) = Na(i) * Sa(i);
		totBio(i) = Na(i) * Wa(i);
		expBio(i) = expN(i) * Wa(i);
	}
  
  Type exploitN = sum(Ca) / sum(expN);
	Type exploitB = sum(CWa) / sum(expBio);
  
  tot_nll += nll;
  tot_nll += nll1;
  tot_nll += von_nll;
  tot_nll += prior_M;
  
  // reports ---- --------------------------------------------------------------------
  
  //von b
  
  REPORT(Linf);
  REPORT(kappa);
  REPORT(t0);
  REPORT(yfit);
  REPORT(von_sigma);
  REPORT(von_nll);
  
  REPORT(M);
  //REPORT(Fa);
  REPORT(Mat);
  
  REPORT(alpha);
  REPORT(beta);
  REPORT(sel_sigma);
  REPORT(nll);
  REPORT(nll1);
  REPORT(fit_prop);
  REPORT(Sa);
  REPORT(Na);
  REPORT(UnNa);
  REPORT(Ca);
  REPORT(Fcur);
  REPORT(M);
  REPORT(CWa);
  REPORT(unfished);
  REPORT(fished);
  REPORT(exploitN);
  REPORT(exploitB)
  
  return tot_nll;
}


