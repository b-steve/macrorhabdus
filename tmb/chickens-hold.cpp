// Poisson generalised linear mixed-effects model for chicken data.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Reading in data.
  DATA_INTEGER(n);
  DATA_MATRIX(mo_count);
  DATA_VECTOR(weight);
  DATA_IVECTOR(group);
  DATA_VECTOR(weight_points);
  // Declaring parameters.
  PARAMETER(beta_0);
  PARAMETER(beta_weight);
  PARAMETER(beta_weight_2);
  PARAMETER(beta_g);
  PARAMETER(log_sigma_u);
  PARAMETER(log_sigma_v);
  Type sigma_u = exp(log_sigma_u);
  Type sigma_v = exp(log_sigma_v);
  ADREPORT(sigma_u);
  ADREPORT(sigma_v);
  // Declaring random effects.
  PARAMETER_VECTOR(u);
  PARAMETER_MATRIX(v);
  // Sorting out coefficients for groups.
  vector<Type> beta_group(3);
  beta_group(0) = 0;
  beta_group(1) = beta_g;
  beta_group(2) = beta_g;
  // Declaring log-likelihood.
  Type ll = 0;
  // Declaring conditional expectation.
  Type lambda = 0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < 10; j++){
      // Calculating conditional expectation.
      lambda = exp(beta_0 + beta_weight*weight(i) + beta_weight_2*pow(weight(i), 2) + beta_group(group(i)) + u(i) + v(i,j));
      // Contribution from MO count.
      ll += dpois(mo_count(i,j), lambda, true);
      // Contribution from slide latent variable.
      ll += dnorm(v(i, j), Type(0), sigma_v, true);
    }
    // Conbribution from chicken latent variable.
    ll += dnorm(u(i), Type(0), sigma_u, true);
  }
  // Getting some fitted lines.
  vector<Type> ntc_preds(1000);
  vector<Type> d25_preds(1000);
  vector<Type> d100_preds(1000);
  for (int i = 0; i < 1000; i++){
    ntc_preds(i) = beta_0 + beta_weight*weight_points(i) + beta_weight_2*pow(weight_points(i), 2) + beta_group(0);
    d25_preds(i) = beta_0 + beta_weight*weight_points(i) + beta_weight_2*pow(weight_points(i), 2) + beta_group(1);
    d100_preds(i) = beta_0 + beta_weight*weight_points(i) + beta_weight_2*pow(weight_points(i), 2) + beta_group(2);
  }
  ADREPORT(ntc_preds);
  ADREPORT(d25_preds);
  ADREPORT(d100_preds);
  return -ll;
}
