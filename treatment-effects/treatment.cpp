#include <TMB.hpp>
using namespace density;

// Negative binomal PMF.
template<class Type>
Type nbinom_pmf (const Type &x, const Type &mu, const Type &theta, const int &give_log){
  Type out;
  out = lgamma(x + theta) - lgamma(theta) - lgamma(x + 1) + theta*(log(theta) - log(theta + mu)) + x*log(1 - theta/(theta + mu));
  if (!give_log){
    out = exp(out);
  }
  return out;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Reading in data.
  DATA_INTEGER(n_birds);
  DATA_VECTOR(times);
  DATA_VECTOR(all_times);
  DATA_IVECTOR(time_indices);
  DATA_VECTOR(end_treatment);
  DATA_MATRIX(mo_counts);
  DATA_IVECTOR(n_measurements);
  DATA_VECTOR(weight_start);
  DATA_INTEGER(cov_function);
  // Total number of times.
  int n_all_times = all_times.size();
  // Parameters.
  PARAMETER(beta_0);
  PARAMETER(beta_1);
  PARAMETER(beta_2);
  PARAMETER(beta_3);
  PARAMETER(beta_w);
  PARAMETER(beta_wt);
  PARAMETER(log_theta);
  PARAMETER(log_sigma_s);
  PARAMETER(log_rho_s);
  // Random effects.
  PARAMETER_MATRIX(s);
  //PARAMETER_VECTOR(u);
  // Back-transforming parameters.
  Type theta = exp(log_theta);
  Type sigma_s = exp(log_sigma_s);
  Type rho_s = exp(log_rho_s);
  // Creating a matrix of expected values.
  matrix<Type> log_mu(n_birds, n_all_times);
  // Initialising final beta.
  Type final_beta;
  for (int i = 0; i < n_birds; i++){
    for (int j = 0; j < n_all_times; j++){
      if (end_treatment(i) == 10){
	final_beta = beta_2;
      } else {
	final_beta = beta_3;
      }
      if (all_times(j) < end_treatment(i)){
	log_mu(i, j) = beta_0 - beta_1*all_times(j) + beta_w*weight_start(i) + beta_wt*all_times(j)*weight_start(i) + s(i, j);
      } else {
	log_mu(i, j) = beta_0 - beta_1*end_treatment(i) + beta_w*weight_start(i) + beta_wt*end_treatment(i)*weight_start(i) + final_beta*(all_times(j) - end_treatment(i)) + s(i, j);
      }
    }
  }
  // Initialising expectation.
  Type mu;
  // Initialising variance-covariance matrix for latent variables.
  matrix<Type> sigma_s_mat(n_all_times, n_all_times);
  // Initialising log of the joint density.
  Type ljd = 0;
  // Looping over individuals.
  for (int i = 0; i < n_birds; i++){
    // Looping over time points.
    for (int j = 0; j < n_measurements(i); j++){
      mu = exp(log_mu(i, time_indices(j)));
      ljd += nbinom_pmf(mo_counts(i, j), mu, theta, true);
    }
    // Contribution from latent variables.
    for (int j = 0; j < n_all_times; j++){
      for (int k = j; k < n_all_times; k++){
	if (cov_function == 1){
	  sigma_s_mat(j, k) = pow(sigma_s, 2)*exp(-(all_times(k) - all_times(j))/rho_s);
	  sigma_s_mat(k, j) = pow(sigma_s, 2)*exp(-(all_times(k) - all_times(j))/rho_s);
	} else if (cov_function == 2){
	  sigma_s_mat(j, k) = pow(sigma_s, 2)*exp(-pow(all_times(k) - all_times(j), 2)/pow(rho_s, 2));
	  sigma_s_mat(k, j) = pow(sigma_s, 2)*exp(-pow(all_times(k) - all_times(j), 2)/pow(rho_s, 2));
	}
      }
    }
    ljd -= MVNORM(sigma_s_mat)(s.row(i));
  }
  // ADREPORTS.
  ADREPORT(theta);
  //ADREPORT(rho_s);
  //ADREPORT(sigma_s);
  ADREPORT(log_mu);
  std::cout << "LJD: " << ljd << std::endl;
  ljd *= -1;
  return ljd;
}
