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
  DATA_VECTOR(weight_end);
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
  PARAMETER(beta_wc);
  PARAMETER(alpha_1);
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
  // Initialising overall treatment effects.
  vector<Type> t_effect(n_birds);
  // Initialising log_mu for the typical bird.
  vector<Type> log_mu_typical(n_all_times);
  // Initialising treatment effect for the typical bird.
  Type t_effect_typical;
  for (int i = 0; i < n_birds; i++){
    for (int j = 0; j < n_all_times; j++){
      if (end_treatment(i) == 10){
	final_beta = beta_2;
      } else {
	final_beta = beta_3;
      }
      if (all_times(j) < end_treatment(i)){
	log_mu(i, j) = beta_0 - beta_1*all_times(j) + beta_w*weight_start(i) + beta_wt*all_times(j)*weight_start(i) + s(i, j);
	// Interaction effect of intercept.
	log_mu(i, j) += alpha_1*all_times(j)*log_mu(i, 0);
      } else {
	log_mu(i, j) = beta_0 - beta_1*end_treatment(i) + beta_w*weight_start(i) + beta_wt*end_treatment(i)*weight_start(i) + final_beta*(all_times(j) - end_treatment(i)) + s(i, j);
	// Interaction effect of intercept.
	log_mu(i, j) += alpha_1*end_treatment(i)*log_mu(i, 0);
      }
      if (all_times(j) < 10){
	if (i == 0){
	  log_mu_typical(j) = beta_0 - beta_1*all_times(j);
	  log_mu_typical(j) += alpha_1*all_times(j)*log_mu_typical(0);
	}
      } else {
	if (i == 0){
	  log_mu_typical(j) = beta_0 - beta_1*end_treatment(i) + beta_2*(all_times(j) - end_treatment(i));
	  log_mu_typical(j) += alpha_1*end_treatment(i)*log_mu_typical(0);
	}
      }
    }
    // Saving overall treatment effects.
    t_effect(i) = -beta_1 + beta_wt*weight_start(i) + alpha_1*log_mu(i, 0);
    // Adding in effect of weight change on final expectation.
    log_mu(i, n_all_times - 1) += beta_wc*(weight_end(i) - weight_start(i));
  }
  // Typical treatment effect.
  t_effect_typical = -beta_1 + alpha_1*log_mu_typical(0);
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
    //std::cout << sigma_s_mat << std::endl; exit(2123);
    ljd -= MVNORM(sigma_s_mat)(s.row(i));
  }
  // ADREPORTS.
  ADREPORT(theta);
  ADREPORT(rho_s);
  ADREPORT(sigma_s);
  ADREPORT(log_mu);
  ADREPORT(log_mu_typical);
  ADREPORT(t_effect);
  ADREPORT(t_effect_typical);
  std::cout << "LJD: " << ljd << std::endl;
  ljd *= -1;
  return ljd;
}
