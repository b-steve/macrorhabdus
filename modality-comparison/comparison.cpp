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
  DATA_INTEGER(n_chickens);
  DATA_VECTOR(times);
  DATA_VECTOR(all_times);
  DATA_IVECTOR(time_indices);
  DATA_VECTOR(end_treatment);
  DATA_INTEGER(cov_function);
  DATA_INTEGER(time_relationship);
  DATA_IVECTOR(n_measurements);
  DATA_MATRIX(dwm);
  DATA_MATRIX(ef);
  DATA_MATRIX(egs);
  DATA_MATRIX(emb);
  DATA_MATRIX(gs);
  DATA_IVECTOR(infected_chickens);
  DATA_VECTOR(test_cases);
  int n_test_cases = test_cases.size();
  PARAMETER(beta_0_base);
  PARAMETER_VECTOR(beta_0_diff);
  PARAMETER(beta_1_base);
  PARAMETER_VECTOR(beta_1_diff);
  PARAMETER(alpha_1);
  PARAMETER(alpha_2);
  PARAMETER(log_theta_base);
  PARAMETER_VECTOR(log_theta_diff);
  PARAMETER(log_sigma_s);
  PARAMETER(log_rho_s);
  PARAMETER(log_sigma_u1);
  PARAMETER(log_sigma_u2);
  // Creating some parameter vectors.
  vector<Type> betas_0(5);
  vector<Type> betas_1(5);
  vector<Type> thetas(5);
  betas_0(0) = beta_0_base;
  betas_1(0) = beta_1_base;
  thetas(0) = exp(log_theta_base);
  for (int i = 1; i < 5; i++){
    betas_0(i) = beta_0_base + beta_0_diff(i - 1);
    betas_1(i) = beta_1_base + beta_1_diff(i - 1);
    thetas(i) = exp(log_theta_base + log_theta_diff(i - 1));
  }
  // Constructing some comparisons between treatments.
  matrix<Type> beta_0_comparison(5, 5);
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 5; j++){
      beta_0_comparison(i, j) = betas_0(i) - betas_0(j);
    }
  }
  // Back-transforming some other parameters.
  Type sigma_s = exp(log_sigma_s);
  Type rho_s = exp(log_rho_s);
  Type sigma_u1 = exp(log_sigma_u1);
  Type sigma_u2 = exp(log_sigma_u2);
  // Random effects.
  PARAMETER_MATRIX(s);
  PARAMETER_VECTOR(u_1);
  PARAMETER_VECTOR(u_2);
  // Total number of times.
  int n_all_times = all_times.size();
  // A scalar for end of treatment.
  Type use_end_treatment;
  // Setting up expectation vector of random effects.
  matrix<Type> mu_s(n_chickens, n_all_times);
  for (int i = 0; i < n_chickens; i++){
    for (int j = 0; j < n_all_times; j++){
    if (time_relationship == 1){
      // Linear.
      mu_s(i, j) = alpha_1*all_times(j);
    } else if (time_relationship == 2){
      // Quadratic.
      // This doesn't account for different treatment lengths, so
      // throws an error for now.
      mu_s(i, j) = alpha_1*all_times(j) + alpha_2*pow(all_times(j), 2);
      exit(2222);
    } else if (time_relationship == 3){
      // Piecewise exponential.
      use_end_treatment = end_treatment(1 - infected_chickens(i));
	if (all_times(j) > use_end_treatment){
	  mu_s(i, j) = (alpha_1 + u_1(i))*use_end_treatment + (alpha_2 + u_2(i))*(all_times(j) - use_end_treatment);
	} else {
	  mu_s(i, j) = (alpha_1 + u_1(i))*all_times(j);
	}
      }
    }
  }
  // Initialising matrix of conditional expectations at the final treatment.
  matrix<Type> lambda_y_mat(n_chickens, 5);
  // Initialising matrix of probabilities of zero detected organisms at final treatment.
  Type p_zero;
  matrix<Type> logit_p_zero(n_test_cases, 5);
  // Initialising conditional expectation of MO count.
  Type lambda_y;
  // Initialising count for Type-casting.
  Type y;
  // Initialising variance-covariance matrix for latent variables.
  matrix<Type> sigma_s_mat(n_all_times, n_all_times);
  // Declaring log of the joint density.
  Type ljd = 0;
  // Looping over individuals.
  for (int i = 0; i < n_chickens; i++){
    // Looping over time points.
    for (int j = 0; j < n_measurements(i); j++){
      // Looping over methods.
      for (int k = 0; k < 5; k++){
  	if (k == 0){
  	  y = dwm(i, j);
  	} else if (k == 1){
  	  y = ef(i, j);
  	} else if (k == 2){
  	  y = egs(i, j);
  	} else if (k == 3){
  	  y = emb(i, j);
  	} else if (k == 4){
  	  y = gs(i, j);
  	}
  	// Calculating conditional expectation.
  	lambda_y = exp(betas_0(k) + betas_1(k)*s(i, time_indices(j)));
  	// Contribution to the log of the joint density.
  	ljd += nbinom_pmf(y, lambda_y, thetas(k), true);
  	// Calculating probabilities of zeros for various test-cases.
	if (j > 4 && i == 1){
	  for (int z = 0; z < n_test_cases; z++){
	    Type dummy_zero = 0;
	    p_zero = nbinom_pmf(dummy_zero, exp(betas_0(k) + betas_1(k)*test_cases(z)), thetas(k), false);
	    logit_p_zero(z, k) = log(p_zero/(1 - p_zero));
	  }
	}
      }
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
    vector<Type> s_i(n_all_times);
    // Getting detrended latent variables for ith animal.
    for (int j = 0; j < n_all_times; j++){
      s_i(j) = s.row(i)(j) - mu_s(i, j);
    }
    ljd -= MVNORM(sigma_s_mat)(s_i);
    // Contribution from latent slope variables.
    Type dummy_zero = 0;
    ljd += dnorm(u_1(i), dummy_zero, sigma_u1, 1);
    ljd += dnorm(u_2(i), dummy_zero, sigma_u2, 1);
  }
  // Comparing p_zero within rows.
  vector<Type> p_zero_comparison(25*n_test_cases);
  int p_zero_index = 0;
  for (int i = 0; i < n_test_cases; i++){
    for (int j = 0; j < 5; j++){
      for (int k = 0; k < 5; k++){
	p_zero_comparison(p_zero_index) = logit_p_zero(i, j) - logit_p_zero(i, k);
	p_zero_index++;
      }
    }
  }
  // ADREPORTS.
  ADREPORT(thetas);
  ADREPORT(sigma_s);
  ADREPORT(rho_s);
  ADREPORT(mu_s);
  ADREPORT(logit_p_zero);
  ADREPORT(beta_0_comparison);
  ADREPORT(p_zero_comparison);
  std::cout << "LJD: " << ljd << std::endl;
  ljd *= -1;
  return ljd;
}
