// src/core.cpp
// High-performance simulation step for SIMEX and MC-SIMEX
//
// Contains:
//   - IRLS GLM solver (Gaussian/Poisson/Binomial)
//   - Matrix exponentiation via eigendecomposition (for MC-SIMEX)
//   - Multinomial resampling (for MC-SIMEX)
//   - simex_sim_cpp: SIMEX simulation loop for continuous measurement error
//   - mcsimex_sim_cpp: MC-SIMEX simulation loop for misclassification

#include <RcppEigen.h>
#include <random>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

// =========================================================================
// IRLS GLM solver
// =========================================================================
// dist_code: 1 = Gaussian (identity), 2 = Poisson (log), 3 = Binomial (logit)

// Inverse link functions
inline double mu_gaussian(double eta) { return eta; }
inline double mu_poisson(double eta)  { return std::exp(eta); }
inline double mu_binomial(double eta) { return 1.0 / (1.0 + std::exp(-eta)); }

// d mu / d eta
inline double mu_dot_gaussian(double)    { return 1.0; }
inline double mu_dot_poisson(double eta) { return std::exp(eta); }
inline double mu_dot_binomial(double eta) {
  double p = 1.0 / (1.0 + std::exp(-eta));
  return p * (1.0 - p);
}

// Variance function V(mu)
inline double var_gaussian(double)   { return 1.0; }
inline double var_poisson(double mu) { return mu; }
inline double var_binomial(double mu) { return mu * (1.0 - mu); }

VectorXd irls_fit(const VectorXd& y, const MatrixXd& X,
                  int dist_code, const VectorXd& wt,
                  int max_iter = 25, double tol = 1e-8) {
  int n = X.rows();
  int p = X.cols();

  double (*mu_fn)(double);
  double (*mu_dot_fn)(double);
  double (*var_fn)(double);

  switch (dist_code) {
    case 1: mu_fn = mu_gaussian; mu_dot_fn = mu_dot_gaussian;
            var_fn = var_gaussian; break;
    case 2: mu_fn = mu_poisson;  mu_dot_fn = mu_dot_poisson;
            var_fn = var_poisson;  break;
    case 3: mu_fn = mu_binomial; mu_dot_fn = mu_dot_binomial;
            var_fn = var_binomial; break;
    default: Rcpp::stop("Unknown dist_code: %d", dist_code);
  }

  VectorXd beta = VectorXd::Zero(p);

  // Warm start: initialize intercept
  if (dist_code == 2) {
    double ymean = 0, wtsum = 0;
    for (int i = 0; i < n; i++) { ymean += wt(i) * y(i); wtsum += wt(i); }
    ymean /= wtsum;
    if (ymean > 0) beta(0) = std::log(ymean);
  }
  if (dist_code == 3) {
    double ymean = 0, wtsum = 0;
    for (int i = 0; i < n; i++) { ymean += wt(i) * y(i); wtsum += wt(i); }
    ymean /= wtsum;
    ymean = std::max(0.01, std::min(0.99, ymean));
    beta(0) = std::log(ymean / (1.0 - ymean));
  }

  VectorXd eta(n), w_irls(n), z_irls(n);

  for (int iter = 0; iter < max_iter; iter++) {
    eta = X * beta;
    for (int i = 0; i < n; i++) {
      double mu_i  = mu_fn(eta(i));
      double md    = mu_dot_fn(eta(i));
      double v     = var_fn(mu_i);
      w_irls(i)    = wt(i) * md * md / std::max(v, 1e-10);
      z_irls(i)    = eta(i) + (y(i) - mu_i) / std::max(std::abs(md), 1e-10);
    }

    MatrixXd Xw = X.array().colwise() * w_irls.array();
    MatrixXd XtWX(p, p);
    VectorXd XtWz(p);
    XtWX.noalias() = Xw.transpose() * X;
    XtWz.noalias() = Xw.transpose() * z_irls;

    VectorXd beta_new = XtWX.ldlt().solve(XtWz);

    double change = (beta_new - beta).cwiseAbs().maxCoeff();
    beta = beta_new;
    if (change < tol) break;
  }

  return beta;
}


// =========================================================================
// Matrix exponentiation: Pi^lambda via eigendecomposition
// =========================================================================
// Pi = V * diag(d) * V^{-1}  =>  Pi^lambda = V * diag(d^lambda) * V^{-1}

MatrixXd mat_power(const MatrixXd& Pi, double lambda) {
  int K = Pi.rows();
  if (lambda == 0.0) return MatrixXd::Identity(K, K);
  if (lambda == 1.0) return Pi;

  Eigen::EigenSolver<MatrixXd> es(Pi);
  MatrixXd V  = es.eigenvectors().real();
  VectorXd d  = es.eigenvalues().real();

  VectorXd d_pow(K);
  for (int k = 0; k < K; k++) {
    d_pow(k) = std::pow(std::abs(d(k)), lambda);
    if (d(k) < 0) d_pow(k) = -d_pow(k);
  }

  return V * d_pow.asDiagonal() * V.inverse();
}


// =========================================================================
// Multinomial resampling: sample z from column k of Pi^lambda
// =========================================================================

VectorXi resample_z(const VectorXi& z_hat, const MatrixXd& Pi_lam,
                    std::mt19937& rng) {
  int n = z_hat.size();
  int K = Pi_lam.rows();
  VectorXi z_new(n);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for (int i = 0; i < n; i++) {
    int k = z_hat(i);
    double u = unif(rng);
    double cumprob = 0.0;
    int cat = K - 1;
    for (int j = 0; j < K; j++) {
      cumprob += Pi_lam(j, k);
      if (u <= cumprob) { cat = j; break; }
    }
    z_new(i) = cat;
  }
  return z_new;
}


// =========================================================================
// Build design matrix [dummies(z, K-1), x] for MC-SIMEX
// =========================================================================

MatrixXd build_design(const VectorXi& z, const MatrixXd& x, int K) {
  int n = z.size();
  int s = K - 1;
  int r = x.cols();
  int p = s + r;
  MatrixXd xi(n, p);

  for (int k = 0; k < s; k++) {
    for (int i = 0; i < n; i++) {
      xi(i, k) = (z(i) == k + 1) ? 1.0 : 0.0;
    }
  }
  xi.rightCols(r) = x;
  return xi;
}


// =========================================================================
// MC-SIMEX simulation step (exported to R)
// =========================================================================
// Returns matrix of dimension (B * n_lambda) x p

// [[Rcpp::export]]
Rcpp::NumericMatrix mcsimex_sim_cpp(
    Rcpp::NumericVector y_r,
    Rcpp::IntegerVector z_hat_r,
    Rcpp::NumericMatrix x_r,
    Rcpp::NumericMatrix Pi_r,
    int K,
    int dist_code,
    Rcpp::NumericVector lambda_r,
    int B,
    Rcpp::NumericVector wt_r,
    unsigned int seed) {

  int n = y_r.size();
  int r = x_r.ncol();
  int s = K - 1;
  int p = s + r;
  int n_lambda = lambda_r.size();

  Eigen::Map<VectorXd> y(y_r.begin(), n);
  Eigen::Map<VectorXi> z_hat(z_hat_r.begin(), n);
  Eigen::Map<MatrixXd> x(x_r.begin(), n, r);
  Eigen::Map<MatrixXd> Pi(Pi_r.begin(), K, K);
  Eigen::Map<VectorXd> wt(wt_r.begin(), n);

  // Pre-compute Pi^lambda for each lambda
  std::vector<MatrixXd> Pi_powers(n_lambda);
  for (int l = 0; l < n_lambda; l++) {
    Pi_powers[l] = mat_power(Pi, lambda_r[l]);
  }

  Rcpp::NumericMatrix result(n_lambda * B, p);
  std::mt19937 rng(seed);

  for (int l = 0; l < n_lambda; l++) {
    const MatrixXd& Pi_lam = Pi_powers[l];

    for (int b = 0; b < B; b++) {
      VectorXi z_sim = resample_z(z_hat, Pi_lam, rng);
      MatrixXd xi = build_design(z_sim, x, K);
      VectorXd beta = irls_fit(y, xi, dist_code, wt);

      int row = l * B + b;
      for (int j = 0; j < p; j++) {
        result(row, j) = beta(j);
      }
    }

    Rcpp::checkUserInterrupt();
  }

  return result;
}


// =========================================================================
// Multi-variable MC-SIMEX simulation step (exported to R)
// =========================================================================
// Handles multiple misclassified variables simultaneously.
// z_list:  list of integer vectors (0-based factor codes), one per mc variable
// Pi_list: list of K_j x K_j misclassification matrices
// K_vec:   integer vector of number of levels per mc variable
// x_r:     other covariates matrix (n x r), including intercept
//
// Design matrix: [dummies_z1(K1-1) | dummies_z2(K2-1) | ... | x_other]
// Returns matrix of dimension (B * n_lambda) x p

// [[Rcpp::export]]
Rcpp::NumericMatrix mcsimex_multi_sim_cpp(
    Rcpp::NumericVector y_r,
    Rcpp::List z_list,
    Rcpp::List Pi_list,
    Rcpp::IntegerVector K_vec,
    Rcpp::NumericMatrix x_r,
    int dist_code,
    Rcpp::NumericVector lambda_r,
    int B,
    Rcpp::NumericVector wt_r,
    unsigned int seed) {

  int n = y_r.size();
  int n_mc = z_list.size();
  int r = x_r.ncol();
  int n_lambda = lambda_r.size();

  // Compute total dummy columns: sum(K_j - 1)
  int total_dummies = 0;
  for (int j = 0; j < n_mc; j++) {
    total_dummies += K_vec[j] - 1;
  }
  int p = total_dummies + r;

  Eigen::Map<VectorXd> y(y_r.begin(), n);
  Eigen::Map<MatrixXd> x(x_r.begin(), n, r);
  Eigen::Map<VectorXd> wt(wt_r.begin(), n);

  // Map z vectors and Pi matrices
  std::vector<Eigen::Map<VectorXi>> z_hats;
  std::vector<Eigen::Map<MatrixXd>> Pis;
  z_hats.reserve(n_mc);
  Pis.reserve(n_mc);

  for (int j = 0; j < n_mc; j++) {
    Rcpp::IntegerVector zj = z_list[j];
    Rcpp::NumericMatrix Pj = Pi_list[j];
    z_hats.push_back(Eigen::Map<VectorXi>(zj.begin(), n));
    Pis.push_back(Eigen::Map<MatrixXd>(Pj.begin(), K_vec[j], K_vec[j]));
  }

  // Pre-compute Pi^lambda for each (variable, lambda)
  // Pi_powers[j][l] = Pi_j^{lambda_l}
  std::vector<std::vector<MatrixXd>> Pi_powers(n_mc);
  for (int j = 0; j < n_mc; j++) {
    Pi_powers[j].resize(n_lambda);
    for (int l = 0; l < n_lambda; l++) {
      Pi_powers[j][l] = mat_power(MatrixXd(Pis[j]), lambda_r[l]);
    }
  }

  Rcpp::NumericMatrix result(n_lambda * B, p);
  std::mt19937 rng(seed);

  for (int l = 0; l < n_lambda; l++) {
    for (int b = 0; b < B; b++) {
      // Build design matrix: resample each mc variable, create dummies
      MatrixXd xi(n, p);

      int col_offset = 0;
      for (int j = 0; j < n_mc; j++) {
        int Kj = K_vec[j];
        int sj = Kj - 1;
        VectorXi z_sim = resample_z(z_hats[j], Pi_powers[j][l], rng);

        // Write dummies for variable j
        for (int k = 0; k < sj; k++) {
          for (int i = 0; i < n; i++) {
            xi(i, col_offset + k) = (z_sim(i) == k + 1) ? 1.0 : 0.0;
          }
        }
        col_offset += sj;
      }

      // Copy other covariates
      xi.rightCols(r) = x;

      VectorXd beta = irls_fit(y, xi, dist_code, wt);

      int row = l * B + b;
      for (int jj = 0; jj < p; jj++) {
        result(row, jj) = beta(jj);
      }
    }

    Rcpp::checkUserInterrupt();
  }

  return result;
}


// =========================================================================
// SIMEX simulation step for continuous measurement error (exported to R)
// =========================================================================
// X_r:          full model matrix (n x p) from the naive model
// simex_cols_r: 0-based column indices in X corresponding to SIMEXvariables
// me_r:         measurement error SDs (n x n_simex_vars)
//
// For each (lambda, b):
//   X_sim = X + sqrt(lambda) * N(0,1) * me  (only in simex columns)
//   fit IRLS on (y, X_sim)

// Compute GLM model-based vcov = (X'WX)^{-1} * phi at given beta
// Uses the same link/variance functions as irls_fit
MatrixXd glm_vcov(const VectorXd& beta, const MatrixXd& X,
                  const VectorXd& y, int dist_code,
                  const VectorXd& wt) {
  int n = X.rows();
  int p = X.cols();

  double (*mu_fn)(double);
  double (*mu_dot_fn)(double);
  double (*var_fn)(double);

  switch (dist_code) {
    case 1: mu_fn = mu_gaussian; mu_dot_fn = mu_dot_gaussian;
            var_fn = var_gaussian; break;
    case 2: mu_fn = mu_poisson;  mu_dot_fn = mu_dot_poisson;
            var_fn = var_poisson;  break;
    case 3: mu_fn = mu_binomial; mu_dot_fn = mu_dot_binomial;
            var_fn = var_binomial; break;
    default: mu_fn = mu_gaussian; mu_dot_fn = mu_dot_gaussian;
             var_fn = var_gaussian; break;
  }

  VectorXd eta = X * beta;
  VectorXd w_irls(n);
  double wt_sum = 0;
  double pearson_sum = 0;

  for (int i = 0; i < n; i++) {
    double mu_i = mu_fn(eta(i));
    double md = mu_dot_fn(eta(i));
    double v = var_fn(mu_i);
    w_irls(i) = wt(i) * md * md / std::max(v, 1e-10);
    wt_sum += wt(i);

    double resid = (y(i) - mu_i) / std::max(std::sqrt(v), 1e-10);
    pearson_sum += wt(i) * resid * resid;
  }

  MatrixXd XtWX = (X.array().colwise() * w_irls.array()).matrix().transpose() * X;
  MatrixXd XtWX_inv = XtWX.ldlt().solve(MatrixXd::Identity(p, p));

  // Dispersion: 1 for Poisson/Binomial, Pearson estimate for Gaussian
  double phi = (dist_code == 1) ? pearson_sum / (wt_sum - p) : 1.0;

  return XtWX_inv * phi;
}


// [[Rcpp::export]]
Rcpp::List simex_sim_cpp(
    Rcpp::NumericVector y_r,
    Rcpp::NumericMatrix X_r,
    Rcpp::IntegerVector simex_cols_r,
    Rcpp::NumericMatrix me_r,
    int dist_code,
    Rcpp::NumericVector lambda_r,
    int B,
    Rcpp::NumericVector wt_r,
    unsigned int seed) {

  int n = y_r.size();
  int p = X_r.ncol();
  int n_lambda = lambda_r.size();
  int n_simex = simex_cols_r.size();

  Eigen::Map<VectorXd> y(y_r.begin(), n);
  Eigen::Map<MatrixXd> X(X_r.begin(), n, p);
  Eigen::Map<VectorXd> wt(wt_r.begin(), n);

  // Map measurement error matrix
  Eigen::Map<MatrixXd> me(me_r.begin(), n, n_simex);

  // Copy simex column indices
  std::vector<int> scols(n_simex);
  for (int j = 0; j < n_simex; j++) {
    scols[j] = simex_cols_r[j];
  }

  Rcpp::NumericMatrix result(n_lambda * B, p);
  // Accumulated model vcov per lambda: stored as n_lambda rows x p*p columns
  Rcpp::NumericMatrix vcov_acc(n_lambda, p * p);

  std::mt19937 rng(seed);
  std::normal_distribution<double> rnorm(0.0, 1.0);

  for (int l = 0; l < n_lambda; l++) {
    double sqrt_lam = std::sqrt(lambda_r[l]);
    MatrixXd vcov_sum = MatrixXd::Zero(p, p);

    for (int b = 0; b < B; b++) {
      // Copy model matrix and add noise to SIMEX columns
      MatrixXd X_sim = X;

      for (int j = 0; j < n_simex; j++) {
        int col = scols[j];
        for (int i = 0; i < n; i++) {
          X_sim(i, col) += sqrt_lam * rnorm(rng) * me(i, j);
        }
      }

      VectorXd beta = irls_fit(y, X_sim, dist_code, wt);

      int row = l * B + b;
      for (int jj = 0; jj < p; jj++) {
        result(row, jj) = beta(jj);
      }

      // Accumulate model vcov from this replicate
      vcov_sum += glm_vcov(beta, X_sim, y, dist_code, wt);
    }

    // Store average vcov for this lambda
    MatrixXd vcov_avg = vcov_sum / B;
    for (int j = 0; j < p * p; j++) {
      vcov_acc(l, j) = vcov_avg(j % p, j / p);
    }

    Rcpp::checkUserInterrupt();
  }

  return Rcpp::List::create(
    Rcpp::Named("theta") = result,
    Rcpp::Named("vcov_model") = vcov_acc
  );
}
