// generated with brms 2.15.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] se;  // known sampling error
  int<lower=1> K_rtref;  // number of population-level effects
  matrix[N, K_rtref] X_rtref;  // population-level design matrix
  int<lower=1> K_e;  // number of population-level effects
  matrix[N, K_e] X_e;  // population-level design matrix
  int<lower=1> K_eh;  // number of population-level effects
  matrix[N, K_eh] X_eh;  // population-level design matrix
  int<lower=1> K_th;  // number of population-level effects
  matrix[N, K_th] X_th;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
}
parameters {
  vector[K_rtref] b_rtref;  // population-level effects
  vector[K_e] b_e;  // population-level effects
  vector[K_eh] b_eh;  // population-level effects
  vector[K_th] b_th;  // population-level effects
  real<lower=0> sigma;  // residual SD
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_rtref = X_rtref * b_rtref;
    // initialize linear predictor term
    vector[N] nlp_e = X_e * b_e;
    // initialize linear predictor term
    vector[N] nlp_eh = X_eh * b_eh;
    // initialize linear predictor term
    vector[N] nlp_th = X_th * b_th;
    // initialize non-linear predictor term
    vector[N] mu;
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (exp(nlp_rtref[n]) * exp(exp(nlp_e[n]) / C_1[n] * (1 / C_2[n] - 1 / (C_3[n] + 273.15)))) * (1 / (1 + exp(exp(nlp_eh[n]) / C_1[n] * (1 / (nlp_th[n] + 273.15) - 1 / (C_3[n] + 273.15)))));
    }
    target += normal_lpdf(Y | mu, sqrt(square(sigma) + se2));
  }
  // priors including constants
  target += normal_lpdf(b_rtref | log(0.2), 0.3);
  target += normal_lpdf(b_e | log(0.5), 0.3);
  target += normal_lpdf(b_eh | log(3), 0.4);
  target += normal_lpdf(b_th | 35, 3);
  target += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
generated quantities {
}