# Sharpe-Schoolfield 4p fit with bmrs
# IR Moodie

# clear env ---------------------------------------------------------------

rm(list = ls())


# load required packages --------------------------------------------------

library(brms)
library(rstan)
library(cowplot)
library(tidybayes)
library(ggmcmc)
library(tidyverse)


# options -----------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# script specific functions -----------------------------------------------

prior_dist_gen <- function(length, rtref_prior, e_prior, eh_prior, th_prior){
  return(
    tibble(
      b_rtref_Intercept = exp(rnorm(1:length,log(rtref_prior[1]), rtref_prior[2])),
      b_e_Intercept = exp(rnorm(1:length,log(e_prior[1]), e_prior[2])),
      b_eh_Intercept = exp(rnorm(1:length,log(eh_prior[1]), eh_prior[2])),
      b_th_Intercept = rnorm(1:length,th_prior[1],th_prior[2])) %>%
      gather('param', 'estimate', 1:ncol(.)) %>%
      separate(param, c('blah', 'term', 'blah2'), sep = '_') %>%
      select(-starts_with('blah')) %>%
      mutate(type = "prior")
    )
}

get_topt <- function(eh, th, e){
  return((eh * (th+273.15))/(eh + (8.62e-05 * (th+273.15) * log((eh/e) - 1)))-273.15)
}

ss_4_formula <- function (temp, r_tref, e, eh, th, tref) 
{
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- r_tref * exp(e/k * (1/tref - 1/(temp + 
                                                      273.15)))
  inactivation.term <- 1/(1 + exp(eh/k * (1/(th + 273.15) - 
                                            1/(temp + 273.15))))
  return(boltzmann.term * inactivation.term)
}

extract_posteriors <- function(model){
  return(
    model %>%
    spread_draws(`b_.*`, regex = TRUE) %>%
    mutate(b_rtref_Intercept = exp(b_rtref_Intercept),
           b_e_Intercept = exp(b_e_Intercept),
           b_eh_Intercept = exp(b_eh_Intercept),
           b_topt_Intercept = get_topt(eh = b_eh_Intercept, th = b_th_Intercept, e = b_e_Intercept),
           b_rmax_Intercept = ss_4_formula(temp = b_topt_Intercept, r_tref = b_rtref_Intercept, e = b_e_Intercept, eh = b_eh_Intercept, th = b_th_Intercept, tref = 20)) %>%
    gather('param', 'estimate', 4:ncol(.)) %>% 
    separate(param, c('blah', 'term', 'blah2'), sep = '_') %>%
    select(-starts_with('blah')) %>%
    filter(!is.nan(estimate)) %>%
    mutate(replicate = cur_replicate,
           speciesMT = cur_speciesMT,
           selection_temperature = cur_selection_temperature)
  )
}

extract_param_estimates <- function(posteriors){
  return(
      posteriors %>%
      select(-replicate, -speciesMT, -selection_temperature) %>%
      group_by(term) %>%
      mean_qi() %>%
      mutate(replicate = cur_replicate,
             speciesMT = cur_speciesMT,
             selection_temperature = cur_selection_temperature)
  )
}

extract_predictions <- function(model, tmin, tmax, resolution){
  return(
      data.frame(temperature = seq(tmin, tmax, length.out = resolution),
                 tref = 20 + 273.15,
                 k = 8.62e-05,
                 r0_sd = 1) %>%
      add_fitted_draws(model, re_formula = NA) %>%
      data.frame() %>%
      group_by(temperature) %>%
      mean_qi(estimate = .value) %>%
      mutate(replicate = cur_replicate,
             speciesMT = cur_speciesMT,
             selection_temperature = cur_selection_temperature)
  )
}

plot_postprior <- function(posteriors){
  return(
    posteriors %>%
      select(-replicate, -speciesMT, -selection_temperature) %>%
      filter(estimate < 50) %>%
      group_by(term) %>%
      mutate(median = median(estimate),
             type = "post") %>%
      bind_rows(prior_dist) %>%
      ggplot(aes(x = estimate)) +
      facet_wrap(~term, scales = "free") +
      geom_density(aes(fill = type), alpha = 0.3) +
      geom_vline(aes(xintercept = median), linetype = "dashed") +
      scale_fill_manual(values = c("red", "grey")) +
      theme_cowplot() +
      panel_border() +
      theme(legend.position = "none")
  )
}

plot_traceplot <- function(model){
  return(
    model %>%
      ggs() %>%
      ggs_traceplot() +
      theme_cowplot() +
      panel_border()
  )
}

plot_niche <- function(predictions){
  return(
    ggplot(predictions, aes(x = temperature, y = estimate)) +
      geom_line() +
      geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
      geom_point(data = d, aes(x = temperature, y = r0_mean), size  = 3) +
      geom_errorbar(data = d, aes(x = temperature, y = r0_mean, ymin = r0_min, ymax = r0_max)) +
      labs(x = "Temperature",
           y = "Intrinsic rate of increase, r0",
           title = paste(cur_speciesMT, cur_replicate, cur_selection_temperature, sep = " ")) +
      theme_cowplot() +
      background_grid(major = "xy", minor = "x")
  )
}

# load in required datasets -----------------------------------------------

data <- read_csv("exampledata.csv")


# define treatment groups -------------------------------------------------

treatments <-
  data %>%
  select(speciesMT, replicate, selection_temperature) %>%
  distinct()


# build brms model --------------------------------------------------------

ss_4 <- bf(r0_mean|se(r0_sd, sigma = TRUE) ~ (exp(rtref) * exp(exp(e)/k * (1/tref - 1/(temperature + 273.15))))*(1/(1 + exp(exp(eh)/k * (1/(th + 273.15) - 1/(temperature + 273.15))))),
           rtref ~ 1,
           e ~ 1,
           eh ~ 1,
           th ~ 1,
           nl = TRUE)

# start loop for fitting --------------------------------------------------

for (row in 1:nrow(treatments)) {
  
  cur_speciesMT <- as.character((treatments[row, "speciesMT"]))
  cur_replicate <- as.character((treatments[row, "replicate"]))  
  cur_selection_temperature <- as.character((treatments[row, "selection_temperature"]))
  
  d <- 
    data %>%
    filter(speciesMT == cur_speciesMT & replicate == cur_replicate & selection_temperature == cur_selection_temperature)
  
  print(paste("Fitting selection temp", cur_selection_temperature, "replicate", cur_replicate, "of", cur_speciesMT,  sep = " "))
  
  
  # define species specific priors ------------------------------------------
  
  # priors that don't change
  e_prior <- c(0.5, 0.3)
  eh_prior <- c(3, 0.4)
  th_prior <- c(35,3)
  
  # Tet VII
  if (cur_speciesMT == "Tet_VII") {
    
    prior_dist <- prior_dist_gen(3000, rtref_prior = c(0.2, 0.3), e_prior, eh_prior, th_prior)
    
    ss_4_priors <- c(prior(normal(log(0.2), 0.3), nlpar = "rtref"),
                     prior(normal(log(0.5), 0.3), nlpar = "e"),
                     prior(normal(log(3), 0.4), nlpar = "eh"),
                     prior(normal(35,3), nlpar = "th"))
  }
  
  # Tet I
  if (cur_speciesMT == "Tet_I") {
    
    prior_dist <- prior_dist_gen(3000, rtref_prior = c(0.2, 0.3), e_prior, eh_prior, th_prior)
    
    ss_4_priors <- c(prior(normal(log(0.2), 0.3), nlpar = "rtref"),
                     prior(normal(log(0.5), 0.3), nlpar = "e"),
                     prior(normal(log(3), 0.4), nlpar = "eh"),
                     prior(normal(35,3), nlpar = "th"))
  }
  
  # Par 69
  if (cur_speciesMT == "Par_69") {
    
    prior_dist <- prior_dist_gen(3000, rtref_prior = c(0.1, 0.5), e_prior, eh_prior, th_prior)
    
    ss_4_priors <- c(prior(normal(log(0.1), 0.5), nlpar = "rtref"),
                     prior(normal(log(0.5), 0.3), nlpar = "e"),
                     prior(normal(log(3), 0.4), nlpar = "eh"),
                     prior(normal(35,3), nlpar = "th"))
  }
  
  # Blep I
  if (cur_speciesMT == "Blep_I") {
    
    prior_dist <- prior_dist_gen(3000, rtref_prior = c(0.05, 0.4), e_prior, eh_prior, th_prior)
    
    ss_4_priors <- c(prior(normal(log(0.05), 0.4), nlpar = "rtref"),
                     prior(normal(log(0.5), 0.3), nlpar = "e"),
                     prior(normal(log(3), 0.4), nlpar = "eh"),
                     prior(normal(35,3), nlpar = "th"))
  }
  
  # Eug II
  if (cur_speciesMT == "Eug_II") {
    
    prior_dist <- prior_dist_gen(3000, rtref_prior = c(0.025, 0.3), e_prior, eh_prior, th_prior)
    
    ss_4_priors <- c(prior(normal(log(0.025), 0.3), nlpar = "rtref"),
                     prior(normal(log(0.5), 0.3), nlpar = "e"),
                     prior(normal(log(3), 0.4), nlpar = "eh"),
                     prior(normal(35,3), nlpar = "th"))
  }
  
  # fit model ---------------------------------------------------------------
  
  gc()
  
  fit_ss_4 <- brm(
    formula = ss_4,
    data = d,
    family = gaussian(),
    prior = ss_4_priors,
    control = list(adapt_delta = 0.99),
    chains = 1,
    iter = 3000
  )
  
  print(summary(fit_ss_4))
  
  stancode(fit_ss_4)
  
  print(paste("Finished fitting selection temp", cur_selection_temperature, "replicate", cur_replicate, "of", cur_speciesMT,  sep = " "))
  

# extract data ------------------------------------------------------------

  posteriors <- extract_posteriors(fit_ss_4)
  
  estimates <- extract_param_estimates(posteriors)
  
  predictions <- extract_predictions(fit_ss_4,0, 60, 200)
  

# plot data ---------------------------------------------------------------

  a <- plot_grid(
    plotlist = list(plot_niche(predictions),
                    plot_postprior(posteriors)),
    ncol = 1, nrow = 2)
  
  plot_grid(
    plotlist = list(a, plot_traceplot(fit_ss_4)),
    ncol = 2, nrow = 1
  )

}