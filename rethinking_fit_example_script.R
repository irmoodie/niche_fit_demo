# Sharpe-Schoolfield 4p fit with rethinking
# IR Moodie

# clear environment -------------------------------------------------------

rm(list = ls())


# load packages -----------------------------------------------------------

library(rethinking)
library(tidybayes)
library(tidybayes.rethinking)
library(cowplot)
library(tidyverse)

# options -----------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# load in required datasets -----------------------------------------------

data <- read_csv("exampledata.csv")


# build model with ulam ---------------------------------------------------

datalist <- list(
  r0 = data$r0_median,
  r0_sd = data$r0_sd,
  temperature = data$temperature,
  k = data$k,
  tref = data$tref,
  N = nrow(data)
)

model <- ulam(
  alist(
    r0 ~ dnorm(r0_true, r0_sd),
    vector[N]:r0_true ~ dnorm(est, sigma),
    est <- (exp(rtref) * exp(exp(e)/k * (1/tref - 1/(temperature + 273.15))))*(1/(1 + exp(exp(eh)/k * (1/(th + 273.15) - 1/(temperature + 273.15))))),
    rtref ~ dnorm(log(0.2), 0.3),
    th ~ dnorm(35,3),
    e ~ dnorm(log(0.5), 0.3),
    eh ~ dnorm(log(3), 0.4),
    sigma ~ dexp(1)
  ), data=datalist , iter = 3000, warmup = 1000)


# show stan code ----------------------------------------------------------

stancode(model)


# functions for topt and rmax ---------------------------------------------

get_topt <- function(eh, th, e){
  return((eh * (th+273.15))/(eh + (8.62e-05 * (th+273.15) * log((eh/e) - 1)))-273.15)
}

ss_4_formula <- function (temp, rtref, e, eh, th, tref) 
{
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- rtref * exp(e/k * (1/tref - 1/(temp + 
                                                      273.15)))
  inactivation.term <- 1/(1 + exp(eh/k * (1/(th + 273.15) - 
                                            1/(temp + 273.15))))
  return(boltzmann.term * inactivation.term)
}


# get draws and calc topt and rmax for each draw --------------------------

posteriors <-
  model %>%
  spread_draws(rtref, th, e, eh, sigma) %>%
  mutate(rtref = exp(rtref),
         e = exp(e),
         eh = exp(eh),
         topt = get_topt(eh, th, e),
         rmax = ss_4_formula(topt, rtref, e, eh, th, 20)) %>%
  gather('param', 'estimate', 4:ncol(.))

estimates <-
  posteriors %>%
  group_by(param) %>%
  mean_qi()


# predict from draws to get mean + error ----------------------------------

predictions <-
  tibble(temperature = seq(0, 50, length.out = 200),
           tref = 20 + 273.15,
           k = 8.62e-05) %>%
  add_fitted_draws(model) %>%
  data.frame() %>%
  group_by(temperature) %>%
  mean_qi(estimate = .value)


# generate prior distributions --------------------------------------------

priors <-
  tibble(
  rtref = exp(rnorm(1:2000,log(0.2), 0.3)),
  e = exp(rnorm(1:2000,log(0.5), 0.3)),
  eh = exp(rnorm(1:2000,log(3), 0.4)),
  th = rnorm(1:2000,35,3)) %>%
  gather('param', 'estimate', 1:ncol(.)) %>%
  mutate(type = "prior")


# plot priors and posteriors ----------------------------------------------

postpriorplot <-
  posteriors %>%
  filter(estimate < 50 & param != "sigma") %>%
  group_by(param) %>%
  mutate(median = median(estimate),
         type = "post") %>%
  bind_rows(priors) %>%
  ggplot(aes(x = estimate)) +
  facet_wrap(~param, scales = "free") +
  geom_density(aes(fill = type), alpha = 0.3) +
  geom_vline(aes(xintercept = median), linetype = "dashed") +
  scale_fill_manual(values = c("red", "grey")) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# traceplot if interested -------------------------------------------------

traceplot <- traceplot_ulam(model)


# plot niche function from predictions ------------------------------------

niche <-
  ggplot(predictions, aes(x = temperature, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_point(data = data, aes(x = temperature, y = r0_mean), size  = 3) +
  geom_errorbar(data = data, aes(x = temperature, y = r0_mean, ymin = r0_min, ymax = r0_max)) +
  labs(x = "Temperature",
       y = "Intrinsic rate of increase, r0") +
  theme_cowplot() +
  background_grid(major = "xy", minor = "x")


# plot niche and post/prior dists -----------------------------------------

plot_grid(
  plotlist = list(niche,
                  postpriorplot),
  ncol = 2, nrow = 1)
