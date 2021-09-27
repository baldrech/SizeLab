
# In this file we are setting up the functions for
# - the dynamics of the resources
# - the contribution of the resources to the encounter rate
# and provide a function for setting model with
# benthos and algae.
#
# JB: this set up is based on Gustav's set-up for Asta's model, editing here to test as we wnat to use multiple resource spectra in a bunch of projects

library(tidyverse)
library(mizerExperimental)


background_semichemostat <- function(params, n_other, rates, dt, component,
                                     ...) {
  c <- params@other_params[[component]]
  # name of interaction parameter for this component in species_params
  interaction_component <- paste0("interaction_", component)
  interaction <- params@species_params[[interaction_component]]
  mort <- as.vector(interaction  %*% rates$pred_rate)
  tmp <- c$rate * c$capacity / (c$rate + mort)
  return(tmp - (tmp - n_other[[component]]) * exp(-(c$rate + mort) * dt))
}

background_encounter <- function(params, n, n_pp, n_other, ...) {
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  prey <- outer(params@species_params$interaction_resource, n_pp) +
  #  outer(params@species_params$interaction_aa, n_other$aa) +
  #  outer(params@species_params$interaction_bb, n_other$bb)
    outer(params@species_params$interaction_aa, n_other[[1]]) +
    outer(params@species_params$interaction_bb, n_other[[2]])
  prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
  prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
  avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                     mvfft(base::t(prey)),
                                   inverse = TRUE))) / length(params@w_full)
  avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  avail_energy[avail_energy < 1e-18] <- 0
  
  params@search_vol * avail_energy
}


addResources <- function(params) {

  
  # Add macroalgae 
  kappa <- 16
  lambda <- 1.6
  r <- 2
  n <- 2/3
  max <- 50
  min <- 1e-3
  rate <- r * params@w_full^(n - 1)
  capacity <- kappa * params@w_full^(-lambda)
  capacity[params@w_full > max] <- 0
  capacity[params@w_full < min] <- 0
  params <- setComponent(params = params, component = "aa",
                         initial_value = params@initial_n_pp,
                         dynamics_fun =  "background_semichemostat",
                         component_params = list(rate = rate,
                                                 capacity = capacity))
  
  # Add benthic resource
  kappa <- kappa_ben
  lambda <- lambda_ben
  r <- 1
  n <- 2/3
  max <- 5
  min <- 1e-3
  rate <- r * params@w_full^(n - 1)
  capacity <- kappa * params@w_full^(-lambda)
  capacity[params@w_full > max] <- 0
  capacity[params@w_full < min] <- 0
  params <- setComponent(params = params, component = "bb",
                         initial_value = params@initial_n_pp,
                         dynamics_fun =  "background_semichemostat",
                         component_params = list(rate = rate,
                                                 capacity = capacity))
  
  # Include these extra resources in the encounter rate
  params <- setRateFunction(params, "Encounter", "background_encounter")
  
 
  params
}


## test with northsea params

params <- NS_params

species_params(params)$interaction_resource<-0.5
species_params(params)$interaction_aa<-rep(0.1,12)
species_params(params)$interaction_bb<-rep(0.1,12)
  
  
newparams<- addResources(params)


newsim<- project(newparams, t_max = 10, effort = 0, dt = dt)

# they don't seem to be there - where do you actually set up the structure of these n_other? 


# instead they are lists but you have to refer to them as so, not being picked up in code above?
newsim@n_other[[1]]
newsim@n_other[[2]]


# edited lines 28 and 29 above BUT only seems to be a single time step 
# is it assumed to be fixed through time?
# vectors are same length as:
length(params@w_full)
# but should be time x w_full - we want to keep track of these resources!

# would it not be easier/neater to simply add them to the n_pp array, when they are size-structured??

newsim<- project(newparams, t_max = 10, effort = 0, dt = dt)
plot(newsim)


