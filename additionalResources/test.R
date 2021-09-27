# devtools::install_github("sizespectrum/mizerMR")
library(mizerMR)
library(tibble)
source("additionalResources/utility.R")

resource_params <- tribble(
  ~resource,  ~kappa, ~lambda, ~r_pp, ~w_min, ~w_max,
  "Resource 1",    1e11,    2.13,     4,    NA ,   0.01,
  "Resource 2",    1e11,    2.05,    10,   1e-4,     NA
)

params <- setMultipleResources(NS_params, resource_params)


sim <- project(params, t_max = 2, t_save = 0.2)


plotSpectra(sim)

# editing the resource on the fly

# resource_params(params)$kappa[1] <- 1e18
# 
# # it works
# 
# 
# # own encounter function
# 
# mizerEncounter2 <- function (params, n, n_pp, n_other, t, ...) 
# {
#   idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
#   if (!is.null(comment(params@pred_kernel))) {
#     n_eff_prey <- sweep(params@interaction %*% n, 2, params@w * 
#                           params@dw, "*", check.margin = FALSE)
#     phi_prey_species <- rowSums(sweep(params@pred_kernel[, 
#                                                          , idx_sp, drop = FALSE], c(1, 3), n_eff_prey, "*", 
#                                       check.margin = FALSE), dims = 2)
#     phi_prey_background <- params@species_params$interaction_resource * 
#       rowSums(sweep(params@pred_kernel, 3, params@dw_full * 
#                       params@w_full * n_pp, "*", check.margin = FALSE), 
#               dims = 2)
#     encounter <- params@search_vol * (phi_prey_species + 
#                                         phi_prey_background)
#   }
#   else {
#     prey <- outer(params@species_params$interaction_resource, 
#                   n_pp)
#     prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% 
#       n
#     prey <- sweep(prey, 2, params@w_full * params@dw_full, 
#                   "*")
#     avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
#                                        mvfft(base::t(prey)), inverse = TRUE)))/length(params@w_full)
#     avail_energy <- avail_energy[, idx_sp, drop = FALSE]
#     avail_energy[avail_energy < 1e-18] <- 0
#     encounter <- params@search_vol * avail_energy
#   }
#   for (i in seq_along(params@other_encounter)) {
#     encounter <- encounter + do.call(params@other_encounter[[i]], 
#                                      list(params = params, n = n, n_pp = n_pp, n_other = n_other, 
#                                           component = names(params@other_encounter)[[i]], 
#                                           ...))
#   }
#   return(encounter)
# }
# 
# params <- setRateFunction(params, "Encounter", "mizerEncounter2")
# params <- setParams(params)
# sim <- project(params)


### testing diet

plotDiet(sim@params)
plotDiet2(sim)
getDiet2(sim@params)
