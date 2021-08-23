# updated plotSpectra and plotDiet to handle multiple backgrounds

plotSpectra2 <- function (object, species = NULL, time_range, wlim = c(NA, NA), 
                          ylim = c(NA, NA), power = 1, biomass = TRUE, total = FALSE, 
                          resource = TRUE, background = TRUE, highlight = NULL, return_data = FALSE, 
                          ...) 
{
  if (missing(power)) {
    power <- as.numeric(biomass)
  }
  species <- valid_species_arg(object, species)
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    time_elements <- get_time_elements(object, time_range)
    n <- apply(object@n[time_elements, , , drop = FALSE], 
               c(2, 3), mean)
    n_pp <- apply(object@n_pp[time_elements, , drop = FALSE], 
                  2, mean)
    
    # if n_other has only one column, selecting time_elements lose the name for some reasons
    if(dim(object@n_other)[2] == 1)
    {
      n_other <- list()
      n_other[[dimnames(object@n_other)$component]] <- unlist(object@n_other[time_elements,])
    } else {
      n_other <- object@n_other[time_elements,]
    }
    
    ps <- plot_spectra2(object@params, n = n, n_pp = n_pp, n_other = n_other,
                        species = species, wlim = wlim, ylim = ylim, power = power, 
                        total = total, resource = resource, background = background, 
                        highlight = highlight, return_data = return_data)
    return(ps)
  }
  else if (is(object, "MizerParams")) {
    ps <- plot_spectra2(object, n = object@initial_n, n_pp = object@initial_n_pp, n_other = object@initial_n_other,
                        species = species, wlim = wlim, ylim = ylim, power = power, 
                        total = total, resource = resource, background = background, 
                        highlight = highlight, return_data = return_data)
    return(ps)
  }
  else {
    stop("First argument of `plotSpectra()` needs to be a MizerSim or ", 
         "a MizerParams object.")
  }
}


plot_spectra2 <- function(params, n, n_pp, n_other,
                          species, wlim, ylim, power,
                          total, resource, background, 
                          highlight, return_data) {
  
  
  params <- validParams(params)
  if (is.na(wlim[1])) {
    wlim[1] <- min(params@w) / 100
  }
  if (is.na(wlim[2])) {
    wlim[2] <- max(params@w_full)
  }
  
  if (total) {
    # Calculate total community abundance
    fish_idx <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    total_n <- n_pp
    total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
    total_n <- total_n * params@w_full^power
  }
  species <- valid_species_arg(params, species)
  # Deal with power argument
  if (power %in% c(0, 1, 2)) {
    y_label <- c("Number density [1/g]", "Biomass density",
                 "Biomass density [g]")[power + 1]
  } else {
    y_label <- paste0("Number density * w^", power)
  }
  n <- sweep(n, 2, params@w^power, "*")
  # Select only the desired species
  spec_n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  # Make data.frame for plot
  plot_dat <- data.frame(w = rep(params@w,
                                 each = dim(spec_n)[[1]]),
                         value = c(spec_n),
                         Species = dimnames(spec_n)[[1]],
                         Legend = dimnames(spec_n)[[1]])
  if (resource) {
    resource_sel <- (params@w_full >= wlim[1]) & 
      (params@w_full <= wlim[2])
    # Do we have any resource to plot?
    if (sum(resource_sel) > 0) {
      w_resource <- params@w_full[resource_sel]
      plank_n <- n_pp[resource_sel] * w_resource^power
      
      plot_dat <- rbind(plot_dat,
                        data.frame(w = w_resource,
                                   value = c(plank_n),
                                   Species = "Resource",
                                   Legend = "Resource"))
      
      if(length(n_other) > 0)
      {
        for(iBackground in 1:length(n_other))
        {
          other_resource <- params@w_full[resource_sel]
          other_n <- n_other[iBackground][[1]][resource_sel] * other_resource^power
          
          plot_dat <- rbind(plot_dat,
                            data.frame(w = w_resource,
                                       value = c(other_n),
                                       Species = names(n_other)[iBackground],
                                       Legend = names(n_other)[iBackground]))
        }
        
      }
      
      
    }
  }
  if (total) {
    plot_dat <- rbind(plot_dat,
                      data.frame(w = params@w_full,
                                 value = c(total_n),
                                 Species = "Total",
                                 Legend = "Total")
    )
  }
  if (background && anyNA(params@A)) {
    back_n <- n[is.na(params@A), , drop = FALSE]
    plot_dat <- 
      rbind(plot_dat,
            data.frame(w = rep(params@w,
                               each = dim(back_n)[[1]]),
                       value = c(back_n),
                       Species = as.factor(dimnames(back_n)[[1]]),
                       Legend = "Background")
      )
  }
  # lop off 0s and apply wlim
  plot_dat <- plot_dat[(plot_dat$value > 0) & 
                         (plot_dat$w >= wlim[1]) &
                         (plot_dat$w <= wlim[2]), ]
  # Impose ylim
  if (!is.na(ylim[2])) {
    plot_dat <- plot_dat[plot_dat$value <= ylim[2], ]
  }
  if (is.na(ylim[1])) {
    ylim[1] <- 1e-20
  }
  plot_dat <- plot_dat[plot_dat$value > ylim[1], ]
  
  if (return_data) return(plot_dat) 
  
  plotDataFrame(plot_dat, params, xlab = "Size [g]", ylab = y_label,
                xtrans = "log10", ytrans = "log10", 
                highlight = highlight, legend_var = "Legend")
}





getDiet2 <- function (params, n = initialN(params), n_pp = initialNResource(params), 
                      n_other = initialNOther(params), proportion = TRUE) 
{
  params <- validParams(params)
  species <- params@species_params$species
  no_sp <- length(species)
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)
  no_other <- length(params@other_encounter)
  other_names <- names(params@other_encounter)
  # assert_that(identical(dim(n), c(no_sp, no_w)), length(n_pp) == 
  # no_w_full)
  diet <- array(0, dim = c(no_sp, no_w, no_sp + 1 + no_other), 
                dimnames = list(predator = species, w = dimnames(params@initial_n)$w, 
                                prey = c(as.character(species), "Resource", other_names)))
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  if (length(params@ft_pred_kernel_e) == 1) {
    ae <- matrix(params@pred_kernel[, , idx_sp, drop = FALSE], 
                 ncol = no_w) %*% t(sweep(n, 2, params@w * params@dw, 
                                          "*"))
    diet[, , 1:no_sp] <- ae
    diet[, , no_sp + 1] <- rowSums(sweep(params@pred_kernel, 
                                         3, params@dw_full * params@w_full * n_pp, "*"), dims = 2)
  }
  else {
    prey <- matrix(0, nrow = no_sp + 1, ncol = no_w_full)
    prey[1:no_sp, idx_sp] <- sweep(n, 2, params@w * params@dw, 
                                   "*")
    prey[no_sp + 1, ] <- n_pp * params@w_full * params@dw_full
    ft <- array(rep(params@ft_pred_kernel_e, times = no_sp + 
                      1) * rep(mvfft(t(prey)), each = no_sp), dim = c(no_sp, 
                                                                      no_w_full, no_sp + 1))
    ft <- matrix(aperm(ft, c(2, 1, 3)), nrow = no_w_full)
    ae <- array(Re(mvfft(ft, inverse = TRUE)/no_w_full), 
                dim = c(no_w_full, no_sp, no_sp + 1))
    ae <- ae[idx_sp, , , drop = FALSE]
    ae <- aperm(ae, c(2, 1, 3))
    ae[ae < 1e-18] <- 0
    diet[, , 1:(no_sp + 1)] <- ae
  }
  inter <- cbind(params@interaction, params@species_params$interaction_resource)
  diet[, , 1:(no_sp + 1)] <- sweep(sweep(diet[, , 1:(no_sp + 
                                                       1), drop = FALSE], c(1, 3), inter, "*"), c(1, 2), params@search_vol, 
                                   "*")
  for (i in seq_along(params@other_encounter)) {
    diet[, , no_sp + 1 + i] <- do.call(params@other_encounter[[i]], 
                                       list(params = params, n = n, n_pp = n_pp, n_other = n_other, 
                                            component = names(params@other_encounter)[[i]]))
  }
  f <- getFeedingLevel(params, n, n_pp, n_other)
  fish_mask <- n > 0
  diet <- sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
  if (proportion) {
    total <- rowSums(diet, dims = 2)
    diet <- sweep(diet, c(1, 2), total, "/")
    diet[is.nan(diet)] <- 0
  }
  return(diet)
}

plotDiet2 <- function (sim, species = NULL, time_range, xlim = c(1, NA), returnData = F) 
{
  
  if (missing(time_range)) time_range <- max(as.numeric(dimnames(sim@n)$time))
  
  time_elements <- get_time_elements(sim, time_range)
  n <- apply(sim@n[time_elements, , , drop = FALSE], 
             c(2, 3), mean)
  n_pp <- apply(sim@n_pp[time_elements, , drop = FALSE], 
                2, mean)
  # if n_other has only one column, selecting time_elements lose the name for some reasons
  if(dim(sim@n_other)[2] == 1)
  {
    n_other <- list()
    n_other[[dimnames(sim@n_other)$component]] <- unlist(sim@n_other[time_elements,])
  } else {
    n_other <- sim@n_other[time_elements,]
  }
  
  params <- sim@params
  diet <- getDiet2(params,n = n, n_pp = n_pp, n_other = n_other)
  plot_dat <- melt(diet)
  plot_dat <- plot_dat[plot_dat$value > 0, ]
  colnames(plot_dat) <- c("Predator", "size", "Prey", "Proportion")
  if (is.null(species)) 
    p <- ggplot(plot_dat) + facet_wrap(. ~ Predator, scales = "free")
  else p <- ggplot(filter(plot_dat, Predator == species))
  p <- p + geom_area(aes(x = size, y = Proportion, fill = Prey)) + 
    scale_x_continuous(limits = c(1, NA), name = "Size [g]", 
                       trans = "log10") + scale_fill_manual(values = sim@params@linecolour) + 
    theme(legend.position = "right", legend.key = element_rect(fill = "black"), 
          panel.background = element_blank(), panel.grid.minor = element_line(color = "gray"), 
          strip.background = element_blank())
  if (returnData) 
    return(plot_dat)
  else return(p)
}


getDietComp<- function(sim)
{
  # initialisation
  object <- sim@params
  pred_kernel <- getPredKernel(object)
  n = sim@n[dim(sim@n)[1],,]
  n_pp = sim@n_pp[dim(sim@n_pp)[1],]
  no_sp <- dim(object@species_params)[1]
  no_w <- length(object@w)
  no_w_full <- length(object@w_full)
  
  
  no_other <- length(object@other_encounter)
  other_names <- names(object@other_encounter)

  feedinglevel <- getFeedingLevel(object, n, n_pp, n_other = sim@n_other[dim(sim@n_other)[1],])

  # diet_comp <- array(0, dim = c(no_sp, no_w, no_sp + 1 + no_other, no_w_full), 
  #               dimnames = list(predator = as.character(object@species_params$species), 
  #                               w = dimnames(object@initial_n)$w, 
  #                               prey = c(as.character(object@species_params$species), "Resource", other_names),
  #                               prey_size = object@w_full))
  
  diet_comp<-array(0, c(no_sp, no_w, no_sp + 1, no_w_full),
                   dimnames=list( predator=as.character(object@species_params$species), pred_size = object@w,
                                  prey = c(as.character(object@species_params$species), "background"),
                                  prey_size = object@w_full))
  
  # Biomass by species
  n_total_in_size_bins<- sweep(n, 2, object@dw , "*")
  b_tot <- sweep(n_total_in_size_bins, 2, object@w , "*")
  
  # Index of predator size classes
  idx_sp<- object@w_full %in% object@w
  
  inter <- cbind(params@interaction, params@species_params$interaction_resource)
  
  
  #  pred_kernel * interaction matrix
  for(iW in 1:no_w){
    for(iSpecies in 1:no_sp){
      diet_comp[iSpecies,iW,,idx_sp]<- sweep(sweep( b_tot, c(1), inter[iSpecies, ], "*"), c(2),
                                                    pred_kernel[iSpecies,iW,idx_sp], "*")
    }
  }
  # Search rate *  feeding level * prey biomass
  diet_comp<- sweep(sweep(sweep(diet_comp, c(1,2), object@search_vol,"*"),
                                      c(1,2),1-feedinglevel,"*"),
                                c(1,2),b_tot,"*")  # Prey eaten: total g prey/ year  (given predator biomass density)
  
  # no interaction matrix for background spectrum
  # b_background <- (sweep(pred_kernel[,,], c(3), object@dw_full*object@w_full*n_pp, "*"))
  # #Search rate *  feeding level * predator biomass
  # b_background<- sweep(b_background, c(1,2), object@search_vol,"*") #Scale up by search volume
  # b_background<- sweep(b_background, c(1,2), feedinglevel,"*") # Scale according to feeding level. Prey eaten: g prey / year / g predator
  # b_background_tot<-sweep(b_background,c(1,2), b_tot, "*") # Prey eaten: total g prey/ year  (given predator biomass density)
  # 
  # # Store background eaten
  # diet_comp[,,no_sp+1,]<- b_background_tot
  
  return(diet_comp)
}

getDietComp2<- function(sim)
{
  # initialisation
  object <- sim@params
  pred_kernel <- getPredKernel(object)
  n = sim@n[dim(sim@n)[1],,]
  n_pp = sim@n_pp[dim(sim@n_pp)[1],]
  no_sp <- dim(object@species_params)[1]
  no_w <- length(object@w)
  no_w_full <- length(object@w_full)
  
  
  no_other <- length(object@other_encounter)
  other_names <- names(object@other_encounter)
  
  feedinglevel <- getFeedingLevel(object, n, n_pp, n_other = sim@n_other[dim(sim@n_other)[1],])
  
  diet_comp <- array(0, dim = c(no_sp, no_w, no_sp + 1 + no_other, no_w_full), 
                     dimnames = list(predator = as.character(object@species_params$species), 
                                     w = dimnames(object@initial_n)$w, 
                                     prey = c(as.character(object@species_params$species), "Resource", other_names),
                                     prey_size = object@w_full))
  
  
  # Biomass by species
  n_total_in_size_bins<- sweep(n, 2, object@dw , "*")
  b_tot <- sweep(n_total_in_size_bins, 2, object@w , "*")
  
  # Index of predator size classes
  idx_sp<- object@w_full %in% object@w
  
  # inter <- cbind(params@interaction, params@species_params$interaction_resource)
  
  
  #  pred_kernel * interaction matrix
  for(iW in 1:no_w){
    for(iSpecies in 1:no_sp){
      diet_comp[iSpecies,iW,1:no_sp,idx_sp]<- sweep(sweep( b_tot, c(1), object@interaction[iSpecies, 1:no_sp], "*"), c(2),
                                                    pred_kernel[iSpecies,iW,idx_sp], "*")
    }
  }
  # Search rate *  feeding level * prey biomass
  diet_comp[,,1:no_sp,]<- sweep(sweep(sweep(diet_comp[,,1:no_sp,], c(1,2), object@search_vol,"*"),
                                      c(1,2),1-feedinglevel,"*"),
                                c(1,2),b_tot,"*")  # Prey eaten: total g prey/ year  (given predator biomass density)
  
  # no interaction matrix for background spectrum
  b_background <- (sweep(pred_kernel[,,], c(3), object@dw_full*object@w_full*n_pp, "*"))
  #Search rate *  feeding level * predator biomass
  b_background<- sweep(b_background, c(1,2), object@search_vol,"*") #Scale up by search volume
  b_background<- sweep(b_background, c(1,2), feedinglevel,"*") # Scale according to feeding level. Prey eaten: g prey / year / g predator
  b_background_tot<-sweep(b_background,c(1,2), b_tot, "*") # Prey eaten: total g prey/ year  (given predator biomass density)
  
  # Store background eaten
  diet_comp[,,no_sp+1,]<- b_background_tot
  
  return(diet_comp)
}