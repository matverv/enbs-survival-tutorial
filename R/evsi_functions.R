# Copyright (c) 2023 Mathyn Vervaart, Eline Aas, Karl Claxton, Anna Heath, Mark Strong, Nicky Welton, Torbjørn Wisløff 
# Licensed under the MIT License

#####################################################################################
# Install and load packages 
#####################################################################################
list_of_packages <- c("here", "MASS", "ggplot2", "mgcv", "drc", "earth")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, library, character.only = TRUE)


##############################################################################################
# Function for computing EVSI for overall survival
##############################################################################################
evsi_os_fun <- function (m_nb, l_os, l_atrisk_times_os, max_add_fu, ncyc_y, l_dropout=NULL, l_enroll=NULL, ...) {
  
  # stop conditions
  if(!is.numeric(m_nb)) {stop("The net benefits is not a numeric vector.")}
  if(!is.list(l_os)) {stop("Overall survival matrices are not stored in a list.")}
  if(!is.list(l_atrisk_times_os)) {stop("Vectors of observed follow-up times for patients at risk for overall survival are not stored in a list.")}
  if(!is.numeric(max_add_fu)) {stop("Maximum additional follow-up time (in months) is not numeric.")}
  if(max_add_fu < 6) {stop("Maximum additional follow-up time (in months) should be at least 6.")}
  if(!is.null(l_enroll) & mean(unlist(l_atrisk_times_os))>0) {stop("The list of enrollment rates should be NULL when the start times > 0.")}
  
  surv_evsi_fun(
    m_nb = m_nb,                            # net benefits
    l_os = l_os,                            # overall survival (OS) probabilities over discrete model cycles for each treatment 
    l_atrisk_times_os = l_atrisk_times_os,  # individual observed follow-up times for OS for each treatment
    add_fu = add_fu_fun(max_add_fu),        # additional follow-up times in months
    ncyc_y = ncyc_y,                        # number of model cycles per year
    l_dropout = l_dropout,                  # random censoring (dropout) rates for each treatment
    l_enroll = l_enroll,                    # enrollment rates for each treatment (optional)
    ...
  )
}


##############################################################################################
# Function for computing EVSI for overall survival and progression-free survival
##############################################################################################
evsi_os_pfs_fun <- function (m_nb, l_os, l_pfs, l_atrisk_times_os, l_atrisk_times_pfs, add_fu, ncyc_y, l_dropout=NULL, l_enroll=NULL, ...) {
  
  # stop conditions
  if(!is.numeric(m_nb)) {stop("The net benefits is not numeric.")}
  if(!is.list(l_os)) {stop("Overall survival matrices are not stored in a list.")}
  if(!is.list(l_pfs)) {stop("Progression-free survival matrices are not stored in a list.")}
  if(!is.list(l_atrisk_times_os)) {stop("Vectors of observed follow-up times for patients at risk for overall survival are not stored in a list.")}
  if(!is.list(l_atrisk_times_pfs)) {stop("Vectors of observed follow-up times for patients at risk for progression-free survival are not stored in a list.")}
  if(!is.numeric(max_add_fu)) {stop("Maximum additional follow-up time (in months) is not numeric.")}
  if(max_add_fu < 6) {stop("Maximum additional follow-up time (in months) should be at least 6.")}
  if(!is.null(l_enroll) & mean(unlist(l_atrisk_times_os))>0) {stop("The list of enrollment rates should be NULL when the start times > 0.")}
  
  surv_evsi_fun(
    m_nb = m_nb,                                  # net benefits
    l_os = l_os,                                    # overall survival (OS) probabilities over discrete model cycles for each treatment 
    l_pfs = l_pfs,                                  # progression-free survival (PFS) probabilities over discrete model cycles for each treatment 
    l_atrisk_times_os = l_atrisk_times_os,          # individual observed follow-up times for OS for each treatment
    l_atrisk_times_pfs = l_atrisk_times_pfs,        # individual observed follow-up times for PFS for each treatment
    add_fu = add_fu_fun(max_add_fu),                # additional follow-up times in months
    ncyc_y = ncyc_y,                                # number of model cycles per year
    l_dropout = l_dropout,                          # random censoring (dropout) rates for each treatment
    l_enroll = l_enroll,                            # enrollment rates for each treatment (optional)
    ...
  )
}

##############################################################################################
# Function for specifying a vector of evenly spaced additional follow-up times
##############################################################################################
add_fu_fun <- function (max_add_fu) {
  
  add_fup <- c(2, max_add_fu)
  fu_range <- max(add_fup) - min(add_fup)
  add_fup <- c(min(add_fup), 
              min(add_fup) + (1/3) * fu_range,
              min(add_fup) + (2/3) * fu_range,
              max(add_fup))
  add_fup <- round(add_fup)
  return(add_fup)
}


##############################################################################################
# function for generating survival data and computing EVSI for OS and/or PFS
##############################################################################################
surv_evsi_fun <- function (
    m_nb,                           # net benefits
    l_os=NULL,                      # overall survival probabilities over discrete model cycles for each treatment 
    l_pfs=NULL,                     # progression-free survival probabilities over discrete model cycles for each treatment 
    l_atrisk_times_os=NULL,         # individual observed follow-up times for OS for each treatment
    l_atrisk_times_pfs=NULL,        # individual observed follow-up times for PFS for each treatment
    add_fu,                         # additional follow-up times in months
    ncyc_y,                         # number of model cycles per year
    l_dropout = NULL,               # random censoring (dropout) rates for each treatment (optional)
    l_enroll = NULL,                # enrollment rates for each treatment (optional)
    data_method = "interpolation",  # "discrete" or "interpolation"
    evsi_method = "gam",            # "gam" or "mars"
    fast = TRUE,
    seed = NULL)                    # only applicable for discrete method
{
  
  #if(!exists("seed")){seed <- 1}
  if(is.null(seed)){seed <- 1}
 
  cycle2month <- 12 / ncyc_y

  m_evsi <- sapply(add_fu / cycle2month, function (t) { 
    
    ######### Overall survival only #########
    if(is.null(l_pfs)) {
      
      arm_indic <- rep(1:length(l_os), 1, each = 2)              # trial arm indicator for each survival curve
      arms <- unique(arm_indic)                                  # number of trial arms
      l_atrisk_times_os <- lapply(l_atrisk_times_os, function (x) {x / cycle2month}) # convert start times to model time unit
      
      # generate survival data
      if(data_method == "interpolation") {
        summ_stat <- interpol_os_data_fun(l_os, l_atrisk_times_os, l_dropout, l_enroll, t, seed, cycle2month) 
      } 
     
      if(data_method == "discrete") {
        summ_stat <- discrete_os_data_fun(l_os, l_atrisk_times_os, l_dropout, l_enroll, t, seed, fast = fast, cycle2month)
      } 
      
    } 
    
    ######### Overall survival and progression-free survival #########
    if(is.null(l_os) == FALSE & is.null(l_pfs) == FALSE) { #& is.null(l_atrisk_times_os) == FALSE
      
      arm_indic <- rep(1:length(l_os), 2, each = 2)              # trial arm indicator for each survival curve
      arms <- unique(arm_indic)                                       # number of trial arms
      l_atrisk_times_os <- lapply(l_atrisk_times_os, function (x) {x / cycle2month}) # convert start times to model time unit
      l_atrisk_times_pfs <- lapply(l_atrisk_times_pfs, function (x) {x / cycle2month}) # convert start times to model time unit
      
      # generate survival data
      if(data_method == "interpolation") {
        summ_stat <- interpol_os_pfs_data_fun(l_os, l_pfs, l_atrisk_times_os, l_atrisk_times_pfs, l_dropout, l_enroll, t, seed, cycle2month) 
      } 
      
      if(data_method == "discrete") {
        summ_stat <- discrete_os_pfs_data_fun(l_os, l_pfs, l_atrisk_times_os, l_atrisk_times_pfs, l_dropout, l_enroll, t, seed, fast = fast, cycle2month)
      } 
   
    }
 
    ######### Compute EVSI using GAM #########
    if(evsi_method == "gam") {
      
      regr_model <- reg_mod_fun(summ_stat, arms, arm_indic)
      evsi <- gam_evsi_fun(m_nb, summ_stat, regr_model, arm_indic)
      
      # output for OS only
      if(is.null(l_pfs)) {
        
        # generate the sequence of column positions for the treatment arms
        arm_cols <- seq(from = 1, to = length(arms) * 2, by = 2)
        
        # extract the os_events and os_events_se for all treatment arms
        os_events <- sum(colMeans(summ_stat[, arm_cols]))
        os_events_se <- sd(colMeans(summ_stat[, arm_cols]))
        
        print(c("add. follow-up" = round(t*cycle2month,4), "evsi" = round(evsi[[1]],3), "se" = round(evsi[[2]],3), "lower" = round(evsi[[3]],3), "upper" = round(evsi[[4]],3), 
                "os_events" = round(os_events, 0),  "os_events_se" = round(os_events_se, 3)  ))
        
        return(c("add. follow-up" = round(t*cycle2month,8), "evsi" = round(evsi[[1]],8), "se" = round(evsi[[2]],8), "lower" = round(evsi[[3]],8), "upper" = round(evsi[[4]],8), 
                 "os_events" = round(os_events, 0),  "os_events_se" = round(os_events_se, 8) )) 
        
      }
      
      # output for OS + PFS 
      if(is.null(l_os) == FALSE & is.null(l_pfs) == FALSE) {
        
        # generate the sequence of column positions for the treatment arms
        arm_cols <- seq(from = 1, to = length(arms) * 2, by = 2)
        
        # extract the os_events and os_events_se for all treatment arms
        os_events <- sum(colMeans(summ_stat[, arm_cols]))
        os_events_se <- sd(colMeans(summ_stat[, arm_cols]))
        
        # extract the pfs_events and pfs_events_se for all treatment arms
        pfs_events <- sum(colMeans(summ_stat[, arm_cols + length(arms) * 2]))
        pfs_events_se <- sd(colMeans(summ_stat[, arm_cols + length(arms) * 2]))
        
        print(c("add. follow-up" = round(t*cycle2month,3), "evsi" = round(evsi[[1]],3), "se" = round(evsi[[2]],3), "lower" = round(evsi[[3]],3), "upper" = round(evsi[[4]],3), 
                "os_events" = round(os_events, 0),  "os_events_se" = round(os_events_se, 3),
                "pfs_events" = round(pfs_events, 0),  "pfs_events_se" = round(pfs_events_se, 3) ))
        
        return(c("add. follow-up" = round(t*cycle2month,6), "evsi" = round(evsi[[1]],8), "se" = round(evsi[[2]],8), "lower" = round(evsi[[3]],8), "upper" = round(evsi[[4]],8), 
                 "os_events" = round(os_events, 0),  "os_events_se" = round(os_events_se, 8), 
                 "pfs_events" = round(pfs_events, 0),  "pfs_events_se" = round(pfs_events_se, 8) ))

        }
      
    } # close gam calculation
    
    
    ######### Compute EVSI using MARS #########
    if(evsi_method == "mars") {
      
      evsi <- mars_evsi_fun(m_nb, summ_stat, arm_indic)
      
      # output for OS only
      if(is.null(l_pfs)) {
        print(c("add. follow-up" = round(t*cycle2month, 0), "evsi" = round(evsi[[1]],3), 
                "os_events" = round(mean(summ_stat[,1]) + mean(summ_stat[,3]), 0) ))
        
        return(c("add. follow-up" = round(t*cycle2month, 0), "evsi" = round(evsi[[1]],3), 
                 "os_events" = round(mean(summ_stat[,1]) + mean(summ_stat[,3]), 0) ))
      }
      
      # output for OS + PFS 
      if(is.null(l_os) == FALSE & is.null(l_pfs) == FALSE) {
        
        print(c("add. follow-up" = round(t*cycle2month, 0), "evsi" = round(evsi[[1]],3), 
                "os_events" = round(mean(summ_stat[,1]) + mean(summ_stat[,3]), 0), "pfs_events" = round(mean(summ_stat[,5]) + mean(summ_stat[,7]), 0) ))
        
        return(c("add. follow-up" = round(t*cycle2month, 6), "evsi" = round(evsi[[1]],6), 
                 "os_events" = round(mean(summ_stat[,1]) + mean(summ_stat[,3]), 0), "pfs_events" = round(mean(summ_stat[,5]) + mean(summ_stat[,7]), 0) ))
      }
      
    } # close mars calculation
    
  })
  
  # interpolate the EVSI using asymptotic regression
  df_evsi <- evsi_ar_fun(m_evsi, add_fu)
  df_evsi <- subset(df_evsi, select = -c(group))
  return(df_evsi)
}


##############################################################################################
# Function to define a GAM regression model
##############################################################################################
reg_mod_fun <- function (summ_stat, arms, arm_indic) {
  
  var_names <- lapply(arms, function (x) {
    colnames(summ_stat[arm_indic == x])
  })
  
  regr_model <- lapply(var_names, function (x) {
    paste("te(", paste(x,collapse = ","), ",k=4)", sep = "")
  })
  
  regr_model <- lapply(2:length(regr_model), function (x) {
    paste(unlist(c(regr_model[[1]], regr_model[[x]])), collapse = "+")
  })
  
  return(regr_model)
}


##############################################################################################
# Function for computing EVSI using GAM
##############################################################################################
gam_evsi_fun <- function (m_nb, summ_stat, regr_model, arm_indic) {
  
  print(paste("Estimating posterior INB"))
    
  # matrix of incremental net benefits
  m_inb <- matrix(m_nb[,2:ncol(m_nb)] - m_nb[,1], ncol = (ncol(m_nb)-1) )

  # number of GAM parameter samples for computation of SE and CI
  nsim_gam <- 2000 
  
  ### GAM regression
  l_gam_est <- lapply(1:length(regr_model), function (i) { 

    f <- update(formula(m_inb[,i] ~ .), formula(paste(".~", regr_model[[i]])))
    #mod <- gam(f, data = data.frame(summ_stat)) #summ_stat[,arm_indic==i]
    mod <- bam(f, data = data.frame(summ_stat), discrete=F) #summ_stat[,arm_indic==i]
    
    # extract fitted values                  
    g_hat <- mod$fitted

    ### compute standard errors

    # extract the basis function values
    Xstar <- model.matrix(mod)
    
    # extract coefficients
    beta <- mod$coef
    
    # covariance matrix
    v_cov <- vcov(mod)
    
    # sample from the parameter distributions
    set.seed(123)
    parameter_draws <- t(mvrnorm(nsim_gam, beta, v_cov))
    
    # # from these draws, calculate draws from fitted values
    fitted_draws <- Xstar %*% parameter_draws
    
    return(list(g_hat = g_hat, fitted_draws = fitted_draws))
  })
  

  # extract g_hat and fitted draws
  g_hat <- sapply(l_gam_est, function (x) x[[1]])
  g_hat <- cbind(0, g_hat)
  fitted_draws <- c(list(matrix(0, nrow = nrow(g_hat), ncol = nsim_gam)),
                    lapply(l_gam_est, function (x) x[[2]]) )
  rm(l_gam_est)

  ### compute evsi
  v_inb_max <- apply(g_hat, 1, function (x) max(x)) # maximum INB with sample information
  evsi <- mean(v_inb_max) - max(colMeans(g_hat)) 
  
  ### compute standard error

  # compute EVSI for each sample
  v_nb_max <- colMeans(do.call(pmax, fitted_draws)) # maximum INB with sample information
  v_max_means <- lapply(fitted_draws, colMeans)
  v_max_means <- do.call(pmax, v_max_means) 
  v_evsi <- v_nb_max - v_max_means

  # standard error and ci
  se <- sd(v_evsi)
  ci <- quantile(v_evsi, c(0.025, 0.975))
  #ci <- c(pmax(evsi - se * qnorm(0.975), 0), evsi + se * qnorm(0.975))
  
  return(list("evsi"= round(evsi,6), "se"=round(se,6), "lower" = round(ci[1],6), "upper" = round(ci[2],6) ))
}



##############################################################################################
# Function for computing EVSI using MARS
##############################################################################################
mars_evsi_fun <- function (m_nb, summ_stat, arm_indic) {
  
  print(paste("Estimating posterior NB.."))
  
  m_mars_post <- sapply(1:ncol(m_nb), function (i) {
    mars_mod <- earth::earth(m_nb[,i] ~ ., data = data.frame(summ_stat[,arm_indic==i, drop = FALSE]),
                             nk = 1e3, fast.k = 20, degree = 10, thresh = 1e-4)
    posterior <- mars_mod$fitted.values
    return(posterior)
  })
  
  ### compute evsi
  v_nb_max <- apply(m_mars_post, 1, function (x) max(x)) # compute the pairwise maximum
  evsi <- mean(v_nb_max) - max(colMeans(m_mars_post)) 
  
  return(evsi)
  
}


##############################################################################################
# Function for computing partial EVPI using GAM
##############################################################################################
surv_pevpi_fun <- function (m_nb, l_os=NULL, l_pfs=NULL) {
  
  ######### Overall survival only #########
  if(is.null(l_pfs)) {
    summ_stat <- data.frame(sapply(l_os, function (x) {
      colSums(x)
    }))
    arm_indic <- rep(1:length(l_os), 1, each = 1)         # trial arm indicator
    arms <- unique(arm_indic)
  }
  
  ######### Progression-free survival only #########
  if(is.null(l_os)) {
    summ_stat <- data.frame(sapply(l_pfs, function (x) {
      colSums(x)
    }))
    arm_indic <- rep(1:length(l_pfs), 1, each = 1)         # trial arm indicator
    arms <- unique(arm_indic)
  }
  
  ######### Overall survival and progression-free survival #########
  if(is.null(l_os) == FALSE & is.null(l_pfs) == FALSE) {
    curves <- c(l_os, l_pfs)
    summ_stat <- data.frame(sapply(curves, function (x) {
      colSums(x)
    }))
    arm_indic <- rep(1:length(l_os), 2, each = 1)         # trial arm indicator
    arms <- unique(arm_indic)
  }
  
  ######### Compute partial EVPI using GAM #########
  regr_model <- reg_mod_fun(summ_stat, arms, arm_indic)
  pevpi <- gam_evsi_fun(m_nb, summ_stat, regr_model, arm_indic)
  names(pevpi) <- c("pevpi", "se", "lower", "upper")
  
  ######### Compute partial EVPI using MARS #########
  #pevpi <- mars_evsi_fun(m_nb, summ_stat)
  
  print(round(unlist(pevpi),3))
  
  return(pevpi)
}


#####################################################################################
# Function for removing evenly spaced out elements from an ordered vector
#####################################################################################
split_fun <- function (times, n) {
  
  # moving average function to compute the values between the quantiles
  ma <- function(x, n = 2) as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))[-length(x)]
  
  quantile_values <- ma(seq(0, 1, length.out = n + 1))
  q <- quantile(times, probs = quantile_values)
  values <- as.numeric(times[findInterval(q, times) ] )
  index <- match(values, times)
  return(times[-index])
  
}


#####################################################################################
# Function for finding the intersect between two vectors
#####################################################################################
intersect_fun <- function (x, y, x2, y2, step) {
  v1 <- approx(x, y, xout = seq.int(min(x), max(x), step))
  v2 <- approx(x2, y2, xout = seq.int(min(x), max(x), step))
  
  intersect <- v1$x[which.min(abs(v1$y - v2$y))]
  return(intersect)
}


####################################################################################
# function for interpolating the EVSI estimates using asymptotic regression
####################################################################################
evsi_ar_fun <- function (evsi, add_fu) {
 
  library(drc)
  
  # fit asymptotic regression to EVSI point estimates
  df <- as.data.frame(cbind("y" = evsi[2,], "x" = evsi[1,]))
  ar_mod_evsi <- drm(y ~ x, data = df, fct = AR.3(fixed = c(NA, NA, NA)), pmodels = list(~1, ~1,~1))
  new_times <- data.frame("x" = seq.int(min(add_fu), max(add_fu), 1)) #seq.int(min(add_fu)
  evsi_ar <- data.frame("times" = new_times, "evsi" = predict(ar_mod_evsi, newdata = new_times), group = 1)
  
  # fit asymptotic regression to lower limit
  df <- as.data.frame(cbind("y" = evsi[4,], "x" = evsi[1,]))
  ar_mod_lower <- drm(y ~ x, data = df, fct = AR.3(fixed = c(NA, NA, NA)), pmodels = list(~1, ~1,~1)) #ar_mod_evsi$fit$par[3]
  evsi_ar$lower <- predict(ar_mod_lower, newdata = new_times)

  # fit asymptotic regression to upper limit
  df <- as.data.frame(cbind("y" = evsi[5,], "x" = evsi[1,]))
  ar_mod_upper <- drm(y ~ x, data = df, fct = AR.3(fixed = c(NA, NA, NA)), pmodels = list(~1, ~1,~1)) #ar_mod_evsi$fit$par[3]
  evsi_ar$upper <- predict(ar_mod_upper, newdata = new_times)
  
  # fit spline model to SE for EVSI values
  df <- as.data.frame(cbind("y" = evsi[3,], "x" = evsi[1,]))
  evsi_ar$se <- spline(x = df$x, y =df$y, xout =  new_times$x)$y
  
  # OS events
  evsi_ar$os_events <- round(spline(evsi[1,], evsi[6,], xout = seq.int(min(add_fu), max(add_fu), 1))$y)
  evsi_ar$os_events_se <- round(spline(evsi[1,], evsi[7,], xout = seq.int(min(add_fu), max(add_fu), 1))$y, 3)
  
  # PFS events
  if(nrow(evsi) > 7) {
    evsi_ar$pfs_events <- round(spline(evsi[1,], evsi[8,], xout = seq.int(min(add_fu), max(add_fu), 1))$y)
    evsi_ar$pfs_events_se <- round(spline(evsi[1,], evsi[9,], xout = seq.int(min(add_fu), max(add_fu), 1))$y, 3)
  }
  
  colnames(evsi_ar)[1] <- "time"
  
  return(evsi_ar)
}


####################################################################################
# function for plotting the EVSI
####################################################################################
evsi_plot_fun <- function (evsi_ar, pevpi = NULL) {
  
  evsi_ar$group <- "EVSI"
  
  if(!is.null(pevpi)) {
    
    # create a dataframe with partial EVPI values
    pevpi_temp <- evsi_ar
    pevpi_temp$evsi <- pevpi$pevpi
    pevpi_temp$upper <- pevpi$upper
    pevpi_temp$lower <- pevpi$lower
    pevpi_temp$group <- "Partial EVPI"
    
    # merge the EVSI and partial EVPI dataframes
    df_evsi_temp <- rbind(evsi_ar, pevpi_temp)
    
  } else {
    df_evsi_temp <- evsi_ar
  }
    
  # select colors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  col <- gg_color_hue(2) 
  
  # plot the EVSI and partial EVPI + confidence bounds
  ggplot(data=df_evsi_temp, aes(x=time, y=evsi, group = group)) +
    theme_light() + 
    geom_line(aes(color = group), linewidth = 0.5) + # aes(color=group),
    xlab("Additional follow-up (months)") +
    ylab("Per Patient Net Benefit") +
    theme(legend.position = "top", legend.title = element_blank()) +
    ggtitle("") + scale_size(range=c(0.1, 2), guide="none") + 
    scale_x_continuous(breaks = seq(0,max(df_evsi_temp$time),12) ) + # , limits = c(min(add_fu),max(add_fu)) +
    #ylim(0,NA) +
    theme(axis.text.x = element_text(color="black")) +
    theme(axis.text.y = element_text(color="black")) +
    theme(text=element_text(size=18)) +
    geom_ribbon(aes(ymin=lower,ymax=upper, fill = group),alpha=0.2) #+
    #geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.8) 
}


#####################################################################################
# Population Expected Net Benefit of Sampling
#####################################################################################
enbs_fun <- function (evsi_ar, m_nb, c_fix, c_var,  c_var_time = NULL, c_var_event = NULL, c_rev, t_lag, inc_pop, prev_pop, dec_th, dr_voi,  reversal=1) {
  
  # stop conditions
  if(dec_th<max(evsi_ar[,1]) | is.null(dec_th)) {stop("The decision relevance horizon must be equal or greater than the maximum additional follow-up time.")}
  
  if(dr_voi>=1 | dr_voi<0 | is.null(dr_voi)) {stop("The annual discount rate is not below 1.")}
  
  if(min(c_var)<0 | is.null(c_var)) {stop("Monthly trial costs are not a positive.")}
  if(length(c_var) > 2) {stop("Monthly trial costs are not a range.")}
  
  if(min(c_fix)<0 | is.null(c_fix)) {stop("Fixed trial costs are not positive.")}
  if(length(c_fix) > 2) {stop("Fixed trial costs are not a range.")}
  
  if(min(c_rev)<0 | is.null(c_rev)) {stop("Decision reversal costs are not positive.")}
  if(length(c_rev) > 2) {stop("Decision reversal costs are not a range.")}
  
  if(min(inc_pop)<0 | is.null(inc_pop)) {stop("Monthly incident population is not positive.")}
  if(length(inc_pop) > 2) {stop("Monthly incident population is not a range.")}
  
  if(min(prev_pop)<0 | is.null(inc_pop)) {stop("Total prevalent population is not positive.")}
  if(length(prev_pop) > 2) {stop("Total prevalent population is not a range.")}
  
  if(!exists("evsi_ar") | is.null("evsi_ar")) {stop("The interpolated EVSI 'evsi_ar' must first be computed using the 'evsi_ar_fun' function.")}
  
  if(!is.null(c_var_time) & !is.null(c_var_event)) {stop("Only c_var_time or c_var_event can be specified, not both.")}
  
  # if no time or events are defined after which variable cost are incurred, set c_var_time = 0
  if(is.null(c_var_time) & is.null(c_var_event)) {c_var_time <- 0}

  
  ### Time lag (time between end of follow-up and decision making) ###
  
  # vector of times for calculations
  x_times <- evsi_ar[,1] #seq.int(min(add_fu), max(add_fu), 1) # times until maximum follow-up
  x_times_desc <- sort(x_times, decreasing = T) +  (dec_th - max(evsi_ar[,1]) - min(evsi_ar[,1]))
  
  # EVSI adjusted for time lag
  t_lag_sigma <- (max(t_lag) - min(t_lag)) / (2 * qnorm(0.975)) # SE for lag in reporting time
  lag_temp <- round(rnorm(5000, mean(t_lag), t_lag_sigma)) # sample reporting delay times
  df_x_times_desc <-  sapply(lag_temp, function (x){ # adjust decision relevance horizons
    pmax(x_times_desc - x, 0)
    #x_times_desc - x
  })
  df_evsi_temp <- apply(df_x_times_desc, 2, function (x) {evsi_ar$evsi * x})  # multiply adjusted decision relevance horizons with EVSI
  df_se_temp <- apply(df_x_times_desc, 2, function (x) {evsi_ar$se * x})  # multiply adjusted decision relevance horizons with SE
  
  # mean EVSI adjusted for delay in reporting
  evsi_delay <- rowMeans(df_evsi_temp) #rowMeans(df_x_times_desc) * evsi_ar$evsi 
 
  # SE due to uncertainty in GAM estimate
  se_gam_temp <- rowMeans(df_se_temp) 
  
  # SE due to uncertainty in lag time (time between end of follow-up and decision making)
  upper_temp <- apply(df_evsi_temp, 1, function (x) {quantile(x, 0.975) }) # upper range due to lag time
  lower_temp <- apply(df_evsi_temp, 1, function (x) {quantile(x, 0.025) }) # lower range due to lag time
  se_delay_temp <- (upper_temp - lower_temp) / (2 * qnorm(0.975)) # SE due to lag time

  # combine SE due to reporting delay with SE for GAM estimate of EVSI
  se_delay <- sqrt(se_delay_temp^2 + se_gam_temp^2)
  
  ### Population EVSI ###
  
  # mean population EVSI (accounting for both incidence and prevalence)
  evsi_pop <- evsi_delay * mean(inc_pop) + evsi_ar$evsi * mean(prev_pop) * rowMeans(df_x_times_desc)

  # SE for incident population EVSI
  upper_temp <- evsi_delay * max(inc_pop) 
  lower_temp <- evsi_delay * min(inc_pop)
  se_pop <- (upper_temp - lower_temp) / (2 * qnorm(0.975)) 
  se_pop <-  sqrt(se_pop^2 + (se_delay * mean(inc_pop))^2)
  
  # SE for incident + prevalent population EVSI
  upper_temp <- evsi_ar$upper * max(prev_pop) * ifelse(rowMeans(df_x_times_desc)>0, 1, 0)
  lower_temp <- evsi_ar$lower * min(prev_pop) * ifelse(rowMeans(df_x_times_desc)>0, 1, 0)
  se_prev <- (upper_temp - lower_temp) / (2 * qnorm(0.975)) 
  se_pop <-  sqrt(se_pop^2 + se_prev^2)
  
  # dataframe with population EVSI and SE
  pop_evsi <- data.frame(times = x_times,
                         evsi = evsi_pop, 
                         se = se_pop) # overall SE accounting for uncertainty in incidence, delayed reporting and GAM estimator

  ### ENBS AWR ###
  
  # variable trial costs
  if(!is.null(c_var_time)) {
    c_trial_upper <- max(c_var) * pmax(0, x_times - c_var_time)
    c_trial_lower <- min(c_var) * pmax(0, x_times - c_var_time)
    c_trial_var <- mean(c_var) * pmax(0, x_times - c_var_time)
    c_trial_var_sigma <- (c_trial_upper - c_trial_lower) / (2 * qnorm(0.975))  # SE for variable costs
  }
  
  if(!is.null(c_var_event)) {
    c_trial_upper <- cumsum(max(c_var) * pnorm(c_var_event, evsi_ar$os_events, evsi_ar$os_events_se, lower.tail = F))
    c_trial_lower <- cumsum(min(c_var) * pnorm(c_var_event, evsi_ar$os_events, evsi_ar$os_events_se, lower.tail = F))
    c_trial_var <- cumsum(mean(c_var) * pnorm(c_var_event, evsi_ar$os_events, evsi_ar$os_events_se, lower.tail = F))
    c_trial_var_sigma <- (c_trial_upper - c_trial_lower) / (2 * qnorm(0.975))  # SE for variable costs
  }
  
  # ENBS for AWR 
  c_enbs_awr <- reversal * pop_evsi$evsi - c_trial_var
  c_enbs_awr_sigma <- sqrt((reversal * pop_evsi$se)^2 + c_trial_var_sigma^2)

  # subtract reversal costs
  c_rev_sigma <- (max(c_rev) - min(c_rev)) / (2 * qnorm(0.975)) # SE for reversal costs
  
  c_enbs_awr <- c_enbs_awr - mean(c_rev)
  c_enbs_awr_sigma <- sqrt(c_enbs_awr_sigma^2 + c_rev_sigma^2)
  
  # apply discounting
  v_dr_voi <- 1/(1+dr_voi)^(x_times/12) # monthly discount rate 
  
  c_enbs_awr <- c_enbs_awr * v_dr_voi
  c_enbs_awr_sigma <- c_enbs_awr_sigma * v_dr_voi
  
  # subtract fixed trial setup costs
  c_fix_sigma <- (max(c_fix) - min(c_fix)) / (2 * qnorm(0.975)) # SE for fixed costs
  
  c_enbs_awr <- c_enbs_awr - mean(c_fix)
  c_enbs_awr_sigma <- sqrt(c_enbs_awr_sigma^2 + c_fix_sigma^2)
  
  # create dataframe for AWR
  df_awr <- data.frame("time" = x_times, "value" = c_enbs_awr, 
                       "lower" = (c_enbs_awr - qnorm(0.975) * c_enbs_awr_sigma), 
                       "upper" = (c_enbs_awr + qnorm(0.975) * c_enbs_awr_sigma), 
                       "group" = "AWR")
  
  ### ENBS OIR ###
  
  # cost due to withholding access
  v_inb <- max(colMeans(m_nb)) - colMeans(m_nb)[1] # incremental net benefits
  c_wait <-  v_inb * mean(inc_pop) * x_times
 
  # ENBS for OIR 
  c_enbs_oir <- pop_evsi$evsi - c_trial_var - c_wait
  c_enbs_oir_sigma <- sqrt(pop_evsi$se^2 + c_trial_var_sigma^2)
  
  # apply discounting
  c_enbs_oir <- c_enbs_oir * v_dr_voi
  c_enbs_oir_sigma <- c_enbs_oir_sigma * v_dr_voi
  
  # subtract fixed costs
  c_enbs_oir <- c_enbs_oir - mean(c_fix)
  c_enbs_oir_sigma <- sqrt(c_enbs_oir_sigma^2 + c_fix_sigma^2)
  
  # create dataframe for OIR
  df_oir <- data.frame("time" = x_times, "value" = c_enbs_oir, 
                       "lower" = (c_enbs_oir - qnorm(0.975) * c_enbs_oir_sigma), 
                       "upper" = (c_enbs_oir + qnorm(0.975) * c_enbs_oir_sigma), 
                       "group" = "OIR")
  
  ### Merge datasets ###  
  df_enbs <- rbind(df_awr, df_oir)
  df_enbs$group <- as.factor(df_enbs$group)
  
  ### compute maximum ENBS for AWR and OIR 
  
  # AWR
  temp <- cbind(subset(df_enbs, group == "AWR", select = c(time, value)))
  x_1 <- temp$time[which.max(temp$value)]
  y_1 <- max(temp$value)
  fu1 <- x_1 - mean(t_lag)
  
  if(y_1<=0) {
    x_1 <- 0
    fu1 <- 0
  }
  
  # OIR
  temp <- cbind(subset(df_enbs, group == "OIR", select = c(time, value)))
  x_2 <- temp$time[which.max(temp$value)]
  y_2 <- max(temp$value)
  fu2 <- x_2 - mean(t_lag)
  if(y_2<=1) {
    x_2 <- 0
    fu2 <- 0
  }
  
  ### plot ENBS ###
  anno1   <- deparse(bquote(paste(~ italic("t"^"*") == .(round(x_1,0)) )) )
  anno2   <- deparse(bquote(paste(~ italic("t"^"*") == .(round(x_2,0)) )) )
  
  y_loc <- pmax((max(df_enbs$value) - min(df_enbs$value)),  max(df_enbs$value)) * 0.06 #(max(df_enbs$value) - min(df_enbs$value)) * 0.08
  
  pmax((max(df_enbs$value) - min(df_enbs$value)),  max(df_enbs$value))
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  col <- gg_color_hue(2) 
  
  # plot ENBS predictions
  p <- ggplot(data=df_enbs, aes(x=time, y=value, group = group)) +
    theme_light() + 
    geom_line(aes(color = group), linewidth = 0.5) + # aes(color=group),
    xlab("Additional follow-up (months)") +
    ylab("Population ENBS") +
    theme(legend.position = "top", legend.title = element_blank()) +
    ggtitle("") + scale_size(range=c(0.1, 2), guide="none") + 
    scale_x_continuous(breaks = seq(0,max(df_enbs$time),12) ) + # , limits = c(min(add_fu),max(add_fu)) +
    #ylim(0,NA) +
    theme(axis.text.x = element_text(color="black")) +
    theme(axis.text.y = element_text(color="black")) +
    theme(text=element_text(size=18)) +
    geom_ribbon(aes(ymin=lower,ymax=upper, fill = group),alpha=0.2) +
    geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.8) 
  
  if(x_1 < max(evsi_ar[,1]) & x_1 > 0 & y_1 > 0) {p <- p +
    annotate("text", label = anno1, x = x_1, y = y_1 + y_loc, size = 6, parse = T, color = col[1]) +
    geom_point(aes(x = x_1, y = y_1), colour = "black", size = 1)}
  
  if(x_2 < max(evsi_ar[,1]) & x_2 > 0 & y_2 > 0) {p <- p +
    annotate("text", label = anno2, x = x_2, y = y_2 + y_loc, size = 6, parse = T, color = col[2]) +
    geom_point(aes(x = x_2, y = y_2), colour = "black", size = 1)}
  
  p
  
}
