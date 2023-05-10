# Copyright (c) 2023 Mathyn Vervaart
# Licensed under the MIT License

#--------------------------------------------------------------------------------------------#
# OS functions
#--------------------------------------------------------------------------------------------#

##############################################################################################
# function for generating overall survival datasets using spline interpolation
##############################################################################################
interpol_os_data_fun <- function (curves, start_times, l_dropout, l_enroll, t, seed1, cycle2month, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))   # number of patients at risk at t1
  cycles <-  0:(length(curves[[1]][,1])-1)                # number of model cycles
  n_sim <- ncol(curves[[1]])                              # number of PA simulations
  arm_indic <- rep(1:length(curves), 2, each = 2)         # trial arm indicator
  arms <- unique(arm_indic)
    
  set.seed(seed1) # set the seed1

  ### sample survival times
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(curves), function (y)  {
    apply(curves[[y]], 2, function (x) runif(atrisk[[y]], min(x), spline(cycles, x, xout=(start_times[[y]]), ties = max, method="hyman")$y)) 
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(curves), function (y) {
    sapply(1:n_sim, function (x) {
      spline(round(curves[[y]][,x], 100), cycles, xout = rand_num[[y]][, x], ties = max, method="hyman")$y
    })
  })

 
  ##### incorporate administrative censoring, random censoring and trial enrollment, if applicable
  
  # sample random censoring times
  if (!is.null(l_enroll)) {
    
    # sample random censoring rates
    l_enroll <- lapply(l_enroll, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    # set censoring time as min(random censoring time, administrative censoring time)
    enroll_times <- lapply(1:length(l_enroll), function (y) {
      sapply(l_enroll[[y]], function (x)
        rexp(atrisk[[y]], x))
    })
  }
  
  # sample enrollment times
  if (!is.null(l_dropout)) {
    
    # sample random censoring rates
    l_dropout <- lapply(l_dropout, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    rand_cens_times <- lapply(1:length(l_dropout), function (y) {
      sapply(l_dropout[[y]], function (x)
        rexp(atrisk[[y]], x))
    })
    
  }
  
  if (!is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(random censoring time, administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        start_times[[y]] + pmin(rand_cens_times[[y]][,x], pmax(t -  enroll_times[[y]][,x], 0))
      })
    })
  }
  
  if (!is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        start_times[[y]] + pmax(t -  enroll_times[[y]][,x], 0)
      })
    })
  }
  
  if (is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        start_times[[y]] + pmin(rand_cens_times[[y]][,x], t)
      })
    })
  }
  
  if (is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time at administrative censoring time 
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) 
        start_times[[y]] + t
      )
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      ifelse(surv_times[[y]][,x] > cens_times[[y]][,x], 0, 1)
    )
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      pmin(surv_times[[y]][,x], cens_times[[y]][,x])
    )
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times[[y]])))
    })
  })
  summ_stat <- as.data.frame(t(do.call(rbind, summ_stat)))
  
} 

##############################################################################################
# function for generating overall survival datasets using discrete sampling
##############################################################################################
discrete_os_data_fun <- function (curves, start_times, l_dropout, l_enroll, t, seed1, fast = TRUE, cycle2month, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk <- sapply(start_times, function (x) length(x))   # number of patients at risk at t1
  cycles <-  0:(length(curves[[1]][,1])-1)                # number of model cycles
  n_sim <- ncol(curves[[1]])                              # number of PA simulations
  arm_indic <- rep(1:length(curves), 2, each = 2)         # trial arm indicator
  arms <- unique(arm_indic)
   
  # function to find the half-cycle time corresponding to the uniform value rounded down to the nearest survival probability
  minpositive_fun = function(x) which.min(x[x > 0]) - 0.5 # minimum of positive values
  #minabs_fun = function(x) which.min(abs(x)) - 0.5 # minimum of absolute values
  
  if(fast == TRUE) {
   
    set.seed(seed1) # set the seed1
    
    # sample model cycles (i.e. survival times) using the densities
    surv_times <- lapply(1:length(curves), function (y)  {

      start_t_temp <- (round(mean(start_times[[y]])) + 0.5) #

      apply(curves[[y]], 2, function (x)
        sample(
          start_t_temp:(max(cycles) - 0.5),
          size = atrisk[[y]],
          replace = T,
          prob = abs(diff(x))[ceiling(start_t_temp):length(abs(diff(x)))] #
        ))
    })
    
  } else {

    set.seed(seed1) # set the seed1
    
    # Sample values from a uniform distribution for each survival curve
    rand_num <- lapply(1:length(curves), function (y)  {
      apply(curves[[y]], 2, function (x) runif(atrisk[[y]], 
                                             min(x), 
                                             x[round(start_times[[y]])+1])) 
    })
    
    # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
    surv_times <- lapply(1:length(curves), function (y) {
      sapply(1:n_sim, function (x) {
        curves_temp <- curves[[y]][,x]
        sapply(1:atrisk[[y]], function (n) {
          minpositive_fun(curves_temp - rand_num[[y]][, x][n])
        })
      })
    })
    
  }
  
  ##### incorporate administrative censoring, random censoring and trial enrollment, if applicable
  
  # sample random censoring times
  if (!is.null(l_enroll)) {
    
    # sample random censoring rates
    l_enroll <- lapply(l_enroll, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    # set censoring time as min(random censoring time, administrative censoring time)
    enroll_times <- lapply(1:length(l_enroll), function (y) {
      sapply(l_enroll[[y]], function (x)
        rexp(atrisk[[y]], x))
    })
  }
  
  # sample enrollment times
  if (!is.null(l_dropout)) {
    
    # sample random censoring rates
    l_dropout <- lapply(l_dropout, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    rand_cens_times <- lapply(1:length(l_dropout), function (y) {
      sapply(l_dropout[[y]], function (x)
        rexp(atrisk[[y]], x))
    })
    
  }
  
  if (!is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(random censoring time, administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        start_times[[y]] + pmin(rand_cens_times[[y]][,x], pmax(t -  enroll_times[[y]][,x], 0))
      })
    })
  }
  
  if (!is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        start_times[[y]] + pmax(t -  enroll_times[[y]][,x], 0)
      })
    })
  }
  
  if (is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        start_times[[y]] + pmin(rand_cens_times[[y]][,x], t)
      })
    })
  }
  
  if (is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time at administrative censoring time 
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) 
        start_times[[y]] + t
      )
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      ifelse(surv_times[[y]][,x] > cens_times[[y]][,x], 0, 1)
    )
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      pmin(surv_times[[y]][,x], cens_times[[y]][,x])
    )
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(start_times[[y]])))
    })
  })
  summ_stat <- as.data.frame(t(do.call(rbind, summ_stat)))
  
} 


#--------------------------------------------------------------------------------------------#
# OS and PFS functions
#--------------------------------------------------------------------------------------------#


############################################################################################################
# function for generating overall survival and progression-free survival datasets using spline interpolation
############################################################################################################
interpol_os_pfs_data_fun <- function (l_os, l_pfs, l_start_os, l_start_pfs, l_dropout, l_enroll, t, seed1, cycle2month, ...) {
  
  gc()
  
  print(paste("Generating datasets.."))
  
  atrisk_os <- sapply(l_start_os, function (x) length(x))    # number of patients at risk at t1 for OS
  atrisk_pfs <- sapply(l_start_pfs, function (x) length(x))  # number of patients at risk at t1 for OS
  cycles <-  0:(length(l_os[[1]][,1])-1)                     # number of model cycles
  n_sim <- ncol(l_os[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(l_os), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  
  set.seed(seed1) # set the seed
  
  ### Overall survival
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(l_os), function (y)  {
    apply(l_os[[y]], 2, function (x) runif(atrisk_os[[y]], min(x), spline(cycles, x, xout=(l_start_os[[y]]), ties = max, method="hyman")$y)) 
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(l_os), function (y) {
    sapply(1:n_sim, function (x) {
      spline(round(l_os[[y]][,x], 100), cycles, xout = rand_num[[y]][, x], ties = max, method="hyman")$y
    })
  })
  
  ##### incorporate administrative censoring, random censoring and trial enrollment, if applicable
 
  # sample enrollment times
  if (!is.null(l_enroll)) {
    
    # sample enrollment times
    l_enroll <- lapply(l_enroll, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    # set censoring time as min(random censoring time, administrative censoring time)
    enroll_times <- lapply(1:length(l_enroll), function (y) {
      sapply(l_enroll[[y]], function (x)
        rexp(atrisk_os[[y]], x))
    })
  }
  
  # sample random censoring times
  if (!is.null(l_dropout)) {
    
    # sample random censoring rates
    l_dropouts <- lapply(l_dropout, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    rand_cens_times <- lapply(1:length(l_dropouts), function (y) {
      sapply(l_dropouts[[y]], function (x)
        rexp(atrisk_os[[y]], x))
    })
    
  }
  
  if (!is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(random censoring time, administrative censoring time minus enrollment time)
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_os[[y]] + pmin(rand_cens_times[[y]][,x], pmax(t -  enroll_times[[y]][,x], 0))
      })
    })
  }
  
  if (!is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_os[[y]] + pmax(t -  enroll_times[[y]][,x], 0)
      })
    })
  }
  
  if (is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_os[[y]] + pmin(rand_cens_times[[y]][,x], t)
      })
    })
  }
  
  if (is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time at administrative censoring time 
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) 
        l_start_os[[y]] + t
      )
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      ifelse(surv_times[[y]][,x] > cens_times[[y]][,x], 0, 1)
    )
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      pmin(surv_times[[y]][,x], cens_times[[y]][,x])
    )
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_os <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(l_start_os[[y]])))
    })
  })
  
  ### Progression-free survival
  
  set.seed(seed1+1) # set the seed
  
  # Sample values from a uniform distribution for each survival curve
  rand_num <- lapply(1:length(l_pfs), function (y)  {
    apply(l_pfs[[y]], 2, function (x) runif(atrisk_pfs[[y]], min(x), spline(cycles, x, xout=(l_start_pfs[[y]]), ties = max, method="hyman")$y)) 
  })
  
  # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
  surv_times <- lapply(1:length(l_pfs), function (y) {
    sapply(1:n_sim, function (x) {
      spline(round(l_pfs[[y]][,x], 100), cycles, xout = rand_num[[y]][, x], ties = max, method="hyman")$y
    })
  })

  ##### incorporate administrative censoring, random censoring and trial enrollment, if applicable
  
  # sample random censoring times
  if (!is.null(l_dropout)) {
    
    # generate random censoring times
    rand_cens_times <- lapply(1:length(l_dropouts), function (y) {
      sapply(l_dropouts[[y]], function (x)
        rexp(atrisk_pfs[[y]], x))
    })
    
  }
  
  if (!is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(random censoring time, administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
      l_start_pfs[[y]] + pmin(rand_cens_times[[y]][,x], pmax(t -  enroll_times[[y]][,x], 0))
      })
    })
  }
  
  if (!is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_pfs[[y]] + pmax(t -  enroll_times[[y]][,x], 0)
      })
    })
  }
  
  if (is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_pfs[[y]] + pmin(rand_cens_times[[y]][,x], t)
      })
    })
  }
  
  if (is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time at administrative censoring time 
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) 
        l_start_pfs[[y]] + t
      )
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      ifelse(surv_times[[y]][,x] > cens_times[[y]][,x], 0, 1)
    )
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      pmin(surv_times[[y]][,x], cens_times[[y]][,x])
    )
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_prog <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(l_start_pfs[[y]])))
    })
  })
  
  # combine summary statistics for OS and for progression in one dataframe
  summ_stat <- as.data.frame(t(do.call(rbind,  c(summ_stat_os, summ_stat_prog))))

} # end data generation function


#########################################################################################################
# function for generating overall survival and progression-free survival datasets using discrete sampling
#########################################################################################################
discrete_os_pfs_data_fun <- function (l_os, l_pfs, l_start_os, l_start_pfs, l_dropout, l_enroll, t, seed1, fast = TRUE, cycle2month, ...) {
  
  gc()
 
  print(paste("Generating datasets.."))
  
  atrisk_os <- sapply(l_start_os, function (x) length(x))     # number of patients at risk at t1 for OS
  atrisk_pfs <- sapply(l_start_pfs, function (x) length(x))   # number of patients at risk at t1 for OS
  cycles <-  0:(length(l_os[[1]][,1])-1)                     # number of model cycles
  n_sim <- ncol(l_os[[1]])                                   # number of PA simulations
  arm_indic <- rep(1:length(l_os), 2, each = 2)              # trial arm indicator
  arms <- unique(arm_indic)
  
  # function to find the half-cycle time corresponding to the uniform value rounded down to the nearest survival probability
  minpositive_fun <- function(x) which.min(x[x > 0]) - 0.5 # minimum of positive values
  #minabs_fun <- function(x) which.min(abs(x)) - 0.5 # minimum of absolute values
  
  ### Overall survival
  
  if(fast == TRUE) {
    
    set.seed(seed1) # set the seed1
    
    # sample model cycles (i.e. survival times) using the densities
    surv_times <- lapply(1:length(l_os), function (y)  {
      
      start_t_temp <- (round(mean(l_start_os[[y]])) + 0.5) #
      
      apply(l_os[[y]], 2, function (x)
        sample(
          start_t_temp:(max(cycles) - 0.5), #0.5:(max(cycles) - 0.5),
          size = atrisk_os[[y]],
          replace = T,
          prob = abs(diff(x))[ceiling(start_t_temp):length(abs(diff(x)))] # prob = abs(diff(x))
        ))
    })
    
  } else {
    
    set.seed(seed1) # set the seed1
    
    # Sample values from a uniform distribution for each survival curve
    rand_num <- lapply(1:length(l_os), function (y)  {
      apply(l_os[[y]], 2, function (x) runif(atrisk_os[[y]], 
                                                  min(x), 
                                                  x[round(l_start_os[[y]])+1])) 
    })
    
    
    # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
    surv_times <- lapply(1:length(l_os), function (y) {
      sapply(1:n_sim, function (x) {
        curves_temp <- l_os[[y]][,x]
        sapply(1:atrisk_os[[y]], function (n) {
          minpositive_fun(curves_temp - rand_num[[y]][, x][n])
        })
      })
    })
    
  }
  
  ##### incorporate administrative censoring, random censoring and trial enrollment, if applicable
  
  # sample enrollment times
  if (!is.null(l_enroll)) {
    
    # sample enrollment rates
    l_enroll <- lapply(l_enroll, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    # set censoring time as min(random censoring time, administrative censoring time)
    enroll_times <- lapply(1:length(l_enroll), function (y) {
      sapply(l_enroll[[y]], function (x)
        rexp(atrisk_os[[y]], x))
    })
  }
  
  # sample random censoring times
  if (!is.null(l_dropout)) {
    
    # sample random censoring rates
    l_dropouts <- lapply(l_dropout, function (x) rgamma(n_sim, x[1], x[2]/cycle2month))
    
    # generate random censoring times
    rand_cens_times <- lapply(1:length(l_dropouts), function (y) {
      sapply(l_dropouts[[y]], function (x)
        rexp(atrisk_os[[y]], x))
    })
    
  }
  
  if (!is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(random censoring time, administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_os[[y]] + pmin(rand_cens_times[[y]][,x], pmax(t -  enroll_times[[y]][,x], 0))
      })
    })
  }
  
  if (!is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_os[[y]] + pmax(t -  enroll_times[[y]][,x], 0)
      })
    })
  }
  
  if (is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_os[[y]] + pmin(rand_cens_times[[y]][,x], t)
      })
    })
  }
  
  if (is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time at administrative censoring time 
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) 
        l_start_os[[y]] + t
      )
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      ifelse(surv_times[[y]][,x] > cens_times[[y]][,x], 0, 1)
    )
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      pmin(surv_times[[y]][,x], cens_times[[y]][,x])
    )
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_os <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(l_start_os[[y]])))
    })
  })
  
  ### Progression-free survival
  if(fast == TRUE) {
    
    set.seed(seed1+1) # set the seed1
    
    # sample model cycles (i.e. survival times) using the densities
    surv_times <- lapply(1:length(l_pfs), function (y)  {
      
      start_t_temp <- (round(mean(l_start_pfs[[y]])) + 0.5) #
      
      apply(l_pfs[[y]], 2, function (x)
        sample(
          start_t_temp:(max(cycles) - 0.5), #0.5:(max(cycles) - 0.5),
          size = atrisk_pfs[[y]],
          replace = T,
          prob = abs(diff(x))[ceiling(start_t_temp):length(abs(diff(x)))] # prob = abs(diff(x))
        ))
    })
    
  } else {
    
    set.seed(seed1+1) # set the seed1
    
    # Sample values from a uniform distribution for each survival curve
    rand_num <- lapply(1:length(l_pfs), function (y)  {
      apply(l_pfs[[y]], 2, function (x) runif(atrisk_pfs[[y]], 
                                                   min(x), 
                                                   x[round(l_start_pfs[[y]])+1])) 
    })
    
    # Look up the sampled values in the vector of survival probabilities and record the corresponding survival times (i.e. model cycles)
    surv_times <- lapply(1:length(l_pfs), function (y) {
      sapply(1:n_sim, function (x) {
        curves_temp <- l_pfs[[y]][,x]
        sapply(1:atrisk_pfs[[y]], function (n) {
          minpositive_fun(curves_temp - rand_num[[y]][, x][n])
        })
      })
    })
    
  }
  
  ##### incorporate administrative censoring, random censoring and trial enrollment, if applicable
  
  # sample random censoring times
  if (!is.null(l_dropout)) {
    
    # generate random censoring times
    rand_cens_times <- lapply(1:length(l_dropouts), function (y) {
      sapply(l_dropouts[[y]], function (x)
        rexp(atrisk_pfs[[y]], x))
    })
    
  }
  
  if (!is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(random censoring time, administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_pfs[[y]] + pmin(rand_cens_times[[y]][,x], pmax(t -  enroll_times[[y]][,x], 0))
      })
    })
  }
  
  if (!is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_pfs[[y]] + pmax(t -  enroll_times[[y]][,x], 0)
      })
    })
  }
  
  if (is.null(l_enroll) & !is.null(l_dropout)) {
    
    # set censoring time as min(administrative censoring time minus enrollment time)
    cens_times <-lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) {
        l_start_pfs[[y]] + pmin(rand_cens_times[[y]][,x], t)
      })
    })
  }
  
  if (is.null(l_enroll) & is.null(l_dropout)) {
    
    # set censoring time at administrative censoring time 
    cens_times <- lapply(1:length(surv_times), function (y) {
      sapply(1:n_sim, function (x) 
        l_start_pfs[[y]] + t
      )
    })
    
  }
  
  # compute censoring indicator for each dataset
  cens <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      ifelse(surv_times[[y]][,x] > cens_times[[y]][,x], 0, 1)
    )
  })
  
  # censor the survival times
  surv_times <- lapply(1:length(surv_times), function (y) {
    sapply(1:n_sim, function (x)
      pmin(surv_times[[y]][,x], cens_times[[y]][,x])
    )
  })
  
  # create datasets
  datasets <- lapply(1:length(surv_times), function (y) {
    lapply(1:n_sim, function (x) data.frame(cbind(surv_times[[y]][,x], cens[[y]][,x])))
  })
  
  # Compute a summary statistic using the number of observed events and time at risk for each dataset
  summ_stat_prog <- lapply(1:length(datasets), function (y) {
    sapply(datasets[[y]], function (x) {
      par <- c(sum(x[, 2]), (sum(x[, 1]) - sum(l_start_pfs[[y]])))
    })
  })
  
  # combine summary statistics for OS and for progression in one dataframe
  summ_stat <- as.data.frame(t(do.call(rbind,  c(summ_stat_os, summ_stat_prog))))
  
} # end data generation function


##############################################################################################
# Function for identify implausible combinations of OS-PFS datasets
##############################################################################################
check_os_pfs_fun <- function (summ_stat, l_start_os, l_start_pfs) {
  
  trunc_times <- c(l_start_os, l_start_pfs)
  
  even_num_os <- seq(2, (ncol(summ_stat) / 2), 2)
  uneven_num_os <- seq(1, (ncol(summ_stat) / 2), 2)
  even_num_pfs <- seq(((ncol(summ_stat) / 2)) + 2, ncol(summ_stat) , 2)
  uneven_num_pfs <- seq(((ncol(summ_stat) / 2) + 1), ncol(summ_stat) , 2)
  
  # identify implausible datasets where the sum of observed progression times exceeds total time at risk for OS
  err_trisk <- lapply(1:length(even_num_os), function (i) {
    which(summ_stat[, even_num_os[i]] < summ_stat[, even_num_pfs[i]])
  }) 
  
  # merge error indices
  err_index <- lapply(1:length(err_trisk), function (i) err_trisk[[i]])
  err_index <- rep(err_index, 2, each = 1)
  
}
