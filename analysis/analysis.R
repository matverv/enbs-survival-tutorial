# Copyright (c) 2023 Mathyn Vervaart
# Licensed under the MIT License

# ------------------------------------------------------------------
# Probabilistic Analysis
# ------------------------------------------------------------------
# Load the cost-effectiveness model functions
source("R/ce_model_functions_pembro.R")

# Run a probabilistic analysis of size "K"
set.seed(123)            # set the seed for reproducibility
K    <- 10e3             # number of simulations
l_pa <- run_pa_pembro(K) # run the probabilistic analysis
str(l_pa)                # inspect the structure of the PA output

# Calculate net benefits
thresh <- 30000                          # willingness-to-pay threshold
m_nb   <- l_pa$m_e - (l_pa$m_c / thresh) # net health benefits 

# ------------------------------------------------------------------
# EVPPI
# ------------------------------------------------------------------

# Compute mean survival 
df_mean_surv <- data.frame(os1 = colSums(l_pa$l_surv$m_os1), 
                           os2 = colSums(l_pa$l_surv$m_os2))

# Compute EVPPI
df_evppi <- voi::evppi(outputs = m_nb, inputs = df_mean_surv, 
                       pars = names(df_mean_surv), se = TRUE)
df_evppi

# ------------------------------------------------------------------
# Data simulation settings
# ------------------------------------------------------------------

# List of matrices with overall survival (OS) probabilities
l_surv <- list(os1 = (l_pa$l_surv$m_os1), os2 = (l_pa$l_surv$m_os2))

# Load the reconstructed OS data
load("data/ipd_os.RData")
str(ipd_os) # inspect the structure of the data

# Extract observed times for patients without event for each treatment
v_timeRisk_os1 <- with(ipd_os, times[treat == 1 & event == 0])
v_timeRisk_os2 <- with(ipd_os, times[treat == 2 & event == 0])

# Number of trial dropouts
dropouts1 <- 20
dropouts2 <- 15

# Function for removing dropouts from the observed follow-up times 
remove_dropouts <- function(times, dropouts) {
  # Calculate indices to remove dropouts by evenly spacing them out
  indices_to_remove <- seq(from = floor(length(times) / (dropouts + 1)), 
                           by = floor(length(times) / (dropouts + 1)), 
                           length.out = dropouts)
  return(times[-indices_to_remove])
}

# Apply the function to remove dropouts from observed times
l_atrisk_times <- list(remove_dropouts(v_timeRisk_os1, dropouts1),
                       remove_dropouts(v_timeRisk_os2, dropouts2))

# Gamma hyperparameters for exponential dropout rates (number of dropouts, time at risk in months)
l_dropout_pars <- list(c(dropouts1, with(ipd_os, sum(times[treat == 1]))),
                       c(dropouts2, with(ipd_os, sum(times[treat == 2]))))

# Specify 4 additional follow-up times from 3 to 60 months, evenly spaced on the log scale
v_add_fu <- round(exp(seq(log(3), log(60), length.out = 4)))

# ------------------------------------------------------------------
# Survival data simulation
# ------------------------------------------------------------------

# Survival data generation function for an ongoing trial
gen_surv_data_ongoing <- function(surv_probs, atrisk_times, dropout_pars, add_fu_months, cycles_per_year) {
  # Arguments:
    # surv_probs:      Matrix of survival probabilities for a given treatment and survival outcome. The number of rows 
    #                  corresponds to the model cycles, and the number of columns corresponds to the PA simulations.
    # atrisk_times:    Vector of observed follow-up times in months for the patients at risk corresponding to the treatment 
    #                  and survival outcome in `surv_probs`.
    # dropout_pars:    Gamma hyperparameters for the exponential dropout rate (number of dropouts, time at risk in months).
    # add_fu_months:   Vector of additional follow-up periods in months for applying administrative censoring.
    # cycles_per_year: Number of model cycles per year; used to convert time between months and model cycles.
  
  # Define variables
  atrisk <- length(atrisk_times)                      # number of patients at risk
  cycle2month <- 12 / cycles_per_year                 # cycle to month conversion
  v_cycles <- 0:(nrow(surv_probs) - 1) * cycle2month  # vector of model cycles in months
  n_sim <- ncol(surv_probs)                           # number of PA simulations

  # Generate random numbers uniformly distributed between survival probabilities 
  # corresponding to 'at risk' times and the model time horizon
  m_rand_num <- apply(surv_probs, 2, function(surv_probs) { 
    runif(atrisk, min(surv_probs),
          spline(v_cycles, surv_probs, xout = atrisk_times, method = "hyman")$y)
  })
  
  # Interpolate the vectors of survival probabilities at the random uniform numbers using 
  # monotone cubic splines and record the interpolated cycle times
  m_surv_times <- sapply(seq_len(n_sim), function(k) { 
    spline(round(surv_probs[, k], digits = 100),
           v_cycles, xout = m_rand_num[, k], method = "hyman")$y
  })

  # Sample dropout rates and then dropout times from exponential distributions
  v_dropout_rates <- rgamma(n_sim, shape = dropout_pars[1], rate = dropout_pars[2])
  m_dropout_times <- sapply(v_dropout_rates, function(dropout_rate) rexp(atrisk, dropout_rate))

  # Initialize summary statistics list
  l_summ_stat <- vector("list", length(add_fu_months))  
  for (i in seq_along(add_fu_months)) {
    # Censoring times: minimum of dropout time and additional follow-up
    m_cens_times <- atrisk_times + pmin(m_dropout_times, add_fu_months[i])
    # Censoring indicators: 0 for censored, 1 for event
    m_cens_ind <- ifelse(m_surv_times > m_cens_times, 0, 1)
    # Apply censoring to survival times
    m_surv_times_censored <- pmin(m_surv_times, m_cens_times)
    # Compute summary statistics: events and time at risk
    df_summ_stat <- data.frame(
      events = colSums(m_cens_ind),
      timeRisk = colSums(pmin(m_surv_times, m_cens_times) - atrisk_times)
    )
    # Store summary statistics dataframe in the list
    l_summ_stat[[i]] <- df_summ_stat 
  }
  return(l_summ_stat) # Return the results
}

# Generate survival data for each treatment and survival outcome using corresponding list elements
l_summ_stat <- Map(f = gen_surv_data_ongoing, 
                   l_surv, l_atrisk_times, l_dropout_pars, 
                   MoreArgs = list(add_fu_months = v_add_fu, cycles_per_year = 52)) 

# Combine data frames for each follow-up time across treatments into a single data frame
l_summ_stat <- lapply(seq_len(length(l_summ_stat[[1]])), function(i) do.call(cbind, lapply(l_summ_stat, "[[", i)))
head(l_summ_stat[[1]])

# ------------------------------------------------------------------
# EVSI computation
# ------------------------------------------------------------------

# Function to compute EVSI using the evppi function from the 'voi' package
compute_evsi <- function(outputs, summ_stat, ...) {
  # Arguments:
    # outputs:   Either a matrix/dataframe of net benefits (rows=samples, columns=options), 
    #            or a list with 'c' and 'e' as cost/effect matrices and 'k' as a willingness-to-pay vector. 
    # summ_stat: A dataframe with summary statistics (rows=samples, columns=summary statistics)
    # ...:       Additional arguments passed onto 'voi::evppi'
  
  # Compute the EVSI using flexible nonparametric regression and format output
  evsi_output <- voi::evppi(outputs = outputs, inputs = summ_stat, 
                            pars = colnames(summ_stat), se = TRUE, ...)
  names(evsi_output)[names(evsi_output) == "evppi"] <- "evsi" 
  return(evsi_output) # Return the results
}

# Apply the EVSI function to the list of summary statistic dataframes and format output
l_evsi <- lapply(X = l_summ_stat, FUN = compute_evsi, outputs = m_nb)
df_evsi <- do.call(rbind, l_evsi) 
df_evsi$add_fu <- v_add_fu

# Function to interpolate the EVSI and standard errors
interpolate_evsi <- function(evsi_out) {
  # Argument:
    # evsi_out: A dataframe containing 'add_fu', 'evsi', and 'se' generated by 'compute_evsi()'
  
  # Asymptotic regression for EVSI and smoothing splines for se's
  ar_mod_evsi  <- drc::drm(evsi ~ add_fu, data = evsi_out, fct = drc::AR.3())
  se_spline    <- smooth.spline(evsi_out$add_fu, log(evsi_out$se))
  
  # New times for predictions
  v_new_times  <- seq(min(evsi_out$add_fu), max(evsi_out$add_fu), by = 1)
  
  # Interpolate the EVSI and se's
  v_preds_evsi <- predict(ar_mod_evsi, newdata = data.frame(add_fu = v_new_times))
  v_pred_se    <- exp(predict(se_spline, data.frame(add_fu = v_new_times))$y)[,1] 
  
  # Create a dataframe and add confidence bounds
  evsi_out_ar <- data.frame(
    add_fu = v_new_times,
    evsi   = v_preds_evsi,
    se     = v_pred_se,
    lower  = v_preds_evsi - qnorm(0.975) * v_pred_se,
    upper  = v_preds_evsi + qnorm(0.975) * v_pred_se
  )
  return(evsi_out_ar) # Return the results
}

# Interpolate the EVSI and standard errors
df_evsi_ar <- interpolate_evsi(evsi_out = df_evsi)

# Plot the EVSI and confidence bounds
ggplot2::ggplot(df_evsi_ar, aes(add_fu)) +
  theme_minimal() + 
  geom_line(aes(y = evsi), color = "#00BFC4") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00BFC4", alpha = 0.1) +
  geom_point(data = df_evsi, aes(y = evsi), color = "#C00000") +
  labs(x = "Additional follow-up (months)", y = "Expected Value of Sample Information") +
  theme(text=element_text(size=14))


# ------------------------------------------------------------------
# ENBS calculation
# ------------------------------------------------------------------

# Function to compute ENBS across additional follow-up times
enbs_fun <- function (evsi_ar, nb, c_fix, c_var, t_lag_awr, t_lag_oir,
                      inc_pop, prev_pop, dec_th, prob_rev, dr) {
  # Arguments:
    # evsi_ar:   Dataframe with 'add_fu' and 'evsi', generated by `interpolate_evsi()`
    # nb:        Net benefits matrix, rows = samples, columns = treatments (1st column = standard care).
    # c_fix:     Fixed trial cost (range or value).
    # c_var:     Monthly trial conduct cost (range or value).
    # t_lag_awr: Time between end of follow-up and change of practice for AWR (range or value).
    # t_lag_oir: Time between end of follow-up and change of practice for OIR (range or value).
    # inc_pop:   Monthly incidence of patients affected by the decision (range or value).
    # prev_pop:  Prevalent patient population affected by the decision (range or value).
    # dec_th:    Decision relevance horizon in months (range or value).
    # prob_rev:  Probability that approval can be reversed if AWR is chosen.
    # dr:        Annual discount rate.

  # Helper function to calculate effective decision horizons
  decision_horizons <- function(add_fu, t_lag, dec_th) {
    dh <- sort(add_fu, decreasing = TRUE) + (mean(dec_th) - sum(range(add_fu))) - mean(t_lag)
    return(pmax(dh, 0))
  }
  
  # Determine decision horizons for AWR and OIR
  v_dh_awr <- decision_horizons(evsi_ar$add_fu, t_lag_awr, dec_th)
  v_dh_oir <- decision_horizons(evsi_ar$add_fu, t_lag_oir, dec_th)
  
  # Calculate population EVSI for AWR and OIR
  v_pop_evsi_awr <- (mean(prev_pop) * evsi_ar$evsi + mean(inc_pop) * v_dh_awr * evsi_ar$evsi) * prob_rev
  v_pop_evsi_oir <- (mean(prev_pop) * evsi_ar$evsi) + (mean(inc_pop) * v_dh_oir * evsi_ar$evsi)
  
  # Calculate variable trial costs and monthly discount rates
  v_c_var <- mean(c_var) * evsi_ar$add_fu
  v_dr    <- 1 / (1 + dr) ^ (evsi_ar$add_fu / 12)
  
  # Determine the value of access to the optimal new treatment for AWR
  m_inb_new    <- matrix(nb[, -1] - nb[, 1]) # incremental net benefits compared to standard care
  v_inb_access <- pmax(colMeans(m_inb_new)) * (mean(prev_pop) + mean(inc_pop) * evsi_ar$add_fu)
  
  # Calculate ENBS for AWR
  v_enbs_awr <- (v_pop_evsi_awr + v_inb_access - v_c_var) * v_dr - mean(c_fix)
  
  # Calculate ENBS for OIR
  v_enbs_oir <- (v_pop_evsi_oir - v_c_var) * v_dr - mean(c_fix)

  # Compile ENBS values into a dataframe
  df_enbs <- data.frame(
    time  = rep(evsi_ar$add_fu, times = 2),
    value = c(v_enbs_awr, v_enbs_oir),
    group = rep(c("Approval with Research", "Only in Research"), each = length(evsi_ar$add_fu))
  )
  
  return(df_enbs) # Return the results
}

# Calculate the ENBS
df_enbs <- enbs_fun(evsi_ar   = df_evsi_ar,
                    nb        = m_nb,
                    c_fix     = c(100000) / thresh,
                    c_var     = 0,
                    t_lag_awr = c(6, 12) - 6,
                    t_lag_oir = c(4, 8) - 6, 
                    inc_pop   = round(277 * c(0.8, 1.2)), # https://www.nice.org.uk/guidance/ta650
                    prev_pop  = 0,
                    dec_th    = c(60, 120),
                    prob_rev  = 1,
                    dr        = 0.035)

# Find the maximum ENBS for AWR and OIR and corresponding follow-up
max_enbs <- aggregate(value ~ group, data = df_enbs, FUN = max)
optimal_fu <- sapply(split(df_enbs, df_enbs$group), function(x) x$time[which.max(x$value)])
df_max_enbs <- data.frame(group = max_enbs$group, time = optimal_fu, value = max_enbs$value)

# Plot the ENBS values with points showing the times and labels for maximum values
ggplot2::ggplot(df_enbs, aes(x = time, y = value, color = group)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.8) +
  geom_line() + # Draw lines to connect points in the same group
  geom_point(data = df_max_enbs, aes(x = time, y = value), shape = 19, size = 3) +
  geom_text(data = df_max_enbs, aes(x = time, y = value, label = paste0("t = ", time), color = group), 
            size = 4, vjust = -1, show.legend = FALSE) +
  labs(x = "Additional follow-up (months)", y = "Expected Net Benefit of Sampling") +
  theme_minimal() +
  theme(text=element_text(size=14), legend.position = "top", legend.title = element_blank())