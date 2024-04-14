############################################################################################
# Title:  Calculating the Expected Net Benefit of Sampling for Survival Data: 
#         A Tutorial and Case Study (2024)
# 
# This is the main code developed for the analysis presented in the paper.
#
# Authors:
# 
# Please cite the article when using this code.
#
# License: This code is licensed under the MIT License.
#
# Developed using:
#   R: Version 4.2.2 (2022-10-31)
#   RStudio: Version 2023.12.1+402 (2023.12.1+402), RStudio, Inc.
############################################################################################


# ------------------------------------------------------------------
# Install and/or load packages 
# ------------------------------------------------------------------
list_of_packages <- c("voi", "ggplot2", "drc")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, library, character.only = TRUE)


# ------------------------------------------------------------------
# Probabilistic Analysis
# ------------------------------------------------------------------

# Load the cost-effectiveness model functions
source("R/ce_fun_pembro.R")

# Run a probabilistic analysis of size "K"
set.seed(123)            # set the seed for reproducibility
K    <- 5e3              # number of simulations
l_pa <- run_pa_pembro(K) # run the probabilistic analysis (PA)
str(l_pa)                # inspect the structure of the PA output

# Calculate net benefits
thresh <- 30000                          # willingness-to-pay threshold
m_nb   <- l_pa$m_e - (l_pa$m_c / thresh) # net health benefits 


# ------------------------------------------------------------------
# EVPPI computation
# ------------------------------------------------------------------

# List of matrices with survival probabilities of interest (OS and PFS)
l_surv <- list(os1  = l_pa$l_surv$m_os1, 
               os2  = l_pa$l_surv$m_os2,
               pfs1 = l_pa$l_surv$m_pfs1, 
               pfs2 = l_pa$l_surv$m_pfs2)

# Compute mean survival across PA simulations for each element in `l_surv`
df_mean_surv <- data.frame(sapply(l_surv, colSums))

# Compute EVPPI
df_evppi <- voi::evppi(outputs = m_nb, inputs = df_mean_surv, 
                       pars = names(df_mean_surv), se = TRUE)
df_evppi


# ------------------------------------------------------------------
# Data simulation settings
# ------------------------------------------------------------------

# Load the reconstructed OS and PFS data
load("data/ipd_os.RData")
load("data/ipd_pfs.RData")

# Extract observed times for patients without event for each treatment
v_cens_times_os1  <- with(ipd_os, times[treat == 1 & event == 0])
v_cens_times_os2  <- with(ipd_os, times[treat == 2 & event == 0])
v_cens_times_pfs1 <- with(ipd_pfs, times[treat == 1 & event == 0])
v_cens_times_pfs2 <- with(ipd_pfs, times[treat == 2 & event == 0])

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
l_atrisk_times <- list(remove_dropouts(v_cens_times_os1, dropouts1),
                       remove_dropouts(v_cens_times_os2, dropouts2),
                       remove_dropouts(v_cens_times_pfs1, dropouts1),
                       remove_dropouts(v_cens_times_pfs2, dropouts2))

# Gamma hyperparameters for exponential dropout rates (number of dropouts, time at risk in months)
l_dropout_pars <- rep(list(c(dropouts1, with(ipd_os, sum(times[treat == 1]))),
                           c(dropouts2, with(ipd_os, sum(times[treat == 2])))), 
                      times = 2)

# Specify 4 additional administrative censoring times from 3 to 60 months, evenly spaced on the log scale
v_adm_cens_months <- round(exp(seq(from = log(3), to = log(60), length.out = 4)))


# ------------------------------------------------------------------
# Survival data simulation
# ------------------------------------------------------------------

# Survival data generation function for an ongoing trial
sim_surv_data_ongoing <- function(surv_probs, atrisk_times, dropout_pars, 
                                  adm_cens_months, cycles_per_year) {
  # Arguments:
    # surv_probs:      Matrix of survival probabilities for a single treatment and survival 
    #                  outcome. Rows = model cycles, columns = PA simulations.
    # atrisk_times:    Vector of observed follow-up times in months for the patients at risk 
    #                  corresponding to the treatment and survival outcome in `surv_probs`.
    # dropout_pars:    Gamma parameters for dropout rates (dropouts, time at risk in months).
    # adm_cens_months: Vector of additional administrative censoring times in months.
    # cycles_per_year: Model cycles per year, for cycle-month conversion.
  
  # Define variables
  atrisk      <- length(atrisk_times)                    # number of patients at risk
  cycle2month <- 12 / cycles_per_year                    # cycle to month conversion
  v_cycles    <- 0:(nrow(surv_probs) - 1) * cycle2month  # vector of model cycles in months
  n_sim       <- ncol(surv_probs)                        # number of PA simulations
  
  # Generate random numbers uniformly distributed between survival probabilities 
  # corresponding to 'at risk' times and the model time horizon
  m_rand_num <- apply(surv_probs, 2, function(surv_probs) { 
    runif(n   = atrisk, 
          min = min(surv_probs),
          max = spline(x      = v_cycles, 
                       y      = surv_probs, 
                       xout   = atrisk_times, 
                       method = "hyman")$y)
  })
  
  # Interpolate the vectors of survival probabilities at the random uniform numbers using 
  # monotone cubic splines and record the interpolated cycle times
  m_surv_times <- sapply(seq_len(n_sim), function(k) { 
    spline(x      = round(surv_probs[, k], digits = 100),
           y      = v_cycles, 
           xout   = m_rand_num[, k], 
           method = "hyman")$y
  })
  
  # Sample dropout rates and then dropout times from exponential distributions
  v_dropout_rates <- rgamma(n = n_sim, shape = dropout_pars[1], rate = dropout_pars[2])
  m_dropout_times <- sapply(v_dropout_rates, function(dropout_rate) {
    rexp(n = atrisk, rate = dropout_rate)
  })
  
  # Initialize summary statistics list
  l_summ_stat <- vector("list", length(adm_cens_months))  
  for (i in seq_along(adm_cens_months)) {
    # Censoring times: minimum of dropout time and additional follow-up
    m_cens_times <- atrisk_times + pmin(m_dropout_times, adm_cens_months[i])
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

# Map the data simulation function across lists for each treatment and outcome
l_summ_stat <- Map(f = sim_surv_data_ongoing, 
                   l_surv, l_atrisk_times, l_dropout_pars, 
                   MoreArgs = list(adm_cens_months = v_adm_cens_months, cycles_per_year = 52)) 

# Combine data frames for each follow-up time across treatments into a single data frame
l_summ_stat <- lapply(seq_len(length(l_summ_stat[[1]])), function(i) do.call(cbind, lapply(l_summ_stat, "[[", i)))
head(l_summ_stat[[1]])

# Total deaths required per trial protocol minus observed deaths
add_events <- 404 - 156

# Predicted number of events at each additional follow-up time
v_events <- sapply(l_summ_stat, function(df) mean(df$os1.events + df$os2.events))

# Estimated follow-up time in months for additional events
round(spline(v_events, v_adm_cens_months, xout = add_events)$y)

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

# Iterate through the list of summary statistic dataframes and calculate EVSI
# Note: this may take a while using the default Gaussian Process regression method.
# To use a faster alternative method, specify method = "earth" as an additional argument below, 
# although be aware that it might be less accurate.
l_evsi <- lapply(X = l_summ_stat, FUN = compute_evsi, outputs = m_nb)

# Combine the list of EVSI calculations into one dataframe
df_evsi <- do.call(rbind, l_evsi) 

# Append the administrative censoring times to the resulting dataframe and print the results
df_evsi$fu <- v_adm_cens_months
df_evsi

# Function to interpolate the EVSI and standard errors
interpolate_evsi <- function(evsi_out) {
  # Argument:
    # evsi_out: A dataframe containing 'fu', 'evsi', and 'se' generated by 'compute_evsi()'
  
  # Asymptotic regression for EVSI and smoothing splines for log-transformed se's
  ar_mod_evsi  <- drc::drm(evsi ~ fu, data = evsi_out, fct = drc::AR.3())
  se_spline    <- smooth.spline(evsi_out$fu, log(evsi_out$se))
  
  # New times for predictions
  v_new_times  <- seq(min(evsi_out$fu), max(evsi_out$fu), by = 1)
  
  # Interpolate the EVSI and se's
  v_pred_evsi <- predict(ar_mod_evsi, newdata = data.frame(fu = v_new_times))
  v_pred_se    <- exp(predict(se_spline, data.frame(fu = v_new_times))$y)[,1] 
  
  # Create a dataframe and add confidence bounds
  evsi_out_ar <- data.frame(
    fu = v_new_times,
    evsi   = v_pred_evsi,
    se     = v_pred_se,
    lower  = v_pred_evsi - qnorm(0.975) * v_pred_se,
    upper  = v_pred_evsi + qnorm(0.975) * v_pred_se
  )
  return(evsi_out_ar) # Return the results
}

# Interpolate the EVSI and standard errors
df_evsi_ar <- interpolate_evsi(evsi_out = df_evsi)

# Plot the EVSI and confidence bounds
ggplot2::ggplot(df_evsi_ar, aes(fu)) +
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
compute_enbs <- function (evsi_ar, nb, c_fix, c_var, t_lag_awr, t_lag_oir,
                          inc_pop, prev_pop, dec_th, prob_rev, dr) {
  # Arguments:
    # evsi_ar:   Dataframe with 'fu' and 'evsi', generated by `interpolate_evsi()`
    # nb:        Net benefits matrix, rows = samples, columns = treatments (1st column = 
    #            standard care).
    # c_fix:     Fixed trial cost (range or value).
    # c_var:     Monthly trial conduct cost (range or value).
    # t_lag_awr: Time between end of follow-up and change of practice for AWR (range or value).
    # t_lag_oir: Time between end of follow-up and change of practice for OIR (range or value).
    # inc_pop:   Monthly incidence of patients affected by the decision (range or value).
    # prev_pop:  Prevalent patient population affected by the decision (range or value).
    # dec_th:    Decision relevance horizon in months (range or value).
    # prob_rev:  Probability that approval can be reversed if AWR is chosen.
    # dr:        Annual discount rate.
  
  # Helper function to calculate effective decision relevance horizons
  decision_horizons <- function(fu, t_lag, dec_th) {
    dh <- sort(fu, decreasing = TRUE) + (mean(dec_th) - sum(range(fu))) - mean(t_lag)
    return(pmax(dh, 0))
  }
  
  # Calculate effective decision relevance horizons for AWR and OIR
  v_dh_awr <- decision_horizons(evsi_ar$fu, t_lag_awr, dec_th)
  v_dh_oir <- decision_horizons(evsi_ar$fu, t_lag_oir, dec_th)
  
  # Calculate population EVSI for AWR and OIR
  v_pop_evsi_awr <- (mean(prev_pop) + mean(inc_pop) * v_dh_awr) * evsi_ar$evsi * prob_rev
  v_pop_evsi_oir <- (mean(prev_pop) + mean(inc_pop) * v_dh_oir) * evsi_ar$evsi
  
  # Incremental net benefits compared to standard care
  m_inb_new <- matrix(nb[, -1] - nb[, 1]) 
  
  # Health opportunity costs if a suboptimal treatment is adopted
  c_awr <-  pmin(max(colMeans(m_inb_new)), 0) * (mean(prev_pop) + mean(inc_pop) * evsi_ar$fu)
  c_oir <-  pmax(colMeans(m_inb_new), 0) * (mean(prev_pop) + mean(inc_pop) * evsi_ar$fu)
  
  # Calculate variable trial costs and monthly discount rates
  v_c_var <- mean(c_var) * evsi_ar$fu
  v_dr    <- 1 / (1 + dr) ^ (evsi_ar$fu / 12)
  
  # Calculate ENBS for AWR
  v_enbs_awr <- (v_pop_evsi_awr - abs(c_awr) - v_c_var) * v_dr - mean(c_fix)
  
  # Calculate ENBS for OIR
  v_enbs_oir <- (v_pop_evsi_oir - abs(c_oir) - v_c_var) * v_dr - mean(c_fix)
  
  # Compile ENBS values into a dataframe
  df_enbs <- data.frame(
    time  = rep(evsi_ar$fu, times = 2),
    value = c(v_enbs_awr, v_enbs_oir),
    group = rep(c("Approval with Research", "Only in Research"), each = length(evsi_ar$fu))
  )
  
  return(df_enbs) # Return the results
  
} # Close `compute_enbs` function


# Compute the ENBS
df_enbs <- compute_enbs(
  evsi_ar   = df_evsi_ar,               # EVSI dataframe generated by `interpolate_evsi`
  nb        = m_nb,                     # Net benefit matrix generated in the PA
  c_fix     = 100000 / thresh,          # Assumed fixed cost of NICE re-appraisal
  c_var     = 0,                        # Assumed monthly variable cost, set to 0
  t_lag_awr = 12 - 9,                   # Adjusted AWR lag assumption
  t_lag_oir = 9 - 9,                    # Adjusted OIR lag assumption
  inc_pop   = 277,                      # Monthly incidence from TA650 (company submission)
  prev_pop  = 50,                       # Assumed prevalent patient population
  dec_th    = 90,                       # Assumed decision relevance horizon
  prob_rev  = 1,                        # Assumed probability of decision reversal for AWR
  dr        = 0.035                     # Discount rate based on NICE guidelines
)

# Function for plotting the ENBS over the range of follow-up times
plot_enbs <- function(enbs) {
  
  # Find the maximum ENBS for AWR and OIR and corresponding follow-up
  max_enbs <- aggregate(value ~ group, data = enbs, FUN = max)
  optimal_fu <- sapply(split(enbs, enbs$group), function(x) x$time[which.max(x$value)])
  df_max_enbs <- data.frame(group = max_enbs$group, 
                            time = optimal_fu, 
                            value = max_enbs$value,
                            vjust = ifelse(max_enbs$value == min(max_enbs$value), 1.5, -0.5))
  
  # Filter the dataframe for positive maximum ENBS values
  df_max_enbs <- df_max_enbs[df_max_enbs$value > 0, ]
  
  # Plot the ENBS values with points showing the optimal follow-up times
  ggplot(enbs, aes(x = time, y = value, color = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dotted", linewidth = 0.8) +
    geom_line() +
    geom_point(data = df_max_enbs, aes(x = time, y = value), shape = 19, size = 2) +
    geom_text(data = df_max_enbs, 
              aes(x = time, y = value, label = paste0("t = ", time), color = group, vjust = vjust), 
              size = 4, show.legend = FALSE) +
    labs(x = "Additional follow-up (months)", y = "Expected Net Benefit of Sampling") +
    theme_minimal() + 
    coord_cartesian(clip = "off") +
    theme(text=element_text(size=14), legend.position = "top", legend.title = element_blank())

} # Close `plot_enbs` function


# Plot the ENBS
plot_enbs(df_enbs)
