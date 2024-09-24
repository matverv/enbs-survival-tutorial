############################################################################################
# Supplementary Material: Reconstruction of Individual Patient Data (IPD) from 
# vector images (PostScript files) derived from published Kaplan-Meier curves.
#
# This script utilizes an algorithm by Liu et al. (2014) for the reconstruction 
# of survival data in Rini et al. (2019) used as supplementary material for 
# the case study analysis presented in the paper:
#
# Vervaart, M. (2024). Calculating the Expected Net Benefit of Sampling for Survival Data: 
#   A Tutorial and Case Study. Medical Decision Making.
# 
# When using this code, please cite the original paper and also:
# Liu et al. (2014). Recovering the raw data behind a non-parametric survival curve. Systematic Reviews
#
# License: This code is licensed under the MIT License.
# Developed using:
#   R: Version 4.2.2 (2022-10-31)
#   RStudio: Version 2023.12.1+402 (2023-12.1+402), RStudio, Inc.
############################################################################################

# ------------------------------------------------------------------
# # Install and/or load packages and source functions necessary for IPD reconstruction
# ------------------------------------------------------------------

list_of_packages <- c("here", "survminer")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, library, character.only = TRUE)

source(here("R", "SurvivalCurveExtraction.R"))  # Load the functions for extracting raw data from vector images


############################################################################################
# Reconstruction of overall survival (OS) data from:
#  Rini et al. (2019). Pembrolizumab plus axitinib versus sunitinib for advanced renal-cell carcinoma. NEJM.
############################################################################################

# ------------------------------------------------------------------
# Extract the coordinates of the KM curves
# ------------------------------------------------------------------

# Extract the raw data from the PostScript file containing the page with the Kaplan-Meier figure for OS
result <- ExtractCurvesFromLargeFileFunction(here("data","pembrolizumab_os.ps"), 90000, 10000)

# Read and check the lengths of the extracted curves
curve_lengths <- unlist(result[3])
summary(curve_lengths)

# X and Y coordinates
x_coords <- result[1][[1]]
y_coords <- result[2][[1]]

# Define bounds and plot the data
x_min <- min(unlist(x_coords))
x_max <- max(unlist(x_coords))
y_min <- min(unlist(y_coords))
y_max <- max(unlist(y_coords))

# Plot the full page containing the KM figure of interest
plot(c(x_min,x_max), c(y_min,y_max), col="white", 
     xlim=c(x_min,x_max), ylim=c(y_min,y_max), yaxs = "i", xaxs = "i")
for(i in 1:length(x_coords)) lines(x_coords[[i]], y_coords[[i]])

# Narrow range to focus on KM figure
x_range <- c(0, 260)
y_range <- c(-50, 5)

# Plot the coordinates of the KM figure
plot(x_range, y_range, col="white", yaxs = "i", xaxs = "i")
for(i in 1:length(x_coords)) lines(x_coords[[i]], y_coords[[i]])


# ------------------------------------------------------------------
# Use the extracted coordinates to derive the KM curve for pembrolizumab-axitinib
# ------------------------------------------------------------------

# Select the curve index corresponding to pembrolizumab-axitinib
curve_idx <- 817

# Retrieve the X and Y coordinates for the selected curve
x_coords <- result[[1]][[curve_idx]]
y_coords <- result[[2]][[curve_idx]]

# Record the X and Y coordinate of the last event based on digitization (since the axis scale is unknown)
x_last <- 18.305676855895197  # X value of the last known event (approximately, based on digitization)
y_last <- 0.8223645894001165  # Y value that corresponds to x_last

# Calculate ratios by comparing all X and Y coordinates to the last known event coordinates
x_ratios <- x_coords / x_coords[105]  # Normalize X coordinates
y_ratios <- (1 - y_last) / y_coords[105]  # Invert and normalize Y coordinates to get survival probabilities

# Apply the ratios to the last known event's time to estimate the time for all events
time_est <- x_ratios * x_last
# Compute the survival probabilities using the inverted Y coordinate ratios
surv_est <- 1 - (y_ratios * y_coords)

# Plot the KM coordinates
plot(time_est, surv_est,
     type="s",          # Step function, suitable for Kaplan-Meier plots
     lwd=2,             # Line width
     ylim = c(0, 1),     # Y-axis limit
     xlab="Time",       # X-Axis label
     ylab="Estimated Survival Probability", # Y-Axis label
     main="Overall Survival - Pembrolizumab-Axitinib")  # Title


# ------------------------------------------------------------------
# Derive the event times from the steps in the KM curve for pembrolizumab-axitinib
# Expected to yield 59 events
# ------------------------------------------------------------------

# keep only unique values (add right censoring time later)
unique_times <- unique(time_est[1:length(time_est)-1])
unique_surv_probs <- unique(surv_est) #c(unique(Survival), Survival[length(Survival)])

# number of jumps in the KM curve
pts <- length(unique_surv_probs); length(unique_surv_probs) 

# we only observe 54 unique times, while 59 events were observed in the trial
# this means there must be some times with multiple events
jumps <- c(unique_surv_probs[1:(pts-1)]) - unique_surv_probs[2:pts ] ; jumps

# Plot the jumps for visual inspection to identify potential multiple simultaneous events
plot(jumps, ylim=c(0, .012), main="Event Jumps", xlab="Event Index", ylab="Jump in Survival Probability")

# Create a data frame to make it easier to inspect and possibly manually identify the jumps with multiple events
jump_df <- data.frame(index = seq_along(jumps), jumps)

# Manual identification of events with jumps assumed to represent multiple events
double_jump_indices <- c(7, 9, 14, 37, 42, 43)

# Using unique times, find the times at which the assumed multiple events occur
multi_event_times <- unique_times[double_jump_indices + 1]

# Update the event times by including the multiple events times
# We remove the first time point assuming the first reported time is a starting point (time 0) and no event has occurred yet
event_times <- sort(c(unique_times[-1], multi_event_times))

# create a dataframe with the event times and censoring indicator
surv_df <- data.frame(times = event_times, event = rep(1,length(event_times)))
surv_df


# ------------------------------------------------------------------
# Estimate censoring times for pembrolizumab-axitinib
# ------------------------------------------------------------------

# Define time intervals in months and number of patients at risk
time_int <- c(0, 4, 8, 12, 16, 20, max(time_est))
nrisk <- c(432, 417, 378, 256, 136, 18, 0) 
interval <- 1:(length(nrisk) - 1) 

# Compute the total number of censored events in each interval
n_cens <- sapply(interval, function(x) {
  nrisk[x] - nrisk[x + 1] - sum(event_times > time_int[x] & event_times <= time_int[x + 1])
})

# Function to compute evenly spaced censoring times within each interval
ma <- function(x, n = 2) {
  if (n <= 1) return(x)
  filter(x, rep(1 / n, n), sides = 2)[-length(x)]
}

# Generate uniformly distributed censoring times for each interval
time_cens_list <- sapply(interval, function(x) {
  ma(seq(time_int[x], time_int[x+1], length.out = n_cens[x] + 1))
})
time_cens <- unlist(time_cens_list)

# Combine event and censoring times along with their corresponding event indicators
times <- c(event_times, time_cens)
event <- c(rep(1, length(event_times)), rep(0, length(time_cens)))
ipd <- data.frame(times, event, treat = 2)

# Create a survival object
fit <- survfit(Surv(times, event) ~ 1, data = ipd)

# Plot the comparison of digitized (observed) survival probabilities and reconstructed (estimated) survival probabilities
plot(unique_surv_probs[2:length(unique_surv_probs)], type = "p", col = "red", ylim = c(0, 1),
     xlab = "Index", ylab = "Survival Probability", 
     main = "Comparison of Observed and Estimated Survival Probabilities")

# Overlay the estimated Kaplan-Meier survival probabilities at the times where events occurred
lines(fit$surv[fit$n.event > 0], type = "p", col = "black", cex = 0.4)

# store the calibrated IPD for pembrolizumab-axitinib
ipd_os_pembro <- ipd


# ------------------------------------------------------------------
# Use the extracted coordinates to derive the KM curve for sunitinib
# ------------------------------------------------------------------

# Select the curve index corresponding to sunitinib
curve_idx <- 443

# Retrieve the X and Y coordinates for the selected curve
x_coords <- result[[1]][[curve_idx]]
y_coords <- result[[2]][[curve_idx]]

# Record the X and Y coordinate of the last event based on digitization (since the axis scale is unknown)
x_last <- 17.61572052401747   # X value of the last known event (approximately, based on digitization)
y_last <- 0.7291788002329644  # Y value that corresponds to x_last

# Calculate ratios by comparing all X and Y coordinates to the last known event coordinates
x_ratios <- x_coords / x_coords[169]  # Normalize X coordinates
y_ratios <- (1 - y_last) / y_coords[169]  # Invert and normalize Y coordinates to get survival probabilities

# Apply the ratios to the last known event's time to estimate the time for all events
time_est <- x_ratios * x_last
# Compute the survival probabilities using the inverted Y coordinate ratios
surv_est <- 1 - (y_ratios * y_coords)

# Plot the KM coordinates
plot(time_est, surv_est,
     type="s",          # Step function, suitable for Kaplan-Meier plots
     lwd=2,             # Line width
     ylim = c(0, 1),     # Y-axis limit
     xlab="Time",       # X-Axis label
     ylab="Estimated Survival Probability", # Y-Axis label
     main="Overall Survival - sunitinib")  # Title


# ------------------------------------------------------------------
# Derive the event times from the steps in the KM curve for sunitinib
# Expected to yield 97 events
# ------------------------------------------------------------------

# keep only unique values (add right censoring time later)
unique_times <- unique(time_est[1:length(time_est)-1])
unique_surv_probs <- unique(surv_est) #c(unique(Survival), Survival[length(Survival)])

# number of jumps in the KM curve
pts <- length(unique_surv_probs); length(unique_surv_probs) 

# we only observe 86 unique times, while 97 events were observed in the trial
# this means there must be some times with multiple events
jumps <- c(unique_surv_probs[1:(pts-1)]) - unique_surv_probs[2:pts ] ; jumps

# Plot the jumps for visual inspection to identify potential multiple simultaneous events
plot(jumps, ylim=c(0, .012), main="Event Jumps", xlab="Event Index", ylab="Jump in Survival Probability")

# Create a data frame to make it easier to inspect and possibly manually identify the jumps with multiple events
jump_df <- data.frame(index = seq_along(jumps), jumps)

# Manual identification of events with jumps assumed to represent double events
double_jump_indices <- c(5, 10, 48, 53, 65, 66, 77, 80)

# Manual identification of events with jumps assumed to represent triple events
triple_jump_indices <- c(36, 76)

# Using unique times, find the times at which the assumed double events occur
double_event_times <- unique_times[double_jump_indices + 1]

# Because these are double events, each time should be repeated once to signify two events
# (the original time, plus one additional time)
double_event_times_repeated <- rep(double_event_times, each = 1)

# Using unique times, find the times at which the assumed triple events occur
triple_event_times <- unique_times[triple_jump_indices + 1]

# Since these are triple events, each time should be repeated twice to signify three events
# (the original time, plus two additional times)
triple_event_times_repeated <- rep(triple_event_times, each = 2)

# Update the event times by including the double and triple event times
# Combine and sort the times of single, double, and triple events
event_times <- sort(c(unique_times[-1], double_event_times_repeated, triple_event_times_repeated))

# Create a data frame with the all event times and censoring indicator
# Setting the event column to 1 indicates these are all event occurrences
surv_df <- data.frame(times = event_times, event = rep(1, length(event_times)))

# Output the first few rows of surv_df to verify
surv_df

# ------------------------------------------------------------------
# Estimate censoring times for sunitinib
# ------------------------------------------------------------------

# Define time intervals in months and number of patients at risk
time_int <- c(0, 4, 8, 12, 16, 20, max(time_est))
nrisk <- c(429, 401, 341, 211, 110, 20, 0)
interval <- 1:(length(nrisk) - 1) 

# Compute the total number of censored events in each interval
n_cens <- sapply(interval, function(x) {
  nrisk[x] - nrisk[x + 1] - sum(event_times > time_int[x] & event_times <= time_int[x + 1])
})

# Function to compute evenly spaced censoring times within each interval
ma <- function(x, n = 2) {
  if (n <= 1) return(x)
  filter(x, rep(1 / n, n), sides = 2)[-length(x)]
}

# Generate uniformly distributed censoring times for each interval
time_cens <- sapply(interval, function(x) {
  ma(seq(time_int[x], time_int[x+1], length.out = n_cens[x] + 1))
})
time_cens <- unlist(time_cens)

# Combine event and censoring times along with their corresponding event indicators
times <- c(event_times, time_cens)
event <- c(rep(1, length(event_times)), rep(0, length(time_cens)))
ipd <- data.frame(times, event, treat = 1)

# Create a survival object
fit <- survfit(Surv(times, event) ~ 1, data = ipd)

# Plot the comparison of digitized (observed) survival probabilities and reconstructed (estimated) survival probabilities
plot(unique_surv_probs[2:length(unique_surv_probs)], type = "p", col = "red", ylim = c(0, 1),
     xlab = "Index", ylab = "Survival Probability", 
     main = "Comparison of Observed and Estimated Survival Probabilities")

# Overlay the estimated Kaplan-Meier survival probabilities at the times where events occurred
lines(fit$surv[fit$n.event > 0], type = "p", col = "black", cex = 0.4)

# store the calibrated IPD for pembrolizumab-axitinib
ipd_os_suni <- ipd

# combine the IPD for OS for both arms
ipd_os <- rbind(ipd_os_suni, ipd_os_pembro)

# Create a survival object and plot the Kaplan-Meier survival estimate
fit <- survfit(Surv(times, event) ~ treat, data = ipd_os)
survminer::ggsurvplot(fit,
                      data = ipd_os,
                      combine = TRUE,
                      risk.table = "nrisk_cumcensor",
                      palette = c("#BC3C29FF", "#0072B5FF"),
                      xlim = c(0,24),
                      xlab = "Months",
                      break.x.by = 4,
                      tables.height = 0.2,
                      tables.y.text = FALSE,
                      tables.theme = theme_cleantable(),
                      break.time.by = 20,
                      legend = "none",
                      title = "Overall Survival"
)

# save dataset as RData file
# save(ipd_os, file = here("data","ipd_os.RData"))




############################################################################################
# Reconstruction of progression-free survival data from:
#  Rini et al. (2019). Pembrolizumab plus axitinib versus sunitinib for advanced renal-cell carcinoma. NEJM.
############################################################################################

# ------------------------------------------------------------------
# Extract the coordinates of the KM curves
# ------------------------------------------------------------------

# Extract the raw data from the PostScript file containing the page with the Kaplan-Meier figure for PFS
result <- ExtractLinesFromLargeFileFunction(here("data","pembrolizumab_pfs.ps"), 90000, 10000)

# Read and check the lengths of the extracted curves
curve_lengths <- unlist(result[3])
summary(curve_lengths)

# X and Y coordinates
x_coords <- result[1][[1]]
y_coords <- result[2][[1]]

# Define bounds and plot the data
x_min <- min(unlist(x_coords))
x_max <- max(unlist(x_coords))
y_min <- min(unlist(y_coords))
y_max <- max(unlist(y_coords))

# Plot the full page containing the KM figure of interest
plot(c(x_min,x_max), c(y_min,y_max), col="white", 
     xlim=c(x_min,x_max), ylim=c(y_min,y_max), yaxs = "i", xaxs = "i")
for(i in 1:length(x_coords)) lines(x_coords[[i]], y_coords[[i]])

# Narrow range to focus on KM figure
x_range <- c(180,460)
y_range <- c(550,740)

# Plot the coordinates of the KM figure
plot(x_range, y_range, col="white", yaxs = "i", xaxs = "i")
for(i in 1:length(x_coords)) lines(x_coords[[i]], y_coords[[i]])


# ------------------------------------------------------------------
# Use the extracted coordinates to derive the KM curve for pembrolizumab-axitinib
# ------------------------------------------------------------------

# Define constants for curve identification and axis mapping
curve_idx <- 377      # Curve index for pembrolizumab-axitinib
x_zero_idx <- 97      # Index for x-axis 0
x_max_idx <- 103      # Index for max x-axis value
y_zero_idx <- 96      # Index for y-axis 0 (survival)
y_max_idx <- 86       # Index for max y-axis value (survival 100%)

# Time mapping parameters based on x-axis tick marks
t_max <- 24    # Max time corresponding to last x-axis tick mark
dx <- (x_coords[[x_max_idx]][1] - x_coords[[x_zero_idx]][1]) / t_max

# Convert X coordinates to estimated time points
time_est <- (x_coords[[curve_idx]] - x_coords[[x_zero_idx]][1]) / dx

# Survival mapping parameters based on y-axis tick marks
dy <- y_coords[[y_max_idx]][1] - y_coords[[y_zero_idx]][1]

# Convert Y coordinates to estimated survival probabilities
surv_est <- 1 - (y_coords[[y_max_idx]][1] - y_coords[[curve_idx]]) / dy

# Plot the Kaplan-Meier curve for pembrolizumab-axitinib
plot(time_est, surv_est,
     type = "s",    # Step function for the KM plot
     lwd = 2,       # Line width
     ylim = c(0, 1), # Y-axis limits for survival probability
     xlab = "Time", # X-axis label
     ylab = "Estimated Survival Probability", # Y-axis label
     main = "Overall Survival - Pembrolizumab-Axitinib")  # Chart title


# ------------------------------------------------------------------
# Derive the event times from the steps in the KM curve for pembrolizumab-axitinib
# Expected to yield 183 events
# ------------------------------------------------------------------

# keep only unique values (add right censoring time later)
unique_times <- unique(time_est[1:length(time_est)-1])
unique_surv_probs <- unique(surv_est) #c(unique(Survival), Survival[length(Survival)])

# number of jumps in the KM curve
pts <- length(unique_surv_probs); length(unique_surv_probs) 

# we only observe 115 unique times, while 182 events were observed in the trial
# this means there must be some times with multiple events
jumps <- c(unique_surv_probs[1:(pts-1)]) - unique_surv_probs[2:pts ] ; jumps

# Plot the jumps for visual inspection to identify potential multiple simultaneous events
plot(jumps, main="Event Jumps", xlab="Event Index", ylab="Jump in Survival Probability")

# Create a data frame to make it easier to inspect and possibly manually identify the jumps with multiple events
jump_df <- data.frame(index = seq_along(jumps), jumps)

# Manual identification of events with jumps assumed to represent 2 events
double_jump_indices <- c(13, 15,  17,  37,  42,  44,  46,  52,  53,  65,  67,  68,
                         70,  71,  84,  85,  86,  91,  93,  97,  98,  99, 102, 107, 111)
double_event_times <- unique_times[double_jump_indices + 1]
double_event_times_repeated <- rep(double_event_times, each = 1) # Repeat once for 2 events

# Manual identification of events with jumps assumed to represent 3 events
triple_jump_indices <- c(23, 27, 33, 39, 51, 75, 105)
triple_event_times <- unique_times[triple_jump_indices + 1]
triple_event_times_repeated <- rep(triple_event_times, each = 2) # Repeat twice for 3 events

# Manual identification of events with jumps assumed to represent 4 events
quad_jump_indices <- c(19)
quad_event_times <- unique_times[quad_jump_indices + 1]
quad_event_times_repeated <- rep(quad_event_times, each = 3) # Repeat 3 times for 4 events

# Manual identification of events with jumps assumed to represent 5 events
quint_jump_indices <- c(24, 26, 36)
quint_event_times <- unique_times[quint_jump_indices + 1]
quint_event_times_repeated <- rep(quint_event_times, each = 4) # Repeat 4 times for 5 events

# Manual identification of events with jumps assumed to represent 6 events
six_jump_indices <- c(50)
six_event_times <- unique_times[six_jump_indices + 1]
six_event_times_repeated <- rep(six_event_times, each = 5) # Repeat 5 times for 6 events

# Manual identification of events with jumps assumed to represent 11 events
elev_jump_indices <- c(25)
elev_event_times <- unique_times[elev_jump_indices]
elev_event_times_repeated <- rep(elev_event_times, each = 10) # Repeat 10 times for 11 events

# Combine and sort the times of single, double, triple, quad, quint, six, and elev events
# We remove the first time point assuming the first reported time is a starting point (time 0) and no event has occurred yet
event_times <- sort(c(unique_times[-1], 
                      double_event_times_repeated, 
                      triple_event_times_repeated,
                      quad_event_times_repeated, 
                      quint_event_times_repeated, 
                      six_event_times_repeated, 
                      elev_event_times_repeated))

# Create a data frame with all event times and censoring indicator
# Setting the event column to 1 indicates these are all event occurrences
all_surv_df <- data.frame(times = event_times, event = rep(1, length(event_times)))

# Output the first few rows of surv_df to verify
surv_df


# ------------------------------------------------------------------
# Estimate censoring times for pembrolizumab-axitinib
# ------------------------------------------------------------------

# Define time intervals in months and number of patients at risk
time_int <- c(0, 4, 8, 12, 16, 20, max(time_est))
nrisk <- c(432, 357, 251, 140, 42, 3, 0)
interval <- 1:(length(nrisk) - 1) 

# Compute the total number of censored events in each interval
n_cens <- sapply(interval, function(x) {
  nrisk[x] - nrisk[x + 1] - sum(event_times > time_int[x] & event_times <= time_int[x + 1])
})

# Function to compute evenly spaced censoring times within each interval
ma <- function(x, n = 2) {
  if (n <= 1) return(x)
  filter(x, rep(1 / n, n), sides = 2)[-length(x)]
}

# Generate uniformly distributed censoring times for each interval
time_cens_list <- sapply(interval, function(x) {
  ma(seq(time_int[x], time_int[x+1], length.out = n_cens[x] + 1))
})
time_cens <- unlist(time_cens_list)

# Combine event and censoring times along with their corresponding event indicators
times <- c(event_times, time_cens)
event <- c(rep(1, length(event_times)), rep(0, length(time_cens)))
ipd <- data.frame(times, event, treat = 2)

# Create a survival object
fit <- survfit(Surv(times, event) ~ 1, data = ipd)

# Plot the comparison of digitized (observed) survival probabilities and reconstructed (estimated) survival probabilities
plot(unique_surv_probs[2:length(unique_surv_probs)], type = "p", col = "red", ylim = c(0, 1),
     xlab = "Index", ylab = "Survival Probability", 
     main = "Comparison of Observed and Estimated Survival Probabilities")

# Overlay the estimated Kaplan-Meier survival probabilities at the times where events occurred
lines(fit$surv[fit$n.event > 0], type = "p", col = "black", cex = 0.4)


# ------------------------------------------------------------------
# Calibrate the censoring times for last time interval
# ------------------------------------------------------------------

# Function to combine estimated survival data frame with observed survival probabilities
# and calculate survival ratio and error between estimated and observed survival
create_all_surv_df <- function (fit, unique_surv_probs) {
  
  fit_summary <- summary(fit)
  
  # Create est_surv without repetitive summary calls
  est_surv <- data.frame(
    time = fit_summary$time, 
    n.risk = fit_summary$n.risk, 
    n.event = fit_summary$n.event,
    est_surv = fit_summary$surv
  )
  
  # Calculate the length of unique_surv_probs vector once instead of multiple times
  length_usp <- length(unique_surv_probs)
  
  # Combine est_surv and obs_surv, and calculate surv_ratio and error in a vectorized operation
  all_surv <- cbind(
    est_surv, 
    obs_surv = unique_surv_probs[2:length_usp],
    surv_ratio = unique_surv_probs[2:length_usp] / unique_surv_probs[1:(length_usp - 1)],
    error = est_surv$est_surv - unique_surv_probs[2:length_usp]
  )
  return(all_surv)
}

# combine estimated survival data frame with observed survival probabilities
# and calculate survival ratio and error between estimated and observed survival
all_surv <- create_all_surv_df(fit, unique_surv_probs)

# select time interval for the calibration
t <- 5

sub_all_surv <- subset(all_surv, all_surv$time > time_int[[t]] & all_surv$time < time_int[[t+1]] ) # subset the intervals that need to be improved
sub_all_surv$obs_surv <- sub_all_surv$obs_surv

cens_new <- rep(0,length(sub_all_surv$time)+1)
n_events_new <- c(sub_all_surv$n.event)
est_surv_new <- if(t>1) {
  c(min(subset(all_surv$est_surv, all_surv$time < time_int[[t]])), rep(0, length(sub_all_surv$time)))} else {
    c(1, sub_all_surv$obs_surv[1],  rep(0, length(sub_all_surv$time)-1))
  }

# error function to optimize
error_function <- function (x, index) {
  index <- i
  error <-
    (nrisk[t] - sum(n_events_new[0:index]) - sum(cens_new[1:index]) - x) /
    (nrisk[t] - sum(n_events_new[0:(index - 1)]) - sum(cens_new[1:index]) - x) *
    est_surv_new[index] - sub_all_surv$obs_surv[index]
  
  return(abs(error))
}


# optimize the censoring times for the selected time interval
for(i in 1:length(sub_all_surv$time)) {
  if ((n_cens[t] - sum(cens_new)) > 0) {
    out <-
      optim(
        par = 0,
        fn = error_function,
        method = "Brent",
        index = i,
        lower = 0,
        upper = (n_cens[t] - sum(cens_new))
      )
    cens_new[i + 1] <- round(out$par)
    est_surv_new[i + 1] <-
      (nrisk[t] - sum(n_events_new[0:i]) - sum(cens_new[1:i]) - cens_new[i +
                                                                           1]) / (nrisk[t] - sum(n_events_new[0:(i - 1)]) - sum(cens_new[1:i]) - cens_new[i +
                                                                                                                                                            1]) * est_surv_new[i]
  }
}
cens_new <- c(cens_new[2:length(cens_new)], n_cens[t] - sum(cens_new))

t_temp <- c(time_int[t], sub_all_surv$time, max(time_int[t+1]))
t_cens_new <- sapply(1:(length(t_temp)-1), function (x) {
  
  if(cens_new[x] > 0) {
    ma(seq(t_temp[x], t_temp[x+1], length.out = cens_new[x]+1 ))
  }
})
t_cens_new <- unlist(t_cens_new)

# replace the uniform censoring times with the optimized censoring times
time_cens_list[[t]] <- t_cens_new
time_cens <- unlist(time_cens_list)

# Combine event and censoring times along with their corresponding event indicators
times <- c(event_times, time_cens)
event <- c(rep(1, length(event_times)),rep(0, length(time_cens)))
ipd <- as.data.frame(cbind(times, event, treat = 2))

# Create a survival object
fit <- survfit(Surv(times, event) ~ 1, data = ipd)

# Plot the comparison of digitized (observed) survival probabilities and reconstructed (estimated) survival probabilities
plot(unique_surv_probs[2:length(unique_surv_probs)], type = "p", col = "red", ylim = c(0, 1),
     xlab = "Index", ylab = "Survival Probability", 
     main = "Comparison of Observed and Estimated Survival Probabilities")

# Overlay the estimated Kaplan-Meier survival probabilities at the times where events occurred
lines(fit$surv[fit$n.event > 0], type = "p", col = "black", cex = 0.4)

# store the calibrated IPD for pembrolizumab-axitinib
ipd_pfs_pembro <- ipd


# ------------------------------------------------------------------
# Use the extracted coordinates to derive the KM curve for sunitinib
# ------------------------------------------------------------------

# Define constants for curve identification and axis mapping
curve_idx <- 232      # Curve index for sunitinib
x_zero_idx <- 97      # Index for x-axis 0
x_max_idx <- 103      # Index for max x-axis value
y_zero_idx <- 96      # Index for y-axis 0 (survival)
y_max_idx <- 86       # Index for max y-axis value (survival 100%)

# Time mapping parameters based on x-axis tick marks
t_max <- 24    # Max time corresponding to last x-axis tick mark
dx <- (x_coords[[x_max_idx]][1] - x_coords[[x_zero_idx]][1]) / t_max

# Convert X coordinates to estimated time points
time_est <- (x_coords[[curve_idx]] - x_coords[[x_zero_idx]][1]) / dx

# Survival mapping parameters based on y-axis tick marks
dy <- y_coords[[y_max_idx]][1] - y_coords[[y_zero_idx]][1]

# Convert Y coordinates to estimated survival probabilities
surv_est <- 1 - (y_coords[[y_max_idx]][1] - y_coords[[curve_idx]]) / dy

# Plot the Kaplan-Meier curve for sunitinib
plot(time_est, surv_est,
     type = "s",    # Step function for the KM plot
     lwd = 2,       # Line width
     ylim = c(0, 1), # Y-axis limits for survival probability
     xlab = "Time", # X-axis label
     ylab = "Estimated Survival Probability", # Y-axis label
     main = "Overall Survival - sunitinib")  # Chart title


# ------------------------------------------------------------------
# Derive the event times from the steps in the KM curve for sunitinib
# Expected to yield 213 events
# ------------------------------------------------------------------

# keep only unique values (add right censoring time later)
unique_times <- unique(time_est[1:length(time_est)-1])
unique_surv_probs <- unique(surv_est) #c(unique(Survival), Survival[length(Survival)])

# number of jumps in the KM curve
pts <- length(unique_surv_probs); length(unique_surv_probs) 

# we only observe 132 unique times, while 213 events were observed in the trial
# this means there must be some times with multiple events
jumps <- c(unique_surv_probs[1:(pts-1)]) - unique_surv_probs[2:pts ] ; jumps

# Plot the jumps for visual inspection to identify potential multiple simultaneous events
plot(jumps, main="Event Jumps", xlab="Event Index", ylab="Jump in Survival Probability")

# Create a data frame to make it easier to inspect and possibly manually identify the jumps with multiple events
jump_df <- data.frame(index = seq_along(jumps), jumps)

# Manual identification of events with jumps assumed to represent 2 events
double_jump_indices <- c(7,  11,  14,  16,  44,  53, 64, 69,  
                         83,  88,  93,  96, 97, 104, 112, 116, 122)
double_event_times <- unique_times[double_jump_indices + 1]
double_event_times_repeated <- rep(double_event_times, each = 1) # Repeat once for 2 events

# Manual identification of events with jumps assumed to represent 3 events
triple_jump_indices <- c(33, 41, 55, 70, 106)
triple_event_times <- unique_times[triple_jump_indices + 1]
triple_event_times_repeated <- rep(triple_event_times, each = 2) # Repeat twice for 3 events

# Manual identification of events with jumps assumed to represent 4 events
quad_jump_indices <- c(31, 34, 35, 40, 57, 82, 121)
quad_event_times <- unique_times[quad_jump_indices + 1]
quad_event_times_repeated <- rep(quad_event_times, each = 3) # Repeat 3 times for 4 events

# Manual identification of events with jumps assumed to represent 7 events
sev_jump_indices <- c(56, 68)
sev_event_times <- unique_times[sev_jump_indices + 1]
sev_event_times_repeated <- rep(sev_event_times, each = 6) # Repeat 6 times for 7 events

# Manual identification of events with jumps assumed to represent 8 events
eight_jump_indices <- c(36, 37, 38)
eight_event_times <- unique_times[eight_jump_indices + 1]
eight_event_times_repeated <- rep(eight_event_times, each = 7) # Repeat 7 times for 8 events

# Combine and sort the times of single, double, triple, quad, quint, six, and elev events
# We remove the first time point assuming the first reported time is a starting point (time 0) and no event has occurred yet
event_times <- sort(c(unique_times[-1], 
                      double_event_times_repeated, 
                      triple_event_times_repeated,
                      quad_event_times_repeated, 
                      sev_event_times_repeated, 
                      eight_event_times_repeated))

# Create a data frame with all event times and censoring indicator
# Setting the event column to 1 indicates these are all event occurrences
all_surv_df <- data.frame(times = event_times, event = rep(1, length(event_times)))

# Output the first few rows of surv_df to verify
surv_df


# ------------------------------------------------------------------
# Estimate censoring times for sunitinib
# ------------------------------------------------------------------

# Define time intervals in months and number of patients at risk
time_int <- c(0, 4, 8, 12, 16, 20, max(time_est))
nrisk <- c(429, 302, 193, 89, 29, 1, 0)
interval <- 1:(length(nrisk) - 1) 

# Compute the total number of censored events in each interval
n_cens <- sapply(interval, function(x) {
  nrisk[x] - nrisk[x + 1] - sum(event_times > time_int[x] & event_times <= time_int[x + 1])
})

# Function to compute evenly spaced censoring times within each interval
ma <- function(x, n = 2) {
  if (n <= 1) return(x)
  filter(x, rep(1 / n, n), sides = 2)[-length(x)]
}

# Generate uniformly distributed censoring times for each interval
time_cens_list <- sapply(interval, function(x) {
  ma(seq(time_int[x], time_int[x+1], length.out = n_cens[x] + 1))
})
time_cens <- unlist(time_cens_list)

# Combine event and censoring times along with their corresponding event indicators
times <- c(event_times, time_cens)
event <- c(rep(1, length(event_times)), rep(0, length(time_cens)))
ipd <- data.frame(times, event)

# Create a survival object
fit <- survfit(Surv(times, event) ~ 1, data = ipd)

# Plot the comparison of digitized (observed) survival probabilities and reconstructed (estimated) survival probabilities
plot(unique_surv_probs[2:length(unique_surv_probs)], type = "p", col = "red", ylim = c(0, 1),
     xlab = "Index", ylab = "Survival Probability", 
     main = "Comparison of Observed and Estimated Survival Probabilities")

# Overlay the estimated Kaplan-Meier survival probabilities at the times where events occurred
lines(fit$surv[fit$n.event > 0], type = "p", col = "black", cex = 0.4)


# ------------------------------------------------------------------
# Calibrate the censoring times for the last interval
# ------------------------------------------------------------------

# combine estimated survival data frame with observed survival probabilities
# and calculate survival ratio and error between estimated and observed survival
all_surv <- create_all_surv_df(fit, unique_surv_probs)
tail(all_surv)

# select time interval for the calibration
t <- 5

sub_all_surv <- subset(all_surv, all_surv$time > time_int[[t]] & all_surv$time < time_int[[t+1]] ) # subset the intervals that need to be improved
sub_all_surv$obs_surv <- sub_all_surv$obs_surv

cens_new <- rep(0,length(sub_all_surv$time)+1)
n_events_new <- c(sub_all_surv$n.event)
est_surv_new <- if(t>1) {
  c(min(subset(all_surv$est_surv, all_surv$time < time_int[[t]])), rep(0, length(sub_all_surv$time)))} else {
    c(1, sub_all_surv$obs_surv[1],  rep(0, length(sub_all_surv$time)-1))
  }

# optimize censoring
for(i in 1:length(sub_all_surv$time)) {
  if ((n_cens[t] - sum(cens_new)) > 0) {
    out <-
      optim(
        par = 0,
        fn = error_function,
        method = "Brent",
        index = i,
        lower = 0,
        upper = (n_cens[t] - sum(cens_new))
      )
    cens_new[i + 1] <- round(out$par)
    est_surv_new[i + 1] <-
      (nrisk[t] - sum(n_events_new[0:i]) - sum(cens_new[1:i]) - cens_new[i +
                                                                           1]) / (nrisk[t] - sum(n_events_new[0:(i - 1)]) - sum(cens_new[1:i]) - cens_new[i +
                                                                                                                                                            1]) * est_surv_new[i]
  }
}
cens_new <- c(cens_new[2:length(cens_new)], n_cens[t] - sum(cens_new))

t_temp <- c(time_int[t], sub_all_surv$time, max(time_int[t+1]))
t_cens_new <- sapply(1:(length(t_temp)-1), function (x) {
  
  if(cens_new[x] > 0) {
    ma(seq(t_temp[x], t_temp[x+1], length.out = cens_new[x]+1 ))
  }
})
t_cens_new <- unlist(t_cens_new)

# replace the uniform censoring times with the optimized censoring times
time_cens_list[[t]] <- t_cens_new
time_cens <- unlist(time_cens_list)

# create dataframe with time and event indicator
times <- c(event_times, time_cens)
event <- c(rep(1, length(event_times)),rep(0, length(time_cens)))
ipd <- as.data.frame(cbind(times, event, treat = 1))

# Create a survival object
fit <- survfit(Surv(times, event) ~ 1, data = ipd)

# Plot the comparison of digitized (observed) survival probabilities and reconstructed (estimated) survival probabilities
plot(unique_surv_probs[2:length(unique_surv_probs)], type = "p", col = "red", ylim = c(0, 1),
     xlab = "Index", ylab = "Survival Probability", 
     main = "Comparison of Observed and Estimated Survival Probabilities")

# Overlay the estimated Kaplan-Meier survival probabilities at the times where events occurred
lines(fit$surv[fit$n.event > 0], type = "p", col = "black", cex = 0.4)

# store the calibrated IPD for sunitinib
ipd_pfs_suni <- ipd

# combine the IPD for OS for both arms
ipd_pfs <- rbind(ipd_pfs_suni, ipd_pfs_pembro)

# Create a survival object and plot the Kaplan-Meier survival estimate
fit <- survfit(Surv(times, event) ~ treat, data = ipd_pfs)
survminer::ggsurvplot(fit,
                      data = ipd_pfs,
                      combine = TRUE,
                      risk.table = "nrisk_cumcensor",
                      palette = c("#BC3C29FF", "#0072B5FF"),
                      xlim = c(0,24),
                      xlab = "Months",
                      break.x.by = 4,
                      tables.height = 0.2,
                      tables.y.text = FALSE,
                      tables.theme = theme_cleantable(),
                      break.time.by = 20,
                      legend = "none",
                      title = "Overall Survival"
)

# save dataset as RData file
# save(ipd_os, file = here("data","ipd_os.RData"))
