# Copyright (c) 2023 Mathyn Vervaart, Eline Aas, Karl Claxton, Anna Heath, Mark Strong, Nicky Welton, Torbjørn Wisløff 
# Licensed under the MIT License

rm(list=ls(all=TRUE))	# clear workspace

#####################################################################################
# Probabilistic Analysis
#####################################################################################

# note: standard care should always be the first element in the lists of OS and PFS probabilities (l_os and l_pfs)
# and the first column in the matrix  of net benefits (m_nb)

# install and/or load R package "here" and the cost-effectiveness model functions
list_of_packages <- c("here")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
library(here); source(here("R", "ce_model_functions.R"))

# Run a probabilistic analysis of size "K"
set.seed(123)         # set the seed for reproducibility
K <- 5e3              # number of simulations
l_pa <- pa_fun(K)     # run the probabilistic analysis

# OS and PFS probabilities
l_os <- l_pa$l_os    # OS probabilities (first element should be OS for standard care)
l_pfs <- l_pa$l_pfs  # PFS probabilities (first element should be PFS for standard care)

# incremental net benefits
thresh <- 3e4                                        # opportunity cost/willingness to pay threshold 
incr_cost <- l_pa$v_cost2 - l_pa$v_cost1             # incremental costs
incr_qaly <- l_pa$v_qaly2 - l_pa$v_qaly1             # incremental QALYs
m_nb <- cbind(l_pa$v_qaly1  - l_pa$v_cost1 / thresh, # net health benefits (first column should be NB for standard care)
              l_pa$v_qaly2  - l_pa$v_cost2 / thresh) 


# Probabilistic ICER
paste("ICER = £", round(mean(incr_cost) / mean(incr_qaly)), " per QALY gained", sep = "")

# Mean OS and PFS curves plot
trt_labs <- c("Sunitinib", "Pembrolizumab + Axitinib" )  # treatment labels
ncyc_y <- 52                                            # number of model cycles per year
plot_surv_fun(trt_labs, ncyc_y, l_os, l_pfs)            # plot OS and PFS

# Cost-effectiveness scatter plot
plot_ce_scatter_fun(incr_cost, incr_qaly, thresh, "£")


#####################################################################################
# Value of Information Analysis
#####################################################################################

# load data-simulation and EVSI functions
source(here("R", "data_gen_functions.R")); source(here("R", "voi_functions.R"))

#################################################
# Expected Value of (Partial) Perfect Information
#################################################

# Expected Value of Perfect Information
evpi <- mean(apply(m_nb, 1, max)) - max(colMeans(m_nb)); evpi

# Expected Value of Partial Perfect Information
pevpi_os <- surv_pevpi_fun(m_nb, l_os)
pevpi_os_pfs <- surv_pevpi_fun(m_nb, l_os, l_pfs)

##################################################
# Expected Value of Sample Information
##################################################

### Settings

# start (truncation) times for OS and PFS in months
v_atrisk_times_os1 <- subset(l_pa$m_ipd_os, treat == 1 & event == 0, select = tt)$tt
v_atrisk_times_os2 <- subset(l_pa$m_ipd_os, treat == 2 & event == 0, select = tt)$tt
v_atrisk_times_pfs1 <- subset(l_pa$m_ipd_pfs, treat == 1 & event == 0, select = tt)$tt
v_atrisk_times_pfs2 <- subset(l_pa$m_ipd_pfs, treat == 2 & event == 0, select = tt)$tt

# specify Gamma distribution hyperparameters for the monthly trial dropout rates
dropout1 <- 20; trisk1 <- sum(subset(l_pa$m_ipd_os, treat == 1, select = tt)$tt) # observed dropouts and time at risk in months for treatment 1
dropout2 <- 15; trisk2 <- sum(subset(l_pa$m_ipd_os, treat == 2, select = tt)$tt) # observed dropouts and time at risk in months for treatment 2
l_dropout <- list(c(dropout1, trisk1),
                  c(dropout2, trisk2))

# remove the number of trial drop-outs from the patients at risk assuming even spacing 
l_atrisk_times_os <- list(split_fun(v_atrisk_times_os1, dropout1),    # list of observed follow-up times for the patients at risk for OS
                   split_fun(v_atrisk_times_os2, dropout2))
l_atrisk_times_pfs <- list(split_fun(v_atrisk_times_pfs1, dropout1),  # list of observed follow-up times for the patients at risk for PFS
                    split_fun(v_atrisk_times_pfs2, dropout2))

# maximum additional follow-up time in months
max_add_fu <- 60

### Computations

# compute EVSI for OS and interpolate the EVSI estimates across different follow-up times using asymptotic regression
df_evsi_os <- evsi_os_fun(m_nb, l_os, l_atrisk_times_os, max_add_fu, ncyc_y, l_dropout, l_enroll = NULL)  
evsi_plot_fun(df_evsi_os, pevpi_os) # plot the EVSI estimates

# compute EVSI for OS + PFS and interpolate the EVSI estimates across different follow-up times using asymptotic regression
# note: the computation time for OS + PFS is about 2-3 times longer than for OS only
df_evsi_os_pfs <- evsi_os_pfs_fun(m_nb, l_os, l_pfs, l_atrisk_times_os, l_atrisk_times_pfs, max_add_fu, ncyc_y,  l_dropout, l_enroll = NULL)
evsi_plot_fun(df_evsi_os_pfs, pevpi_os_pfs) # plot the EVSI estimates


##################################################
# Expected Net Benefit of Sampling
##################################################

# Predict additional follow-up required until end of follow-up (requires total of 404 observed deaths, 156 have been observed during current follow-up) 
add_events <- 404 - 156                             # additional OS events required for the final analysis 
round(spline(df_evsi_os$os_events, df_evsi_os$time, # estimated time in months until end of follow-up
             ties = min,  xout = add_events)$y) 

# range for the fixed trial setup costs
c_fix <- c(0,0) 

# range for the monthly variable trial costs (estimated from Park et al. 2022, JAMA)
d2p <- 0.7271 * 0.7                                    # exchange rate and purchasing power parities US dollar to GBP in 2021
c_site <- 5000 * 124                                   # monthly site management costs for 124 sites
c_database <- 2500                                     # monthly database management costs
c_fu_pat <- 313 * length(unlist(l_atrisk_times_os))    # monthly follow-up costs for 670 patients at risk 
c_var_mu <- (sum(c_site, c_database, c_fu_pat) * d2p)  # mean total monthly costs
c_var <- c(c_var_mu * 0.9, c_var_mu * 1.1)             # range for the monthly costs

# range for the decision reversal costs
c_rev <- c(0, 0)

# monthly incident population
inc_pop <- ((12600 * 0.8 * 0.75 * 0.44) / 12)  # estimated from company submission in TA650
inc_pop <- c(inc_pop * 0.90, inc_pop * 1.10)   # range for the monthly incidence

# prevalent "catch-up" population
prev_pop <- c(50, 100) # assumption

# other settings
t_lag <- c(3, 6)  # range for the lag time between the end of follow-up and decision making in months
dec_th <- 60      # time horizon at t_1 in months (will be reduced by the additional follow-up + lag time)
dr_voi <- 0.035   # annual discount rate
reversal <- 1     # probability that an approval decision can be reversed

#####  compute the ENBS #####
# costs are converted to health units by dividing by <thresh>
enbs_fun(df_evsi_os, m_nb,  # replace "df_evsi_os" with "df_evsi_os_pfs" (if calculated) to compute the ENBS for OS + PFS
         c_fix = c_fix / thresh, 
         c_var = c_var / thresh,
         c_var_time = NULL,
         c_var_event = add_events,
         c_rev = c_rev / thresh,
         t_lag = t_lag,
         inc_pop = inc_pop,
         prev_pop = prev_pop,
         dec_th = dec_th, 
         dr_voi = dr_voi, 
         reversal = reversal)

