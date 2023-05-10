# Copyright (c) 2023 Mathyn Vervaart
# Licensed under the MIT License

#####################################################################################
# Install and load packages 
#####################################################################################
list_of_packages <- c("here", "flexsurv", "survminer", "MASS", "ggplot2")
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(list_of_packages, library, character.only = TRUE)


#####################################################################################
# Function for selecting ggplot colors
#####################################################################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col <- gg_color_hue(3) # set color scheme


#####################################################################################
# Function for plotting K-M data
#####################################################################################
plot_km_fun <- function (dataset, strategies, x_max, break_x, y_lab = "Survival probability", x_lab = "Time", ...) {
  
  fit_obj <- do.call(survfit, list(formula = Surv(tt, event) ~ treat, data = dataset))
  
  km_plot <- ggsurvplot(
    fit_obj,
    combine = T,
    conf.int = F,
    #title = title,
    pval = F,
    font.submain = c(font_size, "plain"),
    palette = "lancet",
    xlab = x_lab,
    ylab = y_lab,
    xlim = c(0, x_max),
    #ylim = c(0.7,1),
    risk.table = TRUE,
    break.time.by = break_x,
    
    legend = "none", 
    legend.labs = strategies,
    legend.title = "",
    
    ggtheme = theme_light(),
    tables.theme = theme_cleantable() + theme(panel.border = element_blank()),
    tables.y.text = F,
    risk.table.height = 0.2,
    
    font.x = c(font_size, "plain"),
    font.y = c(font_size, "plain"),
    font.legend = c(font_size-2, "plain"),
    font.tickslab = c(font_size-2, "plain", "black"),
    fontsize = 4
  )
  
  km_plot$plot <- km_plot$plot+ 
    ggplot2::annotate("text", 
                      x = (max(fit_obj[1]$time)) - 0.1, y = (min(fit_obj[1]$surv)+(0.042+0.006 )), # x and y coordinates of the text
                      label = strategies[1] , size = 4)
  km_plot$plot <- km_plot$plot+ 
    ggplot2::annotate("text", 
                      x = (max(fit_obj[2]$time)) - 0.1, y = (min(fit_obj[2]$surv)-0.042+0.006), # x and y coordinates of the text
                      label = strategies[2] , size = 4)
  
  return(km_plot)
}  


#####################################################################################
# Function for extracting the parameters of a beta distribution from mean and st. deviation
# Taken from: https://darthworkgroup.com/
#####################################################################################
betaPar <- function(m, s)  # 
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

#####################################################################################
# Function for estimating the standard error from the upper and lower 95% confidence interval
# Taken from: https://darthworkgroup.com/
#####################################################################################
Srange <- function(low, high)  
{
  s = (high - low) / (2 * qnorm(0.975))
  list(s = s)
}

#####################################################################################
# Function for extracting the parameters of a gamma distribution from mean and st. deviation
# Taken from: https://darthworkgroup.com/
#####################################################################################
gammaPar <- function(mu, sigma) {   
  # mu: mean  
  # sigma: standard deviation 
  shape <- mu ^ 2 / sigma ^ 2
  scale <- sigma ^ 2 / mu
  list(shape = shape, scale = scale)
}

##############################################################################################
# function computing meanlog and log-sd using method of moments
# https://devinincerti.com/2018/02/10/psa.html#gamma-and-lognormal-distributions
##############################################################################################
lnorm_mom <- function(mean, sd){
  if (mean > 0){
    sigma2 <- log((sd^2 + mean^2)/mean^2)
    mu <- log(mean) - 1/2 * sigma2
  } else{
    stop("Mean must be positive")
  }
  
  return(list(mu = mu, sigma2 = sigma2))
}

#####################################################################################
# Function to fit survival models
#####################################################################################
survfit_fun <- function (dists, 
                         ipd, 
                         t_horizon) {  
  
  
  # create lists to store model summaries
  surv_summ <- vector(mode = "list", length = 0)
  surv_summ$fit <- vector(mode = "list", length = length(dists)) # list of models fits
  surv_summ$aic <- matrix(NaN, length(dists), 1) # vectors of aic scores
  surv_summ$Aw <- matrix(NaN, length(dists), 1) # vectors of aic weights
  surv_summ$mean_survival <- matrix(NaN, length(dists), 1) # mean survival estimates
  surv_summ$events <- matrix(NaN, length(dists), 1) # event numbers
  surv_summ$trisk <- matrix(NaN, length(dists), 1) # time at risk
  
  names(surv_summ$fit) <-
    names(surv_summ$aic) <-
    names(surv_summ$Aw) <-
    names(surv_summ$mean_survival) <-
    names(surv_summ$events) <-
    names(surv_summ$trisk) <- c(dists)
  
  # fit the models
  for (i in 1:length(dists)) {
    
    if(dists[i] == "exp" || dists[i] == "gengamma") {
      
      model <- flexsurvreg(Surv(tt, event) ~ 1, data = ipd, 
                           dist = dists[i], method = "BFGS")
      
    } else {
      
      model <- flexsurvreg(Surv(tt, event) ~ 1, data = ipd, 
                           dist = dists[i], method = "Nelder-Mead")
      
    }
    
    surv_summ$fit[[i]] <- cbind(model$coefficients, model$cov)
    
    surv_summ$aic[[i]] <- model$AIC
    
    surv_summ$mean_survival[[i]] <- 
      eval(parse(text = paste0("rmst_", dists[i], "(", t_horizon, ",", 
                               paste0(as.vector(model$res[,1]), collapse=","), ")")  )) 
    
    surv_summ$events[[i]] <- model$events
    
    surv_summ$trisk[[i]] <- model$trisk
    
  }
  
  # compute AIC-based model weights
  aic_min <- min(surv_summ$aic) # lapply(surv_summ, function (x) {min(x$aic)} ) # identify the lowest AIC for each dataset
  
  # transform the AIC differences back to the scale of probabilities
  Ak <- exp(-0.5*(surv_summ$aic-aic_min))
  
  # compute weights
  surv_summ$Aw <- apply(Ak, 2, function (y) {
    Ak / sum(y)
  })
  
  names(surv_summ$Aw) <- c(dists)
  
  # function output
  return(surv_summ)
}


#####################################################################################
# Function for sampling from multiple prior multivariate normal distributions
#####################################################################################
mvrnorm_ma_fun <- function (surv_summ, 
                            v_cycle, 
                            n_sim) {
  
  # sample K survival distributions
  dists <- names(surv_summ$fit)
  dist_select <- sample(x = dists, n_sim, replace = T, prob = surv_summ$Aw)
  
  # matrix of transformed prior mean vectors
  mu <- lapply(dist_select, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,1]"))) 
  })
  
  # matrix of transformed prior variance matrices
  cov_matrix <- lapply(dist_select, function (x) { 
    eval(parse(text = paste0("surv_summ$fit$", x, "[,2:","ncol(surv_summ$fit$", x,")]")))
  })
  
  # sample from priors on transformed scale
  prior_samples_trans <- lapply(1:n_sim, function (x) {
    mvrnorm(1, mu = mu[[x]], Sigma = cov_matrix[[x]])
  })
  
  # transform parameter samples back to original scale
  prior_samples <- lapply(1:n_sim, function (x) {
    sapply(1:length(prior_samples_trans[[x]]), function (y)  {
      eval(parse(text = paste0("flexsurv.dists$", dist_select[[x]],
                               "$inv.transforms[[y]](prior_samples_trans[[x]][y])")))
    })
  })
  
  # compute prior survival
  prior_surv <- sapply(1:n_sim, function (x) {
    pdist_parse <- parse(text = paste0("p", dist_select[[x]], "(", "y", 
                                       ",", paste(prior_samples[[x]], collapse = ",", sep = ","), ", lower.tail = F)"))
    
    pdist <- function (y) {eval(pdist_parse)}
    pdist(v_cycle)
  })
  
  
  # return output
  output <- list(distribution = dist_select, theta = prior_samples, prior_surv = prior_surv)
  return(output)
}


#####################################################################################
# Function for sampling from a single multivariate normal distribution
#####################################################################################
mvrnorm_fun <- function (cov_matrix,
                         dist,
                         n_sim) {
  
  # matrix of transformed prior mean vectors
  mu <- cov_matrix[,1]
  
  # matrix of transformed prior variance matrices
  covar <- cov_matrix[,2:ncol(cov_matrix)]
  
  # sample from priors on transformed scale
  prior_samples_trans <- mvrnorm(n_sim, mu = mu, Sigma = covar)
  
  # transform parameter samples back to original scale
  prior_samples <- sapply(1:length(prior_samples_trans[1,]), function (x)  {
    eval(parse(text = paste0("flexsurv.dists$", dist,
                             "$inv.transforms[[x]](prior_samples_trans[,x])")))
  })
  prior_samples <- lapply(1:nrow(prior_samples), function (x) prior_samples[x,])
  
  # compute prior survival
  prior_surv <- sapply(1:n_sim, function (x) {
    pdist_parse <- parse(text = paste0("p", dist, "(", "y",
                                       ",", paste(prior_samples[[x]], collapse = ",", sep = ","), ", lower.tail = F)"))
    
    pdist <- function (y) {eval(pdist_parse)}
    pdist(v_cycle)
  })
  
  
  # return output
  return(prior_surv)
}



#####################################################################################
# Ara and Brazier 2010 age-based utility decrement function. Company used starting age of 73
#####################################################################################
u_age_decr_fun <- function(male, age) {
  (0.9508566 + 0.0212126*male - 0.0002587*age - 0.0000332*age^2) - (0.9508566 + 0.0212126*male - 0.0002587*(age + v_cycle_year) - 0.0000332*(age + v_cycle_year)^2)
}

##############################################################################################
# function for applying waning of the treatment effect
##############################################################################################
wane_fun <- function (x, y, wane_start, wane_end) {
  wane_prop <- seq.int(0, 1, length.out = (wane_end - wane_start) +1)
  x[(wane_start):(wane_end)] <- (1-wane_prop) * x[(wane_start):(wane_end)] + wane_prop * y[(wane_start):(wane_end)]
  x[(1:t_horizon) > wane_end] <- y[(1:t_horizon) > wane_end]
  return(x)
}

#####################################################################################
# Function for plotting expected OS and PFS survival curves (1 figure)
#####################################################################################
plot_surv_fun <- function (trt_labs, ncyc_y, l_os, l_pfs = NULL) { 
  
  # create lists and vectors of survival probabilities
  if(!is.null(l_pfs)) {
    v_mean_surv <- c(as.vector(sapply(l_os, rowMeans)),
                     as.vector(sapply(l_pfs, rowMeans))
    )
    l_surv <- c(l_os, l_pfs)
  } else {
    v_mean_surv <- as.vector(sapply(l_os, rowMeans))
    l_surv <- l_os
  }
  
  # upper and lower bounds
  ul <- sapply(l_surv, function (x) {
    apply(x, 1, function (y) quantile(y, 0.975))
  })
  ul <- as.vector(ul)
  ll <- sapply(l_surv, function (x) {
    apply(x, 1, function (y) quantile(y, 0.025))
  })
  ll <- as.vector(ll)
  
  # create indices
  if(!is.null(l_pfs)) {
    curve_index <- rep(c(" Overall survival.", " Progression-free survival."), each = (length(l_os) * nrow(l_os[[1]])))
    group_index <- rep(1:length(l_os), each = nrow(l_os[[1]]), times = 2)
    group_index <- trt_labs[group_index]
    time_index <- rep(0:(nrow(l_os[[1]])-1), length(l_os) * 2) / (ncyc_y / 12)
  } else {
    curve_index <- rep(" Overall survival.", length = (length(l_os) * nrow(l_os[[1]])))
    group_index <- rep(1:length(l_os), each = nrow(l_os[[1]]))
    group_index <- trt_labs[group_index]
    time_index <- rep(0:(nrow(l_os[[1]])-1), length(l_os)) / (ncyc_y / 12)
  }
  
  df_temp <- data.frame(
    "time" = time_index, 
    "surv" = v_mean_surv,
    "group" = as.factor(group_index),
    "curve" = as.factor(curve_index),
    "ul" = ul,
    "ll" = ll
  )
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  col <- gg_color_hue(2) 
  col <- rep(col, length(unique(df_temp$curve)))
  
  override.linetype <- rep(1:length(l_os), each = length(unique(df_temp$curve))) # c("solid", "dashed", "solid", "dashed")
  
  # Create the ggplot line graph
  ggplot(df_temp, aes(x = time, y = surv, color = interaction(group, curve), linetype = curve)) +
    theme_light() +
    geom_line() +
    scale_color_manual(values = c(col)) +
    #scale_linetype_manual(values =override.linetype) +
    labs(x = "Months", y = "Survival Probability", 
         color = "", linetype = "") +
    guides(colour = guide_legend(nrow=length(l_surv)/2,byrow=TRUE, override.aes = list(linetype = override.linetype))) +
    scale_linetype(guide = "none") +
    theme(legend.position = "top", legend.title = element_blank()) +
    geom_ribbon(aes(ymin=ll,ymax=ul, fill = group), alpha=0.2,  colour = NA,  show.legend = F) +
    theme(axis.text.x = element_text(color="black")) +
    theme(axis.text.y = element_text(color="black")) +
    theme(text=element_text(size=18),  legend.text = element_text(size = 9))
}    


#####################################################################################
# Function for creating a cost-effectiveness scatter plot
#####################################################################################
plot_ce_scatter_fun <- function (incr_cost, incr_qaly, thresh, val = "£") {
  
  CE <- data.frame("Cost" = incr_cost, "QALY" =  incr_qaly)
  
  ggplot(data=CE, aes(x=QALY, y=Cost)) + 
    #geom_hline(yintercept = 0, color = "darkgrey", size = 0.8) + 
    #geom_vline(xintercept = 0, color = "darkgrey", size = 0.8) +
    xlab("Incremental QALY") +
    ylab(paste0("Incremental Cost ", "(", val, ")")) +
    theme_light() + 
    geom_point(aes(color = col[1]), size = 0.7, alpha = 0.5) +
    geom_point(aes(x=mean(QALY), y=mean(Cost)) ,colour="black", size = 2) + 
    geom_abline(intercept =  0, slope = thresh, color = "black", linetype = "dashed", linewidth = 1) +
    stat_ellipse() +
    theme(legend.position = "none") +
    ggtitle("") + scale_size(range=c(0.1, 2), guide="none") + 
    #ylim(c(0,45000)) + 
    #scale_y_continuous(breaks = limits = ) +
    theme(axis.text.x = element_text(color="black")) +
    theme(axis.text.y = element_text(color="black")) +
    theme(text=element_text(size=18))
  #+
  #theme(text = element_text(size = 12)) 
}




#####################################################################################
# General model information 
#####################################################################################
strategies <- c("Pembrolizumab - axitinib", "Sunitinib")  # strategy names 
# n_sim <- nsim # number of simulations
#psa <- 0 # "1" is probabilistic (not in use)

# time horizon
year2week <- 365.25/7 # year to week conversion
month2year <- 12 # month to year conversion
t_horizon_years <- 40 # time horizon in years
t_horizon <- round((t_horizon_years*year2week)) # time horizon in weeks

# model cycles
v_cycle <- 0:t_horizon # vector of model cycles in weeks
v_cycle_year <- seq.int(0, t_horizon_years, (1/year2week)) # vector of model cycles in years
n_cycles <- length(v_cycle)

# discounting
v_dr_c <- 1/(1+0.035)^v_cycle_year # weekly discount rate for effects
v_dr_e <- 1/(1+0.035)^v_cycle_year # weekly discount rate for costs

# patient information
age <- 61.5
male <- 0.735
weight <- 81.7

# opportunity cost/willingness to pay threshold 
#threshold <- 3e4

# assumed rebates to achieve ICER below 30,000 per QALY gained
rebate_pembro <-  0.85
rebate_axi <-  0.85
rebate_suni <-  0.85

# treatment effect waning 
waning <- TRUE
wane_start <- round(2*year2week,0)  # mean point at which waning of treatment effect starts
wane_end <- round(5*year2week,0)    # mean point at which waning of treatment effect ends (expected hazards are equal between arms from this point)

font_size <- 14


#####################################################################################
# run the probabilistic analysis for pembrolizumab+axitinib vs. sunitnib
#####################################################################################
pa_fun <- function(nsim, price = 1) {
  
  n_sim <- nsim # number of simulations
  
  #####################################################################################
  # Data sources
  #####################################################################################
  # Company submission, NICE committee papers: https://www.nice.org.uk/guidance/ta650/documents/committee-papers
  # Norwegian Medicines Agency STA report: https://legemiddelverket.no/Documents/Offentlig%20finansiering%20og%20pris/Metodevurderinger/K/Keytruda+Inlyta_%20nyrecellekarsinom_2020.pdf 
  # Bensimon et al. 2020: https://pubmed.ncbi.nlm.nih.gov/32697113/
  # Rini et al. 2019. Pembrolizumab plus Axitinib versus Sunitinib for Advanced Renal-Cell Carcinoma. New England Journal of Medicine [Internet]. Available from: https://doi.org/10.1056/NEJMoa1816714
  # UK Office for National Statistics: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/singleyearlifetablesuk1980to2018
  
  #####################################################################################
  # Overall survival, progression free survival and time on treatment 
  #####################################################################################
  
  #################################
  # Overall survival
  #################################
  
  # load OS data
  load(here("data", "l_pembrolizumab_os.RData"))
  
  # plot KM data
  #km_os <- plot_km_fun(ipd_os,  c("Pembrolizumab + Axitinib", "Sunitinib"), x_max = 24, break_x = 4)
  #km_os
  
  # convert times from months to weeks
  ipd_os$tt <- ipd_os$tt/month2year*year2week
  
  # choose survival models
  dist_os_pembro <- c("llogis", "exp")
  dist_os_suni <- c("exp")
  
  # fit survival models
  smod_os_pembro <- survfit_fun(dist_os_pembro, subset(ipd_os, ipd_os$treat == 1), t_horizon)
  smod_os_suni <- survfit_fun(dist_os_suni, subset(ipd_os, ipd_os$treat == 2), t_horizon)
  
  # sample from joint parameter distributions and compute the proportion surviving at each model cycle
  os_pembro_list <- mvrnorm_ma_fun(smod_os_pembro, v_cycle, n_sim)
  os_suni_list <- mvrnorm_ma_fun(smod_os_suni, v_cycle, n_sim)
  
  S_os_pembro <- os_pembro_list$prior_surv
  S_os_suni <- os_suni_list$prior_surv
  
  ### constrain OS curves by general population survival survival
  
  # compute hazards
  h_os_pembro <- apply(S_os_pembro, 2, function (x) diff(-log(x)))
  h_os_suni <- apply(S_os_suni, 2, function (x) diff(-log(x)))
  
  # GP mortality per model cycle
  gp_mort_rate <- read.csv(here("data","mort_rate_weekly_uk_2019.csv"))  # load lifetable data
  gp_mort_rate_df <- data.frame(age = gp_mort_rate$Age, mort_r = gp_mort_rate$mort_r_f*(1-male)+gp_mort_rate$mort_r_m*male) # compute average conditional survival given proportion men/women
  #gp_mort_rate_df <- rbind(gp_mort_rate_df, cbind(age = seq.int(101, 110, 1), mort_r = rep(gp_mort_rate_df$mort_r[nrow(gp_mort_rate_df)], 10)) )
  # age_cycle <- data.frame(age = floor(age + v_cycle_year[1:(length(v_cycle_year)-1)]))
  # gp_mort_rate_df <- merge(gp_mort_rate_df, age_cycle, by = "age", sort = FALSE)
  gp_mort_rate_df <- gp_mort_rate_df[gp_mort_rate_df$age>=floor(age),]
  age_cycle <- v_cycle_year[1:(length(v_cycle_year)-1)] + age
  gp_mort_r <- spline(gp_mort_rate_df$age, gp_mort_rate_df$mort_r, xout = age_cycle, ties = max)$y
  
  # adjust the hazards for GP mortality hazards
  h_os_pembro <- apply(h_os_pembro, 2, function (x) pmax(x, gp_mort_r))
  h_os_suni <- apply(h_os_suni, 2, function (x) pmax(x, gp_mort_r))
  
  # compute adjusted survival curve probabilities
  S_os_pembro <- apply(h_os_pembro, 2, function (x) c(1,exp(-cumsum(x))))
  S_os_suni <- apply(h_os_suni, 2, function (x) c(1,exp(-cumsum(x))))
  
  
  #################################
  # Progression-free survival
  #################################
  
  # load pfs data
  load(here("data", "l_pembrolizumab_pfs.RData"))
  
  # plot KM data
  #km_pfs <- plot_km_fun(ipd_pfs,  c("Pembrolizumab + Axitinib", "Sunitinib"), x_max = 24, break_x = 4)
  #km_pfs
  
  # convert times from months to weeks
  ipd_pfs$tt <- ipd_pfs$tt/month2year*year2week
  
  # subset the PFS data for each treatment
  ipd_pfs_pembro <- subset(ipd_pfs, ipd_pfs$treat == 1)
  ipd_pfs_suni <- subset(ipd_pfs, ipd_pfs$treat == 2)
  
  # bootstrap the row numbers of the PFS datasets
  samp_mat_pembro <- t(matrix(sample(1:length(ipd_pfs_pembro$tt), length(ipd_pfs_pembro$tt)*n_sim, replace=T),
                              ncol=length(ipd_pfs_pembro$tt), byrow=T))
  
  samp_mat_suni <- t(matrix(sample(1:length(ipd_pfs_suni$tt), length(ipd_pfs_suni$tt)*n_sim, replace=T),
                            ncol=length(ipd_pfs_suni$tt), byrow=T))
  
  # generate bootstrapped KM estimates for week 0:13
  km_pfs_pembro <- apply(samp_mat_pembro, 2, function (x) {
    summary(survfit(formula = Surv(tt, event) ~ 1, data = ipd_pfs_pembro[x,]), time = c(0:13))$surv
  })
  
  km_pfs_suni <- apply(samp_mat_suni, 2, function (x) {
    summary(survfit(formula = Surv(tt, event) ~ 1, data = ipd_pfs_suni[x,]), time = c(0:13))$surv
  })
  
  # choose survival models
  dist_pfs_pembro <- c("exp")
  dist_pfs_suni <- c("exp")
  
  # fit survival models 
  smod_pfs_pembro <- survfit_fun(dist_pfs_pembro, subset(ipd_pfs, ipd_pfs$treat == 1), t_horizon)
  smod_pfs_suni <- survfit_fun(dist_pfs_suni, subset(ipd_pfs, ipd_pfs$treat == 2), t_horizon)
  
  # sample from joint parameter distributions and compute the proportion surviving at each model cycle
  pfs_pembro_list <- mvrnorm_ma_fun(smod_pfs_pembro, v_cycle, n_sim)
  pfs_suni_list <- mvrnorm_ma_fun(smod_pfs_suni, v_cycle, n_sim)
  
  S_pfs_pembro <- pfs_pembro_list$prior_surv
  S_pfs_suni <- pfs_suni_list$prior_surv
  
  # compute piecewise PFS curves
  S_pfs_pembro <- apply(S_pfs_pembro, 2, function (x) x[2:(n_cycles-13)])
  S_pfs_suni <- apply(S_pfs_suni, 2, function (x) x[2:(n_cycles-13)])
  
  S_pfs_pembro <- sapply(1:n_sim, function (x) {
    S_pfs_pembro[,x] * min(km_pfs_pembro[,x])
  })
  S_pfs_suni <- sapply(1:n_sim, function (x) {
    S_pfs_suni[,x] * min(km_pfs_suni[,x])
  })
  
  S_pfs_pembro <- rbind(km_pfs_pembro, S_pfs_pembro)
  S_pfs_suni <- rbind(km_pfs_suni, S_pfs_suni)
  
  ### constrain PFS curves by OS curves
  S_pfs_pembro <- replace(S_pfs_pembro, S_pfs_pembro > S_os_pembro, S_os_pembro[S_pfs_pembro > S_os_pembro])
  S_pfs_suni <- replace(S_pfs_suni, S_pfs_suni > S_os_suni, S_os_suni[S_pfs_suni > S_os_suni])
  
  
  
  #################################
  # Time on treatment
  #################################
  
  # fit Weibull model to PFS data in order to estimate correlation between shape and scale
  smod_pfs_pembro <- survfit_fun("weibull", subset(ipd_pfs, ipd_os$treat == 1), t_horizon)
  
  # covariance matrix for log Weibull parameters for pembrolizumab ToT, assuming same correlation as for PFS
  tot_pembro_weib_cvar <- matrix(rbind(c(log(1/exp(0.2463)), Srange(0.1209, 0.3716)$s, smod_pfs_pembro$fit$weibull[1,3]),
                                       c(4.6185,  smod_pfs_pembro$fit$weibull[2,2], Srange(4.4190, 4.8181)$s)),  
                                 nrow = 2, ncol = 3) 
  
  # covariance matrix for log exponential parameters for axitinib ToT
  tot_axi_exp_cvar <- matrix(c(log(0.0109),  Srange(log(0.0094),log(0.0125))$s), 
                             nrow = 1, ncol = 2)
  
  # covariance matrix for log exponential parameters for sunitinib ToT
  tot_suni_exp_cvar <- matrix(c(log(0.0155),  Srange(log(0.0135),log(0.0174))$s),
                              nrow = 1, ncol = 2)
  
  # sample from joint parameter distributions and compute the proportion surviving at each model cycle
  S_tot_pembro <- mvrnorm_fun(tot_pembro_weib_cvar, "weibull", n_sim)
  S_tot_axi <- mvrnorm_fun(tot_axi_exp_cvar, "exp", n_sim)
  S_tot_suni <- mvrnorm_fun(tot_suni_exp_cvar, "exp", n_sim)
  
  # treatment stopping rule (2-year max) for pembrolizumab
  S_tot_pembro <- replace(S_tot_pembro, v_cycle > 103, 0)
  
  ### constrain ToT curves by OS curves
  S_tot_pembro <- replace(S_tot_pembro, S_tot_pembro > S_os_pembro, S_os_pembro[S_tot_pembro > S_os_pembro])
  S_tot_axi <- replace(S_tot_axi, S_tot_axi > S_os_pembro, S_os_pembro[S_tot_axi > S_os_pembro])
  S_tot_suni <- replace(S_tot_suni, S_tot_suni > S_os_suni, S_os_suni[S_tot_suni > S_os_suni])
  
  
  # Plot mean OS and PFS
  #plot_surv_fun(S_os_pembro, S_os_suni, S_pfs_pembro, S_pfs_suni)
  
  
  
  #####################################################################################
  # Treatment effect waning
  #####################################################################################
  
  
  #################################
  # Additive hazard approach
  #################################
  
  ### Overall survival
  if (waning == TRUE) {
    
    # specify lognormal distributions for start and end times of the waning period
    wane_start_par <- lnorm_mom(wane_start, wane_start/2) # start of waning period
    wane_end_par <- lnorm_mom(wane_end-wane_start, wane_end-wane_start) # duration of waning period
    
    # sample waning start times
    wane_start_runif <- runif(n_sim,
                              0,
                              plnorm(t_horizon, wane_start_par$mu, sdlog = sqrt(wane_start_par$sigma2))
    )
    wane_start_times <- round(qlnorm(wane_start_runif, wane_start_par$mu, sqrt(wane_start_par$sigma2)))
    
    # sample waning end times
    wane_end_runif <- runif(n_sim,
                            0,
                            plnorm(t_horizon, wane_end_par$mu, sdlog = sqrt(wane_end_par$sigma2))
    )
    wane_end_times <- wane_start_times + round(qlnorm(wane_end_runif, wane_end_par$mu, sqrt(wane_end_par$sigma2)))
    wane_end_times[wane_end_times>t_horizon] <- t_horizon
    
    
    ### Overall survival
    
    # compute mean survival
    mean_os_pembro <- rowMeans(S_os_pembro)
    mean_os_suni <- rowMeans(S_os_suni)
    
    # compute the difference in hazards for mean survival between the treatment arms
    haz_mean_os_pembro <- -log(1-(1-(mean_os_pembro[2:(t_horizon+1)]/mean_os_pembro[1:t_horizon])))
    haz_mean_os_suni <- -log(1-(1-(mean_os_suni[2:(t_horizon+1)]/mean_os_suni[1:t_horizon])))
    diff_haz_mean_os <- haz_mean_os_suni - haz_mean_os_pembro
    
    # compute the hazards for each treatment arm
    haz_os_pembro <- apply(S_os_pembro, 2, function (x) -log(1-(1-(x[2:(t_horizon+1)]/x[1:t_horizon]))))
    haz_os_suni <- apply(S_os_suni, 2, function (x) -log(1-(1-(x[2:(t_horizon+1)]/x[1:t_horizon]))))
    
    # apply treatment effect waning to each sampled vector of hazards
    S_os_pembro_wane <- matrix(nrow = nrow(S_os_pembro), ncol = ncol(S_os_pembro))
    for (i in 1:n_sim) {
      
      # apply waning function
      h_wane_os_pembro <- haz_os_pembro[,i] + punif((0:(n_cycles-2) - wane_start_times[i]) / (wane_end_times[i] - wane_start_times[i]), 0, 1) * pmax(0, diff_haz_mean_os)
      S_os_pembro_wane[,i] <- cumprod(c(1, exp(-h_wane_os_pembro)))
    }
    
    
    ### Progression-free survival
    
    # constrain PFS curves by OS curves
    S_pfs_pembro1 <- replace(S_pfs_pembro, S_pfs_pembro > S_os_pembro_wane, S_os_pembro_wane[S_pfs_pembro > S_os_pembro_wane])
    
    # compute mean survival
    mean_pfs_pembro <- rowMeans(S_pfs_pembro)
    mean_pfs_suni <- rowMeans(S_pfs_suni)
    
    # compute the difference in hazards for mean survival between the treatment arms
    haz_mean_pfs_pembro <- -log(1-(1-(mean_pfs_pembro[2:(t_horizon+1)]/mean_pfs_pembro[1:t_horizon])))
    haz_mean_pfs_suni <- -log(1-(1-(mean_pfs_suni[2:(t_horizon+1)]/mean_pfs_suni[1:t_horizon])))
    diff_haz_mean_pfs <- haz_mean_pfs_suni - haz_mean_pfs_pembro
    
    # compute the hazards for each treatment arm
    haz_pfs_pembro <- apply(S_pfs_pembro, 2, function (x) -log(1-(1-(x[2:(t_horizon+1)]/x[1:t_horizon]))))
    haz_pfs_suni <- apply(S_pfs_suni, 2, function (x) -log(1-(1-(x[2:(t_horizon+1)]/x[1:t_horizon]))))
    
    # apply treatment effect waning to each sampled vector of hazards
    S_pfs_pembro_wane <- matrix(nrow = nrow(S_pfs_pembro), ncol = ncol(S_pfs_pembro))
    for (i in 1:n_sim) {
      
      # apply waning function
      h_wane_pfs_pembro <- haz_pfs_pembro[,i] + punif((0:(n_cycles-2) - wane_start_times[i]) / (wane_end_times[i] - wane_start_times[i]), 0, 1) * pmax(0, diff_haz_mean_pfs)
      S_pfs_pembro_wane[,i] <- cumprod(c(1, exp(-h_wane_pfs_pembro)))
    }
    
    ### constrain PFS curves by OS curves
    S_pfs_pembro_wane <- replace(S_pfs_pembro_wane, S_pfs_pembro_wane > S_os_pembro_wane, S_pfs_pembro_wane[S_pfs_pembro_wane > S_os_pembro_wane])
    
    ### constrain ToT curves by OS curves
    S_tot_pembro_wane <- replace(S_tot_pembro, S_tot_pembro > S_os_pembro_wane, S_os_pembro_wane[S_tot_pembro > S_os_pembro_wane])
    S_tot_axi_wane <- replace(S_tot_axi, S_tot_axi > S_os_pembro_wane, S_os_pembro_wane[S_tot_axi > S_os_pembro_wane])
    
  }
  
  # Plot mean OS and PFS
  #plot_surv_fun(S_os_pembro_wane, S_os_suni, S_pfs_pembro_wane, S_pfs_suni)
  
  
  #####################################################################################
  # Model inputs
  #####################################################################################
  inputs <- list(
    
    # overall survival
    S_os_pembro = split(S_os_pembro, rep(1:ncol(S_os_pembro), each = nrow(S_os_pembro))),
    S_os_suni = split(S_os_suni, rep(1:ncol(S_os_suni), each = nrow(S_os_suni))),
    
    # progression-free survival
    S_pfs_pembro = split(S_pfs_pembro, rep(1:ncol(S_pfs_pembro), each = nrow(S_pfs_pembro))),
    S_pfs_suni = split(S_pfs_suni, rep(1:ncol(S_pfs_suni), each = nrow(S_pfs_suni))),
    
    # time on treatment
    S_tot_pembro =  split(S_tot_pembro, rep(1:ncol(S_tot_pembro), each = nrow(S_tot_pembro))),
    S_tot_axi = split(S_tot_axi, rep(1:ncol(S_tot_axi), each = nrow(S_tot_axi))),
    S_tot_suni = split(S_tot_suni, rep(1:ncol(S_tot_suni), each = nrow(S_tot_suni))),
    
    # utility values based on proximity to death
    ttd_week_index = rep(list(round(c(0, 30, 90, 180, 360)/7)), n_sim),  # thresholds in weeks for applying time to death utility weights
    
    u_ttd_0_29 = rbeta(n_sim, betaPar(0.462, 0.037)$a, betaPar(0.462, 0.037)$b), # Utility time to death <30 days
    u_ttd_30_89 = rbeta(n_sim, betaPar(0.594, 0.023)$a, betaPar(0.594, 0.023)$b),  # Utility time to death days [30,89)
    u_ttd_90_179 = rbeta(n_sim, betaPar(0.750, 0.019)$a, betaPar(0.750, 0.019)$b), # Utility time to death days [90,179) 
    u_ttd_180_359 = rbeta(n_sim, betaPar(0.769, 0.017)$a, betaPar(0.769, 0.017)$b), # Utility time to death days [180,359) 
    u_ttd_360_ = rbeta(n_sim, betaPar(0.824, 0.01)$a, betaPar(0.824, 0.01)$b), # Utility time to death >=360 days 
    
    du_ae = rnorm(n_sim, 0.057, 0.01), # AE-related disutility 
    
    # regimen related costs
    c_pembro = rep(5260, n_sim)*(1-rebate_pembro), # c_simple_chemo1 every 3 weeks
    c_axi = rep(3517, n_sim)*(1-rebate_axi), # orally twice daily continuously
    c_suni = rep(3138.8, n_sim)*(1-rebate_suni), # orally for 4 weeks then 2 weeks off treatment
    
    # first line dose intensity (assuming SE of 5% of the mean)
    dose1_pembro = rnorm(n_sim, 4986.48/5260, 4986.48/5260 * 0.05),
    dose1_axi = rnorm(n_sim, 2975.38/3517, 2975.38/3517 * 0.05),
    dose1_suni = rnorm(n_sim, 2344.68/3138.8, 2344.68/3138.8 * 0.05),
    
    # administration cost for IV
    c_simple_chemo1 = rgamma(n_sim, shape = gammaPar(174.40, Srange(156.96, 191.84)$s)$shape, scale = gammaPar(174.40, Srange(156.96, 191.84)$s)$scale),  # deliver simple parenteral chemotherapy at first attendance
    c_complex_chemo1 = rgamma(n_sim, shape = gammaPar(309.20, Srange((309.20/174.40)*156.96, (309.20/174.40)*191.84)$s)$shape, scale = gammaPar(309.20, Srange((309.20/174.40)*156.96, (309.20/174.40)*191.84)$s)$scale), # Deliver Complex Chemotherapy, including Prolonged Infusional Treatment, at First Attendance 
    c_oral_chemo = rep(0,n_sim),  
    # c_oral_chemo =rgamma(n_sim, shape = gammaPar(131.60, Srange(156.96, 191.84)$s)$shape, scale = gammaPar(174.40, Srange(156.96, 191.84)$s)$scale),  # company estimate
    
    # disease management costs
    c_week_pfs0 = rgamma(n_sim, shape = gammaPar(280.05, Srange(252.05, 308.06)$s)$shape, scale = gammaPar(280.05, Srange(252.05, 308.06)$s)$scale), # weekly cost in PFS state (cycle 0)
    c_week_pfs = rgamma(n_sim, shape = gammaPar(51.05, Srange(45.95, 56.16)$s)$shape, scale = gammaPar(51.05, Srange(45.95, 56.16)$s)$scale), # Weekly cost in progression-free state (subsequent cycles) 
    c_week_pps = rgamma(n_sim, shape = gammaPar(51.05, Srange(45.95, 56.16)$s)$shape, scale = gammaPar(51.05, Srange(45.95, 56.16)$s)$scale), # weekly cost in progressive disease state
    
    # subsequent treatment costs (assuming SE of 10% of the mean)
    # mean subsequent treatment cost for pembro = 10352.28 - (0.3*20884.69 + 0.2 * 18140.54)) + (0.20*63235.08) + (0.2*20884.69) + (0.1*18140.54) = 19096.77
    c_subs_tr_pembro =  rgamma(n_sim,  shape = gammaPar(19096.77, 19096.77*0.1)$shape, scale = gammaPar(19096.77, 19096.77*0.1)$scale), # Subsequent treatment cost (following intervention) # ERG/NICE estimate
    c_subs_tr_suni = rgamma(n_sim, shape = gammaPar(24700.62, 24700.62*0.1)$shape, scale = gammaPar(24700.62, 24700.62*0.1)$scale), # Subsequent treatment cost (following comparator) # ERG/NICE estimate
    # c_subs_tr_pembro =  rgamma(n_sim,  shape = gammaPar(9100.17, 9100.17*0.1)$shape, scale = gammaPar(9100.17, 9100.17*0.1)$scale), # Subsequent treatment cost (following intervention) # company estimate
    # c_subs_tr_suni = rgamma(n_sim, shape = gammaPar(20480.90, 20480.90*0.1)$shape, scale = gammaPar(20480.90, 20480.90*0.1)$scale), # Subsequent treatment cost (following comparator) # company estimate
    
    # end of life costs (assuming SE proportional to reported in company submission (page 204 committee papers))
    c_terminal = rgamma(n_sim, shape = gammaPar(8073, 8073 / 6789.76 * Srange(6110.78, 7468.74)$s)$shape, scale = gammaPar(8073, 8073 / 6789.76 * Srange(6110.78, 7468.74)$s)$scale),  # ERG/NICE estimate
    # c_terminal = rgamma(n_sim, shape = gammaPar(6789.76, Srange(6110.78, 7468.74)$s)$shape, scale = gammaPar(6789.76, Srange(6110.78, 7468.74)$s)$scale), # company estimate
    
    # adverse event management costs
    c_ae_pembro = rgamma(n_sim, shape = gammaPar(379.90, Srange(341.91, 417.89)$s)$shape, scale = gammaPar(379.90, Srange(341.91, 417.89)$s)$scale),
    c_ae_suni =  rgamma(n_sim, shape = gammaPar(348.34, Srange(313.51, 383.17)$s)$shape, scale = gammaPar(348.34, Srange(313.51, 383.17)$s)$scale),
    
    # administration schedules
    three_week_admin = rep(list(rep(c(1,0,0), n_cycles)[1:n_cycles]), n_sim),
    four_week_admin =  rep(list(rep(c(1,0,0,0), n_cycles)[1:n_cycles]), n_sim),
    six_week_admin =  rep(list(rep(c(1,0,0,0,0,0), n_cycles)[1:n_cycles]), n_sim)
    
  )
  
  
  if (waning == TRUE) {
    inputs$S_os_pembro <- split(S_os_pembro_wane, rep(1:ncol(S_os_pembro_wane), each = nrow(S_os_pembro_wane))) 
    inputs$S_pfs_pembro <- split(S_pfs_pembro_wane, rep(1:ncol(S_pfs_pembro_wane), each = nrow(S_pfs_pembro_wane))) 
  }
  
  
  
  #####################################################################################
  # Partitioned survival analysis model
  #####################################################################################
  partSA_fun <- function(inputs) {
    with(as.list(inputs),
         {
           
           #################################
           # Costs pembrolizumab-axitinib arm
           #################################
           
           # on treatment costs
           c_pembro_tr2 <- S_tot_pembro * three_week_admin * (c_simple_chemo1 + c_pembro*price * dose1_pembro)
           c_axi_tr2 <- S_tot_axi * four_week_admin * (c_oral_chemo + c_axi * dose1_axi)
           
           # PFS costs
           c_pfs_tr2 <- c(c_week_pfs0, S_pfs_pembro[2:n_cycles] * c_week_pfs)
           
           # PPS costs
           c_pps_tr2 <- (S_os_pembro - S_pfs_pembro) * c_week_pps
           
           # subsequent treatment costs
           c_subs_tr2 <- c(0, S_pfs_pembro[1:(length(S_pfs_pembro)-1)] - S_pfs_pembro[2:length(S_pfs_pembro)]) * c_subs_tr_pembro
           
           # adverse events costs
           c_ae_tr2 <- c(c_ae_pembro, rep(0, n_cycles-1))
           
           # terminal care costs
           c_terminal_tr2 <- c(0, S_os_pembro[1:(length(S_os_pembro)-1)] - S_os_pembro[2:length(S_os_pembro)]) * c_terminal
           
           # vector of total undiscounted costs
           v_tot_cost_tr2 <- c_pembro_tr2 + c_axi_tr2 + c_pfs_tr2 + c_pps_tr2 + c_subs_tr2 + c_ae_tr2 + c_terminal_tr2
           
           # sum of total discounted costs
           cost_tr2 <- v_tot_cost_tr2 %*% v_dr_c 
           
           #################################
           # Costs sunitinib
           #################################
           
           # on treatment costs
           c_suni_tr1 <- S_tot_suni * six_week_admin * (c_oral_chemo + c_suni * dose1_suni)
           
           # PFS costs
           c_pfs_tr1 <- c(c_week_pfs0, S_pfs_suni[2:n_cycles] * c_week_pfs)
           
           # PPS costs
           c_pps_tr1 <- (S_os_suni - S_pfs_suni) * c_week_pps
           
           # subsequent treatment costs
           c_subs_tr1 <- c(0, S_pfs_suni[1:(length(S_pfs_suni)-1)] - S_pfs_suni[2:length(S_pfs_suni)]) * c_subs_tr_suni
           
           # adverse events costs
           c_ae_tr1 <- c(c_ae_suni, rep(0, n_cycles-1))
           
           # terminal care costs
           c_terminal_tr1 <- c(0, S_os_suni[1:(length(S_os_suni)-1)] - S_os_suni[2:length(S_os_suni)]) * c_terminal
           
           # vector of total undiscounted costs
           v_tot_cost_tr1 <- c_suni_tr1 + c_pfs_tr1 + c_pps_tr1 + c_subs_tr1 +  c_ae_tr1 + c_terminal_tr2
           
           # sum of total discounted costs
           cost_tr1 <- v_tot_cost_tr1 %*% v_dr_c
           
           #################################
           # Life years pembrolizumab-axitinib
           #################################
           
           ly_tr2 <- v_dr_e %*% S_os_pembro  * (1/year2week)
           
           #################################
           # Life years sunitinib
           #################################
           
           ly_tr1 <- v_dr_e %*% S_os_suni  * (1/year2week)
           
           #################################
           # QALYs pembrolizumab-axitinib
           #################################
           
           # compute proportion of patients in TTD states
           ttd_0_29_pembro <- S_os_pembro[ttd_week_index[1]:n_cycles] - S_os_pembro[c(ttd_week_index[2]:n_cycles, rep(n_cycles, ttd_week_index[2] - ttd_week_index[1] - 1) )]
           
           ttd_30_89_pembro <- c(S_os_pembro[ttd_week_index[2]:n_cycles] - S_os_pembro[c(ttd_week_index[3]:n_cycles, rep(n_cycles, ttd_week_index[3] - ttd_week_index[2]) )],
                                 rep(0, ttd_week_index[2] - 1))
           
           ttd_90_179_pembro <- c(S_os_pembro[ttd_week_index[3]:n_cycles] - S_os_pembro[c(ttd_week_index[4]:n_cycles, rep(n_cycles, ttd_week_index[4] - ttd_week_index[3]) )],
                                  rep(0, ttd_week_index[3] - 1))
           
           ttd_180_359_pembro <- c(S_os_pembro[ttd_week_index[4]:n_cycles] - S_os_pembro[c(ttd_week_index[5]:n_cycles, rep(n_cycles, ttd_week_index[5] - ttd_week_index[4]) )],
                                   rep(0, ttd_week_index[4] - 1))
           
           ttd_360_pembro <- 1 - (ttd_0_29_pembro + ttd_30_89_pembro + ttd_90_179_pembro + ttd_180_359_pembro + (1-S_os_pembro) )
           ttd_360_pembro <- c(ttd_360_pembro[1:(n_cycles - ttd_week_index[5])], rep(0, ttd_week_index[5]))
           
           # compute time to death based utility
           u_ttd_0_29_pembro <- ttd_0_29_pembro * u_ttd_0_29
           u_ttd_30_89_pembro <-  ttd_30_89_pembro * u_ttd_30_89
           u_ttd_90_179_pembro <- ttd_90_179_pembro * u_ttd_90_179
           u_ttd_180_359_pembro <- ttd_180_359_pembro * u_ttd_180_359
           u_ttd_360_pembro <- ttd_360_pembro * u_ttd_360_
           
           # vector of total undiscounted QALYs
           v_tot_qaly_tr2 <- u_ttd_0_29_pembro + u_ttd_30_89_pembro + u_ttd_90_179_pembro + u_ttd_180_359_pembro + u_ttd_360_pembro #- (u_age_decr_fun(male, age) * S_os_pembro)
           
           # sum of total discounted QALYs
           qaly_tr2 <- v_tot_qaly_tr2 %*% v_dr_e * (1/year2week) - du_ae
           
           #################################
           # QALYs sunitinib
           #################################
           
           # compute proportion of patients in TTD states
           ttd_0_29_suni <- S_os_suni[ttd_week_index[1]:n_cycles] - S_os_suni[c(ttd_week_index[2]:n_cycles, rep(n_cycles, ttd_week_index[2] - ttd_week_index[1] - 1) )]
           
           ttd_30_89_suni <- c(S_os_suni[ttd_week_index[2]:n_cycles] - S_os_suni[c(ttd_week_index[3]:n_cycles, rep(n_cycles, ttd_week_index[3] - ttd_week_index[2]) )],
                               rep(0, ttd_week_index[2] - 1))
           
           ttd_90_179_suni <- c(S_os_suni[ttd_week_index[3]:n_cycles] - S_os_suni[c(ttd_week_index[4]:n_cycles, rep(n_cycles, ttd_week_index[4] - ttd_week_index[3]) )],
                                rep(0, ttd_week_index[3] - 1))
           
           ttd_180_359_suni <- c(S_os_suni[ttd_week_index[4]:n_cycles] - S_os_suni[c(ttd_week_index[5]:n_cycles, rep(n_cycles, ttd_week_index[5] - ttd_week_index[4]) )],
                                 rep(0, ttd_week_index[4] - 1))
           
           ttd_360_suni <- 1 - (ttd_0_29_suni + ttd_30_89_suni + ttd_90_179_suni + ttd_180_359_suni + (1-S_os_suni) )
           ttd_360_suni <- c(ttd_360_suni[1:(n_cycles - ttd_week_index[5])], rep(0, ttd_week_index[5]))
           
           # compute time to death based utility
           u_ttd_0_29_suni <- ttd_0_29_suni * u_ttd_0_29
           u_ttd_30_89_suni <-  ttd_30_89_suni * u_ttd_30_89
           u_ttd_90_179_suni <- ttd_90_179_suni * u_ttd_90_179
           u_ttd_180_359_suni <- ttd_180_359_suni * u_ttd_180_359
           u_ttd_360_suni <- ttd_360_suni * u_ttd_360_
           
           # vector of total undiscounted QALYs
           v_tot_qaly_tr1 <- u_ttd_0_29_suni + u_ttd_30_89_suni + u_ttd_90_179_suni + u_ttd_180_359_suni + u_ttd_360_suni # - (u_age_decr_fun(male, age) * S_os_suni)
           
           # sum of total discounted QALYs
           qaly_tr1 <- v_tot_qaly_tr1 %*% v_dr_e * (1/year2week) - du_ae
           
           #################################
           # store results in dataframe
           #################################
           results <- data.frame(v_cost1 = cost_tr1, v_cost2 = cost_tr2,
                                 v_ly1 = ly_tr1, v_ly2 = ly_tr2,
                                 v_qaly1 = qaly_tr1, v_qaly2 = qaly_tr2)
           
           return(results)
           
           
         })
  }
  
  # convert IPD from weekly times to months
  ipd_os$tt <- ipd_os$tt / year2week * month2year
  ipd_os$treat <- ifelse(ipd_os$treat == 1, 2, 1)
  ipd_pfs$tt <- ipd_pfs$tt / year2week * month2year
  ipd_os$treat <- ifelse(ipd_pfs$treat == 1, 2, 1)
  
  #####################################################################################
  # run the PartSA model for each parameter sample
  #####################################################################################
  results <- lapply(1:n_sim, function (x) {
    partSA_fun(lapply(inputs, `[[`, x))
  })
  results <- do.call(rbind, results)
  
  return(list(v_cost1 = results$v_cost1, v_cost2 = results$v_cost2, v_ly1 = results$v_ly1, v_ly2 = results$v_ly2, v_qaly1 = results$v_qaly1, v_qaly2 = results$v_qaly2,
              l_os = list(m_os1 = S_os_suni, m_os2 = S_os_pembro_wane), l_pfs = list(m_pfs1 = S_pfs_suni, m_pfs2 = S_pfs_pembro_wane),
              m_ipd_os = ipd_os, m_ipd_pfs = ipd_pfs))
  
} # end pa_mod_fun




