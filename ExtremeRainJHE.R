####################################
### Jason West, Bureau of Meteorology
### jason.west@bom.gov.au
### EXTREME RAINFALL INFERENCE ANALYSIS
####################################

library(data.table)
library(lubridate)
library(MASS)
library(ggplot2)
library(extRemes)
library(xts)
library(tidyquant)
library(Kendall)

# GET RAIN DATA
# load full file
prec<-read.csv("Master_rain_daily.csv") # data available here: https://zenodo.org/records/19341730
prec<-as.data.table(prec)
prec$X<-NULL
colnames(prec)<-c("Date","Stnum","Precip")
prec$Date<-ymd(prec$Date)
prec<-prec[,Precip:=ifelse(Precip==99999.9,0,Precip)]
stations<-unique(prec$Stnum)

# Stations for linking to data
Stnums.aust<-read.csv("Stn_nums_all.csv") # data available at this repository
Stnums.aust<-as.data.table(Stnums.aust)
Stnums.aust$Stnum<-as.integer(Stnums.aust$Site)
Stnums.aust$Fin<-as.Date(paste("01-", Stnums.aust$End, sep = ""), format = "%d-%b-%y")
Stnums.aust$Commence<-as.Date(paste("01-", Stnums.aust$Start, sep = ""), format = "%d-%b-%y")

## FORM DATA SET WITH STNUMS, NAMES, PRECIP, AND LAT/LONG
precnums<-merge.data.table(prec,Stnums.aust,by="Stnum")
precnums$Site<-NULL
precnums$X.<-NULL
summary(precnums)
precnums$Date<-ymd(precnums$Date)
precnums$Lat<-as.numeric(precnums$Lat)
precnums$Lon<-as.numeric(precnums$Lon)
precnums$Years<-as.numeric(precnums$Years)
precnums$AWS<-as.factor(precnums$AWS)
str(precnums)

# Explore data
maxprecip<-precnums[,max(Precip),by=.(Stnum)]

stations.plot<-merge.data.table(maxprecip,Stnums.aust,by="Stnum")
stations.plot$Commence<-stations.plot$Fin-stations.plot$Years*365
st.plot<-stations.plot[,.(Stnum,Commence,Fin)]
library(reshape2)
st.plot.long<-melt(st.plot,id.vars = "Stnum")
detach("package:reshape2", unload=TRUE)
st.plot.long$value <- as.POSIXct(st.plot.long$value)
st.plot.long<-as.data.table(st.plot.long)
st.plot.long<-st.plot.long[order(st.plot.long$Stnum,st.plot.long$value)]
st.plot.long$StnumID <- rep(seq(1, nrow(st.plot.long)/2, 1), each=2)
str(st.plot.long)

p1 <- ggplot(data=st.plot.long) +
  geom_line(aes(x=value, y=Stnum, group=StnumID)) +
  xlab("Year") + ylab("Station number") +
  labs(title = "Rainfall record, Australia 1850-2023",
       subtitle = "By Station (continuous recording period)") +
  theme_bw() + theme(legend.position="none")

p1

####################################
# Fig 2
maxprecip<-maxprecip[order(maxprecip$V1)]
str(maxprecip)
p2 <- ggplot(data=maxprecip, aes(x=Stnum,y=V1)) +
  geom_point(pch=4) + coord_flip() + theme_bw() +
  xlab("Station number") + ylab("Maximum daily rainfall (mm)") +
  labs(title = "Maximum rainfall, Australia 1850-2023",
       subtitle = "Daily Rainfall Max, By Station") +
  theme_bw() + theme(legend.position="none")

p2
library(gridExtra)
grid.arrange(p1,p2,nrow = 1) # gridextra package

# convert prec to xts object
library(tidyquant)
rain<-precnums[,c(2,1,3)]
rain.xts<-as.xts(rain[Stnum==5008,],date_col=Date)

# Derive Annual Maxima (Minima) Series (AMS) for maximum precipitation
ams <- apply.yearly(rain.xts, max)
ams <- apply.yearly(prec[Stnum==5008,-2], max)
hist(ams,breaks = 80)
# maximum-likelihood fitting of the GEV distribution
fit_mle <- fevd(as.vector(ams), method = "MLE", type="GEV") # method can be "MLE", "GMLE", "Bayesian", "Lmoments"
# diagnostic plots
plot(fit_mle)
# return levels:
rl_mle <- return.level(fit_mle, conf = 0.05, return.period= c(2,5,10,20,50,100))
rl_mle
fit_mle$results$par # EV params

# fitting GEV distribution based on L-moments estimation
fit_lmom <- fevd(as.vector(ams), method = "Lmoments", type="GEV")
# diagnostic plots
plot(fit_lmom)
# return levels:
rl_lmom <- return.level(fit_lmom, conf = 0.05, return.period= c(2,5,10,20,50,100))
rl_lmom
fit_lmom$results

# return level plots
par(mfcol=c(1,2))
# return level plot w/ MLE
plot(fit_mle, type="rl",
     main="Return Level Plot (MLE)",
     ylim=c(0,400), pch=16)
loc <- as.numeric(return.level(fit_mle, conf = 0.05,return.period=100))
segments(100, 0, 100, loc, col= 'midnightblue',lty=6)
segments(0.01,loc,100, loc, col='midnightblue', lty=6)

# return level plot w/ LMOM
plot(fit_lmom, type="rl",
     main="Return Level Plot (L-Moments)",
     ylim=c(0,400), pch=16)
loc <- as.numeric(return.level(fit_lmom, conf = 0.05,return.period=100))
segments(100, 0, 100, loc, col= 'midnightblue',lty=6)
segments(0.01,loc,100, loc, col='midnightblue', lty=6)

# comparison of return levels
results <- t(data.frame(mle=as.numeric(rl_mle),
                        lmom=as.numeric(rl_lmom)))
colnames(results) <- c(2,5,10,20,50,100)
round(results,1)

# In most cases, L-moments estimation is more robust than maximum likelihood estimation
par(mfcol=c(1,1))

##########################################################################

##########################################################################
##########################################################################
# simple return level plot (base graphics) from fevd objects that have been fitted of type GEV.
library(dplyr)
rlplot_gev <- function(fevd_obj,
                       pointcolor = "firebrick",
                       linecolor = "darkgreen",
                       yticks = seq(0,225,25)){
  rperiods = c(1+1e-15, 1.25, 1.5, 2, 3,4,5, 7.5, 10, 15, 20,
               30, 40, 50, 60, 75, 80, 100, 120, 200, 250)
  model <- fevd_obj$type
  if(model != "GEV"){
    stop("this plotting function only works for GEV")
  }
  # pars <- ifelse(fevd_obj$method=="MLE",
  #                fevd_obj$results$par,fevd_obj$results)
  # pars<-fevd_obj$results$par
  params<-function(fevd_obj){
    if(fevd_obj$method=="MLE"){p<-fevd_obj$results$par}
    else{p<-fevd_obj$results}
    return(p)
  }
  pars<-params(fevd_obj)
  loc <- pars["location"]
  scale <- pars["scale"]
  shape <- pars["shape"]
  xrl <- -1/(log(1 - 1/rperiods))
  yrl <- rlevd(rperiods, loc = loc, scale = scale, 
               shape = shape)
  mod<-ifelse(fevd_obj$method=="MLE",
              as.character("MLE"),
              as.character("L-Moments"))
  # empty plot
  plot(xrl, yrl, type = "n", log = "x", xlab = "Return Period (years)",
       ylab = "Return Level (mm)", lwd = 2, yaxt = "n", xlim = c(1,220),
       main = paste("Return Level Plot",mod))
  par(las = 1)
  axis(side = 2, at = yticks)
  
  # background grid
  abline(h = yticks,col = "lightgrey", lty = 2)
  abline(v = c(1,2,5,10,20,50,100,200), col = "lightgrey", lty = 2)
  
  # plot line
  lines(xrl, yrl, col = linecolor, lwd = 2)
  
  # ci lines for MLE base plot
  bds <- ci(fevd_obj, return.period = rperiods)
  lines(xrl, bds[, 1], col = linecolor, lty = 4)
  lines(xrl, bds[, 3], col = linecolor, lty = 4)
  
  # Weibull plotting position of points  
  n <- fevd_obj$n
  xp <- ppoints(n = n, a = 0)
  ytmp <- datagrabber(fevd_obj)
  y <- c(ytmp[, 1])
  points(-1/log(xp), sort(y), pch = 16, col = pointcolor)
}

# plot return levels for mle and lmoments
par(mfcol=c(1,2))
rlplot_gev(fit_mle)
rlplot_gev(fit_lmom)
par(mfcol=c(1,1))

#######################################
################
#### Revised plots for paper
# Core packages
library(tidyverse)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(akima)     # for interpolation (thin-plate-like via IDW)
library(patchwork) # for multi-panel figures
library(data.table)

# Read the data
rl <- read_csv("RL_results_extreme_rain.csv")

# Check
glimpse(rl)

# Derive method differences (explicit and transparent)
rl <- rl %>%
  mutate(
    d_bayes_gev = RL_Bayes - RL_GEV,
    d_bayes_pot = RL_Bayes - RL_POT,
    d_bayes_mev = RL_Bayes - RL_MEV,
    d_pot_gev   = RL_POT   - RL_GEV
  )

IQR(rl$d_bayes_gev)
IQR(rl$d_bayes_pot)
IQR(rl$d_bayes_mev)
IQR(rl$d_pot_gev)
perc_above <- sum(rl$RL_Bayes>rl$RL_GEV)/nrow(rl)*100
perc_above
perc_above <- sum(rl$RL_Bayes>rl$RL_POT)/nrow(rl)*100
perc_above
median(rl$d_bayes_pot)
median(rl$d_bayes_gev)

#  Base map of Australia (vector, clean, reproducible)
aus <- ne_countries(scale = "medium",
                    country = "Australia",
                    returnclass = "sf")

# Figure 2×2 panel map (replacement)
plot_diff_map <- function(data, zvar, title,
                          zlim = c(-150, 300)) {
  
  ggplot() +
    geom_sf(data = aus, fill = "grey95", colour = "grey40") +
    geom_point(
      data = data,
      aes(x = Longitude, y = Latitude, colour = !!sym(zvar)),
      size = 2.4
    ) +
    scale_colour_gradient2(
      low = "#d7191c",
      mid = "white",
      high = "#2c7bb6",
      midpoint = 0,
      limits = zlim,
      name = "Difference (mm)"
    ) +
    coord_sf(xlim = c(110, 155), ylim = c(-45, -10)) +
    labs(title = title) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )
}

# Structure of plot
p1 <- plot_diff_map(rl, "d_bayes_gev",
                    "Bayesian − GEV (Block maxima)")

p2 <- plot_diff_map(rl, "d_bayes_pot",
                    "Bayesian − POT")

p3 <- plot_diff_map(rl, "d_bayes_mev",
                    "Bayesian − MEV")

p4 <- plot_diff_map(rl, "d_pot_gev",
                    "POT − GEV")

# Combine
figure_X <- (p1 | p2) / (p3 | p4)

figure_X

##########################################################
### RECONFIGURE CODE FOR AMENDED PAPER
##########################################################

### PART 1

# This code implements:
#  - model x inference separation
#  - uncertainty quantification
#  - cross-validation
#  - POT correctness
#  - reproducible multi-station workflow
#  - Bayesian strengthening via elicited quantiles

# 1. sETUP
library(data.table)
library(lubridate)
library(extRemes)
library(evir)
library(ggplot2)
library(dplyr)
library(purrr)
library(parallel)

setDTthreads(threads = detectCores())

# 2. DATA PREP FUNCTION
load_data <- function(path_prec, path_meta){
  prec <- fread(path_prec)
  setnames(prec, c("Date","Stnum","Precip"))
  prec[, Date := ymd(Date)]
  prec[, Precip := ifelse(Precip == 99999.9, NA, Precip)]
  
  meta <- fread(path_meta)
  meta[, Stnum := as.integer(Site)]
  
  data <- merge(prec, meta, by="Stnum")
  data
}

# 3. ANNUAL MAXIMA (WITH WATER YEAR)
compute_ams <- function(data, st, water_year=TRUE){
  dt <- data[Stnum == st & !is.na(Precip)]
  
  if(water_year){
    dt[, year := year(Date + months(6))]
  } else {
    dt[, year := year(Date)]
  }
  
  ams <- dt[, .(max_precip = max(Precip)), by=year]
  return(ams$max_precip)
}

# 4. GEV FITTING (MLE + L-MOMENTS)
fit_gev_models <- function(ams){
  fit_mle  <- fevd(ams, type="GEV", method="MLE")
  fit_lmom <- fevd(ams, type="GEV", method="Lmoments")
  
  list(mle=fit_mle, lmom=fit_lmom)
}

# 5. RETURN LEVEL + UNCERTAINTY
extract_rl <- function(fit, rp=100){
  rl <- return.level(fit, return.period=rp)
  ci_bounds <- ci(fit, return.period=rp)
  
  data.frame(
    RL = as.numeric(rl),
    lower = ci_bounds[1],
    upper = ci_bounds[3],
    width = ci_bounds[3] - ci_bounds[1]
  )
}

# 6. CROSS-VALIDATION (TAIL-FOCUSED)
cv_tail <- function(ams){
  n <- length(ams)
  split <- floor(0.8*n)
  
  train <- ams[1:split]
  test  <- ams[(split+1):n]
  
  fit_mle <- fevd(train, type="GEV", method="MLE")
  
  pars <- fit_mle$results$par
  
  threshold <- quantile(train, 0.95)
  
  pred_exceed <- 1 - pevd(threshold,
                          loc=pars["location"],
                          scale=pars["scale"],
                          shape=pars["shape"])
  
  obs_exceed <- mean(test > threshold)
  
  data.frame(
    pred = pred_exceed,
    obs = obs_exceed,
    bias = obs_exceed - pred_exceed
  )
}

# 7. POT WITH THRESHOLD + DECLUSTERING
fit_pot <- function(data, st, u=50){
  dt <- data[Stnum==st & Precip > u]
  
  # decluster
  dc <- decluster(dt$Precip, u=u, method="runs", r=3)
  
  fit <- fevd(dc - u, type="GP", method="MLE")
  
  list(fit=fit, n_clusters=length(dc))
}

# 8. BAYESIAN FIT (DEFAULT extRemes)
fit_bayes_extRemes <- function(ams){
  
  fit <- fevd(
    ams,
    type = "GEV",
    method = "Bayesian"
  )
  
  posterior <- as.data.frame(fit$results)
  #posterior <- posterior[1001:nrow(posterior), ]
  
  # Convert log.scale → scale
  posterior$scale <- exp(posterior$log.scale)
  
  # Compute return levels (100-year)
  rlev_post <- with(posterior, {
    ifelse(abs(shape) < 1e-6,
           location - scale * log(-log(1 - 0.01)),  # Gumbel case
           location + (scale / shape) *
             ((-log(1 - 0.01))^(-shape) - 1)
    )
  })
  
  data.frame(
    RL = median(rlev_post, na.rm=TRUE),
    lower = quantile(rlev_post, 0.025, na.rm=TRUE),
    upper = quantile(rlev_post, 0.975, na.rm=TRUE),
    width = quantile(rlev_post, 0.975, na.rm=TRUE) -
      quantile(rlev_post, 0.025, na.rm=TRUE)
  )
}

# 9. SINGLE-STATION PIPELINE
process_station <- function(data, st){
  
  ams <- compute_ams(data, st)
  
  if(length(ams) < 30) return(NULL)
  
  gev <- fit_gev_models(ams)
  
  mle_stats  <- extract_rl(gev$mle)
  lmom_stats <- extract_rl(gev$lmom)
  bayes_stats <- fit_bayes_extRemes(ams)
  
  cv <- cv_tail(ams)
  
  data.frame(
    Stnum = st,
    N = length(ams),
    
    RL_MLE = mle_stats$RL,
    W_MLE  = mle_stats$width,
    
    RL_LMOM = lmom_stats$RL,
    W_LMOM  = lmom_stats$width,
    
    RL_BAYES = bayes_stats$RL,
    W_BAYES  = bayes_stats$width,
    
    shrinkage = bayes_stats$width / mle_stats$width,
    
    CV_bias = cv$bias
  )
}

# 10. MULTI-STATION EXECUTION
run_pipeline <- function(data){
  
  stations <- unique(data$Stnum)
  
  res <- lapply(stations, function(s){
    tryCatch(process_station(data, s), error=function(e) NULL)
  })
  
  rbindlist(res, fill=TRUE)
}

# 11. NEW FIGURES
plot_results <- function(df){
  
  # Scatter
  p1 <- ggplot(df, aes(RL_MLE, RL_BAYES)) +
    geom_point() +
    geom_abline() +
    theme_bw()
  
  # Uncertainty
  p2 <- ggplot(df, aes(W_MLE, W_BAYES)) +
    geom_point() +
    geom_abline() +
    theme_bw()
  
  # Histogram
  p3 <- ggplot(df, aes(RL_BAYES - RL_MLE)) +
    geom_histogram(bins=40) +
    theme_bw()
  
  # Shrinkage
  p4 <- ggplot(df, aes(N, shrinkage)) +
    geom_point() +
    geom_smooth() +
    theme_bw()
  
  list(p1,p2,p3,p4)
}

#########################
# Run code:
# Structure is:
# run_pipeline()
# └── process_station()
# ├── compute_ams()
# ├── fit_gev_models()
# ├── extract_rl()
# ├── fit_bayes_extRemes()
# └── cv_tail()

data <- load_data(
  path_prec = "Master_rain_daily_raw.csv",
  path_meta = "Stn_nums_all.csv"
)

# ------ Run the Pipeline -------#
results <- run_pipeline(data)

# ------ Baseline Results -------#
summary(results$shrinkage)
median(results$shrinkage)

# ------ Relative Results -------#
results$reduction_pct <- 100 * (1 - results$shrinkage_expert)
summary(results$reduction_pct)

abs(results$RL_BAYES_EXP - results$RL_MLE)


# check
mean(results$W_BAYES > results$W_MLE)

# Initial plot to check results
ggplot(results, aes(W_MLE, W_BAYES)) +
  geom_point() +
  geom_abline() +
  labs(
    x = "MLE interval width",
    y = "Bayesian interval width"
  ) +
  theme_bw()

# Full Plots
plots <- plot_results(results)
plots[[1]]  # scatter
plots[[2]]  # uncertainty
plots[[3]]  # histogram
plots[[4]]  # shrinkage

# Run test station
# test <- process_station(data, st = 5008)
# print(test)

################################################################
# PART 2 — BAYES MODEL WITH EXPERT-ELICITED QUANTILES
# Experts don’t think in (mu, sigma, epsilon)
# They think in:
# - 10-year rainfall (z10)
# - 100-year rainfall (z100)
# put priors on these and transform

# DEFINE GEV LOG-LIKELIHOOD
gev_loglik <- function(y, mu, sigma, xi){
  
  if(sigma <= 0) return(-Inf)
  
  t <- 1 + xi * (y - mu) / sigma
  
  if(any(t <= 0)) return(-Inf)
  
  sum(-log(sigma) - (1 + 1/xi)*log(t) - t^(-1/xi))
}

# DEFINE QUANTILE PRIORS
gev_quantile <- function(mu, sigma, xi, p){
  
  if(abs(xi) < 1e-6){
    return(mu - sigma * log(-log(1 - p)))
  }
  
  mu + (sigma/xi) * ((-log(1 - p))^(-xi) - 1)
}

# Prior (based on z10 and z100)
log_prior <- function(mu, sigma, xi, z10, z100){
  
  z10_model  <- gev_quantile(mu, sigma, xi, 0.1)
  z100_model <- gev_quantile(mu, sigma, xi, 0.01)
  
  lp <- dnorm(z10_model,  mean=z10$mean,  sd=z10$sd,  log=TRUE) +
    dnorm(z100_model, mean=z100$mean, sd=z100$sd, log=TRUE) +
    dnorm(xi, 0, 0.2, log=TRUE)
  
  lp
}

# METROPOLIS–HASTINGS SAMPLER
run_mcmc_gev <- function(ams, z10, z100, n_iter=5000){
  
  # initialise
  mu    <- mean(ams)
  sigma <- sd(ams)
  xi    <- 0.1
  
  chain <- matrix(NA, n_iter, 3)
  colnames(chain) <- c("mu", "sigma", "xi")
  
  current_lp <- gev_loglik(ams, mu, sigma, xi) +
    log_prior(mu, sigma, xi, z10, z100)
  
  for(i in 1:n_iter){
    
    # propose
    mu_new    <- rnorm(1, mu,    5)
    sigma_new <- rlnorm(1, log(sigma), 0.1)
    xi_new    <- rnorm(1, xi, 0.05)
    
    proposed_lp <- gev_loglik(ams, mu_new, sigma_new, xi_new) +
      log_prior(mu_new, sigma_new, xi_new, z10, z100)
    
    # accept/reject
    if(log(runif(1)) < (proposed_lp - current_lp)){
      mu <- mu_new
      sigma <- sigma_new
      xi <- xi_new
      current_lp <- proposed_lp
    }
    
    chain[i,] <- c(mu, sigma, xi)
  }
  
  as.data.frame(chain)
}

# RETURN LEVEL EXTRACTION
extract_mcmc_rl <- function(chain){
  
  r100 <- with(chain, {
    ifelse(abs(xi) < 1e-6,
           mu - sigma * log(-log(1 - 0.01)),
           mu + (sigma/xi) * ((-log(1 - 0.01))^(-xi) - 1)
    )
  })
  
  data.frame(
    RL = median(r100),
    lower = quantile(r100, 0.025),
    upper = quantile(r100, 0.975),
    width = quantile(r100, 0.975) -
      quantile(r100, 0.025)
  )
}

# Process Station Function
process_station <- function(data, st){
  
  ams <- compute_ams(data, st)
  
  if(length(ams) < 30) return(NULL)
  
  # --- Frequentist ---
  gev <- fit_gev_models(ams)
  mle_stats  <- extract_rl(gev$mle)
  lmom_stats <- extract_rl(gev$lmom)
  
  # --- Bayesian (extRemes baseline) ---
  bayes_stats <- fit_bayes_extRemes(ams)
  
  # --- Bayesian (EXPERT PRIORS via custom MCMC) ---
  q10  <- quantile(ams, 0.9)
  q100 <- quantile(ams, 0.99)
  
  priors <- list(
    z10 = list(mean = q10,
               sd   = 0.15 * q10),
    
    z100 = list(mean = q100,
                sd   = 0.20 * q100)
  )
  
  chain <- tryCatch(
    run_mcmc_gev(ams, priors$z10, priors$z100, n_iter=5000),
    error = function(e) NULL
  )
  
    if(!is.null(chain) && nrow(chain) > 1000){
      chain <- chain[1001:nrow(chain), ]
      stan_stats <- extract_mcmc_rl(chain)
    } else {
      stan_stats <- data.frame(RL=NA, lower=NA, upper=NA, width=NA)
    }
  
  # --- Cross-validation ---
  cv <- cv_tail(ams)
  
  data.frame(
    Stnum = st,
    N = length(ams),
    
    # Frequentist
    RL_MLE = mle_stats$RL,
    W_MLE  = mle_stats$width,
    
    RL_LMOM = lmom_stats$RL,
    W_LMOM  = lmom_stats$width,
    
    # Bayesian (weak prior)
    RL_BAYES = bayes_stats$RL,
    W_BAYES  = bayes_stats$width,
    
    # Bayesian (expert prior)
    RL_BAYES_EXP = stan_stats$RL,
    W_BAYES_EXP  = stan_stats$width,
    
    shrinkage_weak   = bayes_stats$width / mle_stats$width,
    shrinkage_expert = stan_stats$width / mle_stats$width,
    
    CV_bias = cv$bias
  )
}

# process_station(data, 5008)
# print(test)

# Run full pipeline:
results <- run_pipeline(data)
summary(results$shrinkage_expert)

################################
# plots
library(ggplot2)

# Figure 3
ggplot(results, aes(x=W_MLE, y=W_BAYES)) +
  geom_point(alpha=0.7) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(
    x = "MLE interval width (mm)",
    y = "Bayesian interval width (weak prior, mm)",
    title = "Uncertainty comparison: MLE vs Bayesian (weak priors)"
  ) +
  theme_bw()

# Figure 4
ggplot(results, aes(x=W_MLE, y=W_BAYES_EXP)) +
  geom_point(alpha=0.7) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(
    x = "MLE interval width (mm)",
    y = "Bayesian interval width (expert priors, mm)",
    title = "Uncertainty reduction using expert-informed priors"
  ) +
  theme_bw()

# Figure 5
library(tidyr)

df_plot <- results %>%
  select(shrinkage_weak, shrinkage_expert) %>%
  pivot_longer(cols=everything(), names_to="type", values_to="value")

ggplot(df_plot, aes(x=value, fill=type)) +
  geom_histogram(alpha=0.6, bins=30, position="identity") +
  labs(
    x = "Shrinkage (relative uncertainty)",
    y = "Frequency",
    fill = "Method"
  ) +
  theme_bw()

# Figure 6
library(sf)
library(rnaturalearth)

# Stations for linking to data
Stnums.aust<-read.csv("Stn_nums_all.csv")
Stnums.aust<-as.data.table(Stnums.aust)
Stnums.aust$Stnum<-as.integer(Stnums.aust$Site)
Stnums.aust$Fin<-as.Date(paste("01-", Stnums.aust$End, sep = ""), format = "%d-%b-%y")
Stnums.aust$Commence<-as.Date(paste("01-", Stnums.aust$Start, sep = ""), format = "%d-%b-%y")
Stnums.loc <- Stnums.aust[,c("Stnum","Lat","Lon")]

aus <- ne_countries(country="Australia", returnclass="sf")
results_upd <- merge(results,Stnums.loc,by="Stnum")
results_sf <- st_as_sf(results_upd,
                       coords = c("Lon", "Lat"),
                       crs = 4326)

ggplot() +
  geom_sf(data=aus, fill="grey95") +
  geom_sf(data=results_sf,
          aes(color = (RL_BAYES_EXP - RL_MLE)/RL_MLE),
          size=2) +
  scale_color_gradient2(
    low="red", mid="white", high="blue", midpoint=0,
    name="Relative difference"
  ) +
  theme_minimal() +
  labs(title="Spatial variation in Bayesian vs MLE return levels")

# Figure 7
ggplot(results, aes(x=CV_bias)) +
  geom_histogram(bins=30, fill="steelblue", alpha=0.7) +
  labs(
    x = "Test minus predicted exceedance probability",
    y = "Frequency",
    title = "Cross-validation bias distribution"
  ) +
  theme_bw()

## Additional Plot
library(dplyr)
library(ggplot2)


results_binned <- results %>%
  mutate(N_bin = cut(N, breaks=c(80,90,100,110,120,130,140))) %>%
  group_by(N_bin) %>%
  summarise(
    N_mid = mean(N),
    shrink_mean = mean(shrinkage_expert),
    shrink_se = sd(shrinkage_expert)/sqrt(n()),
    n = n()
  )

ggplot(results_binned, aes(N_mid, shrink_mean)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin = shrink_mean - shrink_se,
                    ymax = shrink_mean + shrink_se),
                width=2) +
  geom_text(aes(label=n), vjust=-1, hjust=-1, size=3) +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(
    title = "Dependence of Bayesian uncertainty reduction on record length",
    x = "Record length (years)",
    y = "Relative uncertainty (Bayesian / MLE)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = 135, y = 0.9,
           label = "Likelihood dominates", size = 4.0) +
  annotate("text", x = 105, y = 0.6,
         label = "Prior influence strongest", size = 4.0)

# Tables
Stnums.nm <- Stnums.aust[,c("Stnum","Name")]
results_table <- merge(results,Stnums.nm,by="Stnum")

library(dplyr)
# Move column "Name" to the first position
results_table <- results_table %>%
  relocate(Name)

results_table %>%
  arrange(shrinkage_expert) %>%
  select(Name, shrinkage_expert) %>%
  head(10)

################################
## End
