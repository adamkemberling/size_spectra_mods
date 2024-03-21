####  Size Spectra Linear Model Selection  
# Taking this out of manuscript file to dig into outputs


####  Packages  ####
#library(EnvCpt)
{
  library(targets)
  library(here)
  library(gmRi)
  library(patchwork)
  library(scales)
  library(emmeans)
  library(tidyverse)
  library(ggpubr)
  library(DHARMa)
  library(nlme)
  library(ggeffects)
  
}

# Package Conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")




#### Load from `spectra_lm_data_prep.R`
model_data <- read_csv(here::here("data/size_and_spectra_model_data.R"))


# Function to grab the correlation data and lag data
get_ccf_vector <- function(x,y, lagmax = 10){
  
  # Run the ccf
  ccf_data <- ccf(x, y, plot= F , lag.max = lagmax)
  
  # Get the signif:
  # https://www.squaregoldfish.co.uk/programming/r_acf_significance.md/
  # Not 100% sure n is the same for ccf as it is for acf, but...
  significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(x)))
  
  data.frame(
    "acf"    = ccf_data$acf,
    "lag"    = ccf_data$lag,
    "sigpos" = significance_level,
    "signeg" = significance_level*-1
  )
}



####___________________#####



# Median Length ~ stuff
# Length is measured for all individuals, weight is estimated

# General/base model
# Median Body Size ~ region * (year + drivers) + GSI

#### Body-Size Models  ####


# Data for the size changes

# Drop NA's relevant to the sst data we want ot use
mod_dat_oisst <- model_data %>% drop_na(avg_len, med_len, sst, landings) %>% 
  arrange(year, survey_area) 

# Drop NA's with ERSST
mod_dat_ersst <- model_data %>% drop_na(avg_len, med_len, ersst_anom, landings, zp_large) %>% 
  arrange(year, survey_area) 

#### a.) Base Model  ####

# Make the global model

# OISST version
length_lm <- lm(
  med_len ~ survey_area * (year + zp_large + landings + sst) + gsi_oisst_resid, 
  data = mod_dat_oisst)

# # ERSST version
# length_lm <- lm(
#   med_len ~ survey_area * (year + zp_large + landings + ersst_anom) + gsi_ersst_resid, 
#   data = mod_dat_ersst)




#####  Model Summary  ####
summary(length_lm)

# Check normality of residuals - iffyy
hist(resid(length_lm))


# Model evaluation library
res <- simulateResiduals(length_lm, plot = T)

#####  VIF  ####
# Check variance inflation, hard to appreciate with intentional interactions
# time * area

# Resource https://stacyderuiter.github.io/s245-notes-bookdown/collinearity-and-multicollinearity.html
# If VIF (or squared scaled GVIF) is greater than 4, then there’s a problem and you should probably try to fix it;
# if VIF (or squared scaled GVIF) is more than 10, then something definitely must be done to correct the problem.
mod_vif <- car::vif(length_lm, type = "predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2) %>% 
  rownames_to_column("predictor")
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)



# Sequential removal:
# From Zuur 2010
# collinearity has been removed by sequentially deleting each 
# variable for which the VIF value was highest until all remaining VIFs were below 3.

# Drop landings
vif_drop_landings <-  lm(
  med_len ~ survey_area * (year + zp_large + sst) + gsi_oisst_resid, 
  data = mod_dat_oisst)
mod_vif <- car::vif(vif_drop_landings, type = "predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2) %>% 
  rownames_to_column("predictor")
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# Drop SSt
vif_drop_sst <-  lm(
  med_len ~ survey_area * (year + zp_large) + gsi_oisst_resid, 
  data = mod_dat_oisst)
mod_vif <- car::vif(vif_drop_sst, type = "predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2) %>% 
  rownames_to_column("predictor")
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# Drop zooplankton
vif_drop_zoo <-  lm(
  med_len ~ survey_area * (year) + gsi_oisst_resid, 
  data = mod_dat_oisst)
mod_vif <- car::vif(vif_drop_zoo, type = "predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2) %>% 
  rownames_to_column("predictor")
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 3)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)



# swap model over so we can keep the code below consistent
length_lm_base <- length_lm
length_lm <- vif_drop_zoo




#####  Paired-down Model Summary  ####
summary(length_lm)





#####  Fit & Residuals  ####

# Exploration of autocorrelation & residuals
length_preds <- broom::augment(length_lm, mod_dat_oisst, se_fit = T) %>% 
  mutate(.fit_hi = .fitted + 1.96*.se.fit,
         .fit_lo = .fitted - 1.96*.se.fit)



# Prediction vs Observed
ggplot(length_preds) +
  geom_ribbon(aes(year, ymin = .fit_lo, ymax = .fit_hi, fill = "Fitted"), alpha = 0.2, show.legend = F) +
  geom_line(
    data = model_data, aes(year, med_len, color = "Observed")) +
  geom_line(aes(year, .fitted, color = "Fitted")) +
  facet_wrap(~survey_area) +
  labs(y = "Median Length (cm)", title = "Median Length ~ Covariates")



# Raw Residuals
ggplot(length_preds) +
  geom_col(aes(year, .resid), fill = gmri_cols("blue"), alpha = 0.4) +
  geom_hline(yintercept = 0) +
  facet_wrap(~survey_area) +
  labs(y = "residuals", title = "Median Length ~ Covariate Model")




# Plot the residuals with acf function
length_preds %>% 
  split(.$survey_area) %>% 
  map_dfr(function(x){
    x <- arrange(x, year)
    get_ccf_vector(x$.resid, y = x$year, lagmax = 10)}, 
    .id = "survey_area") %>% 
  ggplot() +
  geom_col(aes(lag, acf), fill = gmri_cols("blue"), alpha = 0.4) +
  geom_hline(aes(yintercept = signeg), linetype = 2) +
  geom_hline(aes(yintercept = sigpos), linetype = 2) +
  facet_wrap(~survey_area) +
  labs(y = "ACF", title = "Median Length ~ Covariate Model\nACF on Residuals")







# Partial Dependence Plots

# Shared
plot(ggpredict(length_lm,  ~ survey_area))        
plot(ggpredict(length_lm,  ~ year * survey_area)) 
# plot(ggpredict(length_lm,  ~ landings * survey_area))
# plot(ggpredict(length_lm,  ~ zp_large * survey_area))

# OISST
#plot(ggpredict(length_lm,  ~ sst * survey_area))
plot(ggpredict(length_lm,  ~ gsi_oisst_resid))


# #ERSST
# plot(ggpredict(length_lm,  ~ ersst_anom * survey_area))
# plot(ggpredict(length_lm,  ~ gsi_ersst_resid))







#### b.) Auto-Correlative Error Structure  ####


# Going to make a new model with auto-correlative structure

#corAR1(0, ~ year | survey_area)
# its NA sensitive so drop explicitly


# Make new model w/ autoregressive covariance structure

# OISST Model - Max Resolution
length_ar <- gls(
  med_len ~ survey_area * (year + zp_large + landings + sst) + gsi_oisst_resid, 
  data = mod_dat_oisst,
  correlation = corAR1(0, form = ~ 1 | survey_area))

# # ERSST Model - More years, lower sst resolution
# length_ar <- gls(
#   med_len ~ survey_area * (year + zp_large + landings + ersst_anom) + gsi_ersst_resid, 
#   data = model_data %>% drop_na(avg_len, med_len, ersst_anom, landings),
#   correlation = corAR1(0, form = ~ 1 | survey_area))

# # Taking out ZP would give us the most years

summary(length_ar)



##### VIF  ####
mod_vif <- car::vif(length_ar, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)

vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# so basically everything this time
# survey_area:year, survey_area, year, landings, sst, zp_large, survey_area:landings



# Sequential removal: Not going to drop region:year because we need them
# From Zuur 2010
# collinearity has been removed by sequentially deleting each 
# variable for which the VIF value was highest until all remaining VIFs were below 3.

# Drop landings
vif_drop <-  gls(
  med_len ~ survey_area * (year + zp_large + sst) + gsi_oisst_resid, 
  data = mod_dat_oisst,
  correlation = corAR1(0, form = ~ 1 | survey_area))
# re-check vif
mod_vif <- car::vif(vif_drop, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# Drop sst
vif_drop <-  gls(
  med_len ~ survey_area * (year + zp_large) + gsi_oisst_resid, 
  data = mod_dat_oisst,
  correlation = corAR1(0, form = ~ 1 | survey_area))

# re-check vif
mod_vif <- car::vif(vif_drop, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# Drop zp_large
vif_drop <-  gls(
  med_len ~ survey_area * (year) + gsi_oisst_resid, 
  data = mod_dat_oisst,
  correlation = corAR1(0, form = ~ 1 | survey_area))

# re-check vif
mod_vif <- car::vif(vif_drop, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)





##### Paired-down Model Summary  ####

# swap model over so we can keep the code below consistent
length_ar_base <- length_ar
length_ar <- vif_drop


# Summary
summary(length_ar)



##### Fit & Residuals  ####
length_ar_preds <- bind_cols(
  mod_dat_oisst,
  data.frame(
    .fitted = length_ar$fitted,
    .resid = length_ar$residuals,
    .se.fit = attributes(length_ar$residuals)$std)) %>% 
  mutate(
    .resid = med_len - .fitted,
    .fit_hi = .fitted + 1.96*.se.fit,
    .fit_lo = .fitted - 1.96*.se.fit)




# Raw Residuals
ggplot(length_ar_preds) +
  geom_col(aes(year, .resid), fill = gmri_cols("gmri green"), alpha = 0.4) +
  geom_hline(yintercept = 0) +
  facet_wrap(~survey_area) +
  labs(y = "residuals", title = "Median Length ~ Covariate Model + AR1")


# Plot the residuals with acf function
length_ar_preds %>% 
  split(.$survey_area) %>% 
  map_dfr(function(x){
    x <- arrange(x, year)
    get_ccf_vector(x$.resid, y = x$year, lagmax = 10)}, 
    .id = "survey_area") %>% 
  ggplot() +
  geom_col(aes(lag, acf), fill = gmri_cols("gmri green"), alpha = 0.4) +
  geom_hline(aes(yintercept = signeg), linetype = 2) +
  geom_hline(aes(yintercept = sigpos), linetype = 2) +
  facet_wrap(~survey_area) +
  labs(y = "ACF", title = "Median Length ~ Covariate Model + AR\nACF on Residuals")



#### Model Predictions:


# Prediction vs Observed
ggplot(length_ar_preds) +
  geom_ribbon(aes(year, ymin = .fit_lo, ymax = .fit_hi, fill = "Fitted"), alpha = 0.2, show.legend = F) +
  geom_line(aes(year, med_len, color = "Observed")) +
  geom_line(aes(year, .fitted, color = "Fitted")) +
  facet_wrap(~survey_area) +
  labs(y = "Average Length (cm)", title = "Avg Length ~ Covariates")




# Partial Dependence Plots
# Shared
plot(ggpredict(length_ar,  ~ survey_area))        
plot(ggpredict(length_ar,  ~ year * survey_area)) 
# plot(ggpredict(length_ar,  ~ landings * survey_area))
# plot(ggpredict(length_ar,  ~ zp_large * survey_area))

# OISST
# plot(ggpredict(length_ar,  ~ sst * survey_area))
plot(ggpredict(length_ar,  ~ gsi_oisst_resid))


# #ERSST
# plot(ggpredict(spectra_lm,  ~ ersst_anom * survey_area))
# plot(ggpredict(spectra_lm,  ~ gsi_ersst_resid))


#### c.) Rolling SST & Landings  ####


# Make a model with what we think we need:
# ar(1) correlation
# some smoothed averages of SST anomalies
# use ersst so we don't lose years of data
rolling_avg_lm <- gls(
  med_len ~ survey_area * (year + zp_large + land_5 + ersst_5) + gsi_ersst_resid,
  data = mod_dat_ersst,
  correlation = corAR1(0, form = ~ 1 | survey_area))







#####  VIF  ####

# Resource https://stacyderuiter.github.io/s245-notes-bookdown/collinearity-and-multicollinearity.html
# If VIF (or squared scaled GVIF) is greater than 4, then there’s a problem and you should probably try to fix it;
# Drop landings

mod_vif <- car::vif(rolling_avg_lm, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)

vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)
# literally everything!

# Drop landings
vif_drop <-  gls(
  med_len ~ survey_area * (year + zp_large + ersst_5) + gsi_ersst_resid, 
  data = mod_dat_ersst,
  correlation = corAR1(0, form = ~ 1 | survey_area))
# re-check vif
mod_vif <- car::vif(vif_drop, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# Drop sst_5
vif_drop <-  gls(
  med_len ~ survey_area * (year + zp_large) + gsi_ersst_resid, 
  data = mod_dat_ersst,
  correlation = corAR1(0, form = ~ 1 | survey_area))

# re-check vif
mod_vif <- car::vif(vif_drop, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)

# Drop zp_large
vif_drop <-  gls(
  med_len ~ survey_area * (year) + gsi_ersst_resid, 
  data = mod_dat_ersst,
  correlation = corAR1(0, form = ~ 1 | survey_area))

# re-check vif
mod_vif <- car::vif(vif_drop, type = "predictor") %>% 
  as.data.frame() %>% 
  rownames_to_column("predictor") %>% 
  mutate(sq_scaled_vif = `GVIF^(1/(2*Df))`^2)
vif_probs <- mod_vif %>% filter(sq_scaled_vif > 4)  %>% arrange(desc(sq_scaled_vif))
vif_probs %>% pull(predictor)






#####  Paired-down Model Summary  ####

# swap model over so we can keep the code below consistent
length_rolling_base <- rolling_avg_lm
length_all_methods <- vif_drop



# Summary
summary(length_ar)



##### Fit & Residuals  ####
length_rolling_preds <- bind_cols(
  mod_dat_ersst,
  data.frame(
    .fitted = length_all_methods$fitted,
    .resid = length_all_methods$residuals,
    .se.fit = attributes(length_all_methods$residuals)$std)) %>% 
  mutate(
    .resid = med_len - .fitted,
    .fit_hi = .fitted + 1.96*.se.fit,
    .fit_lo = .fitted - 1.96*.se.fit)




# Raw Residuals
ggplot(length_rolling_preds) +
  geom_col(aes(year, .resid), fill = gmri_cols("gmri green"), alpha = 0.4) +
  geom_hline(yintercept = 0) +
  facet_wrap(~survey_area) +
  labs(y = "residuals", title = "Median Length ~ Covariate Model + AR1")


# Plot the residuals with acf function
length_rolling_preds %>% 
  split(.$survey_area) %>% 
  map_dfr(function(x){
    x <- arrange(x, year)
    get_ccf_vector(x$.resid, y = x$year, lagmax = 10)}, 
    .id = "survey_area") %>% 
  ggplot() +
  geom_col(aes(lag, acf), fill = gmri_cols("gmri green"), alpha = 0.4) +
  geom_hline(aes(yintercept = signeg), linetype = 2) +
  geom_hline(aes(yintercept = sigpos), linetype = 2) +
  facet_wrap(~survey_area) +
  labs(y = "ACF", title = "Median Length ~ Covariate Model + AR\nACF on Residuals")



#### Model Predictions:


# Prediction vs Observed
ggplot(length_rolling_preds) +
  geom_ribbon(aes(year, ymin = .fit_lo, ymax = .fit_hi, fill = "Fitted"), alpha = 0.2, show.legend = F) +
  geom_line(aes(year, med_len, color = "Observed")) +
  geom_line(aes(year, .fitted, color = "Fitted")) +
  facet_wrap(~survey_area) +
  labs(y = "Average Length (cm)", title = "Avg Length ~ Covariates")




# Partial Dependence Plots
# Shared
plot(ggpredict(length_all_methods,  ~ survey_area))        
plot(ggpredict(length_all_methods,  ~ year * survey_area)) 



# #ERSST
# plot(ggpredict(spectra_lm,  ~ ersst_anom * survey_area))
# plot(ggpredict(spectra_lm,  ~ gsi_ersst_resid))









####  Figures 2-3 - Partial Dependence Plots  ####

# Take whichever candidate model we feel most confident about
# make the partial dependence figures for the length models





####___________________#####
####___________________#####



####  Size Spectrum Models  ####



##### a.) Base Model  ####
# Spectrum Slope ~ region * (year + drivers)


# Data for models
# Drop NA's relevant to the sst data we want ot use
mod_dat_oisst <- model_data %>% drop_na(avg_len, med_len, sst, landings) %>% 
  arrange(year, survey_area) 

# Drop NA's with ERSST
mod_dat_ersst <- model_data %>% drop_na(avg_len, med_len, ersst_anom, landings, zp_large) %>% 
  arrange(year, survey_area) 



# Base mod - oisst
spectra_lm <- lm(
  b ~ survey_area * (year + zp_large + landings + sst) + gsi_oisst_resid, 
  data = mod_dat_oisst)

# # Base mod - ersst
# spectra_lm <- lm(
#   b ~ survey_area * (year + zp_large + landings + ersst_anom) + gsi_ersst_resid, 
#   data = mod_dat_ersst)


 
# Check variance inflation, hard to appreciate with intentional interactions
# time * area
car::vif(spectra_lm, type = "predictor")

# Check normality of residuals
hist(resid(spectra_lm))


# Model evaluation library
res <- simulateResiduals(spectra_lm)
plot(res)


# Exploration of autocorrelation
spectra_preds <- broom::augment(spectra_lm, mod_dat_oisst, se_fit = T) %>% 
  mutate(.fit_hi = .fitted + 1.96*.se.fit,
         .fit_lo = .fitted - 1.96*.se.fit)



# Prediction vs Observed
ggplot(spectra_preds) +
  geom_ribbon(aes(year, ymin = .fit_lo, ymax = .fit_hi, fill = "Fitted"), alpha = 0.2, show.legend = F) +
  geom_line(
    data = model_data, 
    aes(year, b, color = "Observed")) +
  geom_line(aes(year, .fitted, color = "Fitted")) +
  facet_wrap(~survey_area) +
  labs(y = "Spectra Exponent (b)", title = "Spectra Slope ~ Covariates")



# Raw Residuals
ggplot(spectra_preds) +
  geom_col(aes(year, .resid), fill = gmri_cols("gmri green"), alpha = 0.4) +
  geom_hline(yintercept = 0) +
  facet_wrap(~survey_area) +
  labs(y = "residuals", title = "Spectra Slope ~ Covariate Model")




# Plot the residuals with acf function
spectra_preds %>% 
  split(.$survey_area) %>% 
  map_dfr(function(x){
    x <- arrange(x, year)
    get_ccf_vector(x$.resid, y = x$year, lagmax = 10)}, 
    .id = "survey_area") %>% 
  ggplot() +
  geom_col(aes(lag, acf), fill = gmri_cols("gmri green"), alpha = 0.4) +
  geom_hline(aes(yintercept = signeg), linetype = 2) +
  geom_hline(aes(yintercept = sigpos), linetype = 2) +
  facet_wrap(~survey_area) +
  labs(
    y = "ACF", 
    title = "Spectra Slope ~ Covariate Model\nACF on Residuals")







# Partial Dependence Plots

# Shared
plot(ggpredict(length_ar,  ~ survey_area))        #signif
plot(ggpredict(length_ar,  ~ year * survey_area)) #signif
plot(ggpredict(length_ar,  ~ landings * survey_area))
plot(ggpredict(length_ar,  ~ zp_large * survey_area))

# OISST
plot(ggpredict(length_ar,  ~ sst * survey_area))
plot(ggpredict(length_ar,  ~ gsi_oisst_resid))


#ERSST
plot(ggpredict(length_ar,  ~ ersst_anom * survey_area))
plot(ggpredict(length_ar,  ~ gsi_ersst_resid))







##### b. auto-correlative structures  ####

# Going to make a new model with auto-correlative structure

corAR1(0, ~ year | survey_area)

# its NA sensitive and order sensitive so drop explicitly above

# Make new model w/ autoregressive covariance structure
spectra_ar <- gls(
  b ~ survey_area * (year + zp_large + landings + sst) + gsi_oisst_resid, 
  data = regression_df_ar,
  correlation = corAR1(0, form = ~ 1 | survey_area))



# Look at residuals and CCF




# Save this data for Bart
# write_csv(df_lags, here::here("data/size_spectra_results/spectra_slope_model_data.csv"))


##### c.) Rolling SST & Landings  ####


# OISST
spectra_rolling <- lm(
  b ~ survey_area * (year + zp_large + land_5 + sst_5) + gsi_oisst_resid, 
  data = mod_dat_oisst)

# ERSST
spectra_rolling <- lm(
  b ~ survey_area * (year + zp_large + land_5 + sst_5) + gsi_oisst_resid, 
  data = mod_dat_ersst)


# did we resolve temporal autocorrelation
acf(resid(spectra_rolling))




# Results Exploration
# Shared
plot(ggpredict(spectra_rolling,  ~ survey_area))        #signif
plot(ggpredict(spectra_rolling,  ~ year * survey_area)) #signif
plot(ggpredict(spectra_rolling,  ~ landings * survey_area))
plot(ggpredict(spectra_rolling,  ~ zp_large * survey_area))

# OISST
plot(ggpredict(spectra_rolling,  ~ sst * survey_area))
plot(ggpredict(spectra_rolling,  ~ gsi_oisst_resid))


#ERSST
plot(ggpredict(spectra_rolling,  ~ ersst_anom * survey_area))
plot(ggpredict(spectra_rolling,  ~ gsi_ersst_resid))










#####  Figures 4-5 - Partial Dependence Plots  ####
















