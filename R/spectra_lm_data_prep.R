####  Reshaping the results from size spectra and the results from median sizes
# with regional covariates
# moving all the prep code here to focus on models over at "gom_spectra_driver_lm.R"


####  Packages  ####
#library(EnvCpt)
{
  library(targets)
  library(here)
  library(gmRi)
  library(tidyverse)
  
}

# Package Conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# # Support functions
# source(here("R/support/sizeSpectra_support.R"))

# Resource Path
res_path <- gmRi::cs_path("res")

# Set a seed for reproducing stochastic elements
set.seed(123)


# levels for faceting areas
area_levels <- c("Northeast Shelf", "GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Northeast Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# table to join for swapping shorthand for long-hand names
area_df <- data.frame(
  area = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "All"),
  survey_area = c("SS", "GoM", "GB", "SNE", "MAB", "Northeast Shelf"),
  area_titles = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"))





####  Load Data  ####

##### 1. Abundance Data  ####
tar_load(catch_log2_labelled)  

# rename and format
# Add the area titles
catch_size_bins <- catch_log2_labelled %>% 
  left_join(area_df)


# 3. Total Abundance and Biomass

# Total abundance and biomass for the whole shelf
shelf_summary <- catch_size_bins %>% 
  mutate(area_titles = "Northeast Shelf") %>% 
  group_by(Year, area_titles) %>% 
  summarise(
    lwbio_sum = sum(strat_total_lwbio_s, na.rm = T),
    abund_sum = sum(strat_total_abund_s, na.rm = T),
    .groups = "drop") %>% 
  mutate(
    lwbio_mill = lwbio_sum / 1e6,
    abund_mill = abund_sum / 1e6)



# Total Abundance and Biomass for the regions
region_summary <- catch_size_bins %>% 
  group_by(Year, area_titles) %>% 
  summarise(
    lwbio_sum = sum(strat_total_lwbio_s, na.rm = T),
    abund_sum = sum(strat_total_abund_s, na.rm = T),
    .groups = "drop") %>% 
  mutate(
    lwbio_mill = lwbio_sum / 1e6,
    abund_mill = abund_sum / 1e6,
    area_titles = factor(area_titles, levels = area_levels_long))



# Combine the full area to the regional summaries
all_totals <- bind_rows(shelf_summary, region_summary) %>% 
  mutate(area_titles = factor(area_titles, area_levels_long)) %>% 
  group_by(area_titles) %>% 
  mutate(
    mean_abund = mean(abund_sum),
    mean_lwbio = mean(lwbio_sum),
    sd_abund = sd(abund_sum),
    sd_lwbio = sd(lwbio_sum),
    abund_z = (abund_sum - mean_abund) / sd_abund,
    lwbio_z = (lwbio_sum - mean_lwbio) / sd_lwbio)




##### 2.  Load Community Indices & Spectra Slope Estimates  ####
tar_load(size_spectrum_indices)

# Full Shelf
shelf_indices <- filter(size_spectrum_indices, `group ID` == "single years") %>% 
  mutate(survey_area = "Northeast Shelf")

# Sub-Regions
region_indices <- size_spectrum_indices  %>% 
  filter(`group ID` == "single years * region") %>% 
  mutate(
    yr = as.numeric(as.character(Year)),
    survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB")))


##### 3. Body Size Indices  ####
tar_load(mean_sizes_ss_groups)


# And do it again for body size changes

# Full shelf
shelf_sizes <- mean_sizes_ss_groups  %>% 
  filter(`group ID` == "single years") %>% 
  mutate(
    yr = as.numeric(as.character(Year)),
    survey_area = "Northeast Shelf")

# Regions
region_sizes <- mean_sizes_ss_groups  %>% 
  filter(`group ID` == "single years * region") %>% 
  mutate(
    yr = as.numeric(as.character(Year)),
    survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB")))



# # Export for Bart
# write_csv(region_indices, here::here("data/multispecies_spectra_results.R"))

# Reshape spectra
spectra_long <- bind_rows(shelf_indices, region_indices) %>% 
  left_join(area_df) %>% 
  select(year = Year, 
         area_titles,
         value = b) %>% 
  mutate(metric = "Size Spectra Slope",
         year = as.numeric(year))


# Reshape body sizes
sizes_long <- bind_rows(shelf_sizes, region_sizes) %>% 
  left_join(area_df) %>% 
  select(year = Year, 
         area_titles,
         "Average Weight (kg)" = mean_wt_kg,
         "Average Length (cm)" = mean_len_cm,
         "Median Weight (kg)" = med_wt_kg,
         "Median Length (cm)" = med_len_cm) %>% 
  pivot_longer(
    cols = c(3:6), 
    names_to = "metric", values_to = "value") %>% 
  mutate(year = as.numeric(year))





##### 4.  Spectra Covariate Datasets  ####

# These are the environmental drivers that are hypothesized to matter


# scaled by mean+sd within region
drivers_scaled <- read_csv(here::here("data/env_drivers_scaled.R"))

# Or the raw versions
# This is the unscaled predictor dataset
unscaled_drivers <- read_csv(here::here("data/unscaled_spectra_predictor_df.R"))






####___________________####

####  Prepare Regression Dataframe(s)  ####




##### A. Size Spectrum Slope DF  ####

# Use scaled drivers for the regressions
# Need to drop zooplankton/gsi because they are 
# non-unique with regions since we split the MAB EPU and repeated it
regression_df <- drivers_scaled %>% 
  select(!ends_with(c("zp_small","zp_large", "stream_index")))


# Build back out the data into a useful form:
# Want the size spectra data and the covariates in their respective columns
# Identify what region is associated with both: spectrum values & driver values
regression_df <- regression_df %>% 
  # Pull out spectra info
  pivot_longer(
    names_to = "spectra_param", 
    values_to = "spectra_values", 
    cols = ends_with("slope") | ends_with("int")) %>% 
  # pull out env covariates
  pivot_longer(
    names_to = "driver_var", 
    values_to = "driver_values", 
    cols = -matches("spectra|year")) %>% 
  mutate(
    year = as.numeric(year),
    # A. Flag what the driver type was
    driver_type = case_when(
      str_detect(driver_var, "landings") ~ "landings",
      str_detect(driver_var, "sst") ~ "sst",
      str_detect(driver_var, "index") ~ "gsi",
      TRUE ~ "Missed Something"),
    # B. Flag what the spectrum feature was
    param_feature = case_when(
      # Flag what the Size Distribution Parameter was
      str_detect(spectra_param, "int") ~ "Spectra Intercept",
      str_detect(spectra_param, "isd") ~ "ISD Exponent",
      str_detect(spectra_param, "slope") ~ "Spectra Slope",
      TRUE ~ "Missed Something"),
    # C. These are the areas associated with the Spectra Features
    survey_area = case_when(
      str_detect(spectra_param, "All") ~ "All",
      str_detect(spectra_param, "Georges") ~ "Georges Bank",
      str_detect(spectra_param, "Gulf") ~ "Gulf of Maine",
      str_detect(spectra_param, "Southern") ~ "Southern New England",
      str_detect(spectra_param, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    survey_area = factor(survey_area, levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")),
    
    # D. These are the areas of the drivers
    driver_area = case_when(
      str_detect(driver_var, "All") ~ "All",
      str_detect(driver_var, "Georges") ~ "Georges Bank",
      str_detect(driver_var, "Gulf") ~ "Gulf of Maine",
      str_detect(driver_var, "Southern") ~ "Southern New England",
      str_detect(driver_var, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    driver_area = factor(
      driver_area, 
      levels = c("All", "Gulf of Maine", "Georges Bank", 
                 "Southern New England", "Mid-Atlantic Bight")),
    # E. Add Decade in
    decade = floor_decade(year)) 


# Filter it so predictors are only within matching regions, 
# Don't care about cross-region relationships
regression_df <- regression_df %>% 
  filter(driver_area == survey_area) %>% 
  arrange(year, driver_var)  %>% 
  select(-c(spectra_param, driver_var)) %>% 
  pivot_wider(names_from = "driver_type", values_from = "driver_values") %>% 
  pivot_wider(names_from = "param_feature", values_from = "spectra_values") %>% 
  filter(driver_area != "All")

# GSI is so annoying here, just rejoin it without the area column
gsi <- ecodata::gsi  %>% 
  mutate(year = as.numeric(str_sub(Time, 1, 4))) %>% 
  filter(EPU == "All") %>% 
  # Put on annual scale
  group_by(year) %>% 
  summarise(
    gsi = mean(Value, na.rm = T),
    .groups = "drop")



# Load zooplankon index in from ecodata as well
# Reformat to match the other indices
zp_sli <- ecodata::zoo_sli_anom %>% 
  filter(EPU != "SS") 

# EPUS combine SNE and MAB, we can repeat the values here
zp_index <- bind_rows(
  zp_sli,
  filter(zp_sli, EPU == "MAB") %>% mutate(EPU = "SNE")) %>% 
  mutate(
    survey_area = case_when(
      EPU == "GB"  ~ "Georges Bank",
      EPU == "GOM" ~ "Gulf of Maine",
      EPU == "MAB" ~ "Mid-Atlantic Bight",
      EPU == "SNE" ~ "Southern New England"),
    Units = str_c("zp_", Var)) %>% 
  select(year = Time, survey_area, Value, Units) %>% 
  pivot_wider(names_from = Units, values_from = Value)


# Join stupid zooplankton and GSI back in...
regression_df <- regression_df %>% 
  left_join(gsi, by = c("year") ) %>% 
  mutate(survey_area = fct_drop(survey_area)) %>% 
  left_join(
    zp_index , 
    by = join_by(year, survey_area)) %>% 
  select(
    year,
    survey_area,
    b = `ISD Exponent`,
    sst,
    gsi,
    zp_small,
    zp_large,
    landings)





##### B. Body Size DF ####

# # # Size Data Regression dataframe

# Join the spectra dataframe with body sizes information
size_regression_df <- sizes_long %>% 
  #mutate(metric = ifelse(str_detect(metric, "Weight"), "avg_wt", "avg_len")) %>% 
  mutate(metric = case_when(
    metric ==  "Average Weight (kg)" ~ "avg_wt", 
    metric ==  "Average Length (cm)" ~ "avg_len", 
    metric ==  "Median Weight (kg)" ~ "med_wt", 
    metric ==  "Median Length (cm)" ~ "med_len")) %>% 
  pivot_wider(names_from = "metric", values_from = "value") %>% 
  rename(survey_area = area_titles) %>% 
  left_join(regression_df) %>% 
  filter(survey_area != "Northeast Shelf")





#### bring in ERSST5?  ####

# Yea F it
# Gives us better temporal coverage so we can see when the changes occur
# Coverage gap - gains 5 years
size_regression_df %>% 
  pivot_longer(cols = c(8:12), names_to = "env_covariate", values_to = "vals") %>% 
  mutate(flag_sst = ifelse(env_covariate == "sst", "darkred", "gray70")) %>% 
  filter(survey_area == "Gulf of Maine") %>% 
  ggplot(aes(year, vals)) +
  geom_line(linewidth = 1, aes(color = I(flag_sst))) +
  facet_wrap(~env_covariate, ncol = 1) +
  scale_x_continuous(breaks = seq(1970, 2020, 2))


# Load ersst
ersst_path <- cs_path("res", "ERSSTv5")
shape_names <- str_c(
  ersst_path, "ERSSTv5_anom_nmfs_trawl_", 
  c("gulf_of_maine", "georges_bank", "mid_atlantic_bight", "southern_new_england"), 
  ".csv")

# Load them
ersst_ts <- map(shape_names, ~read_csv(.x)) %>% 
  setNames(c("Gulf of Maine", "Georges Bank", "Mid-Atlantic Bight", "Southern New England")) %>% 
  bind_rows(.id = "survey_area") %>% 
  select(survey_area, time, sst = area_wtd_sst, sst_anom = area_wtd_sst_anom) %>% 
  group_by(year = lubridate::year(time), survey_area) %>% 
  summarise(ersst = mean(sst),
            ersst_anom = mean(sst_anom))


# Add ERSST to regression df
size_regression_df <- size_regression_df %>% 
  left_join(ersst_ts)





#### Multi-Collinearity Checks  ####


# OISST correction
plot(size_regression_df$sst, size_regression_df$gsi, main = "GSI ~ SST")
abline(0,1, col = "darkred")
cor.test(size_regression_df$sst, size_regression_df$gsi)
plot(size_regression_df$zp_large, size_regression_df$gsi, main = "GSI ~ Large Zooplankton")
abline(0,1, col = "darkred")
cor.test(size_regression_df$zp_large, size_regression_df$gsi)


# sst & gsi highly correlated
# we want* to keep sst
mod1 <- lm(gsi ~ sst, data = drop_na(size_regression_df, sst))
mod1_resid <- resid(mod1) # represents variance in GSI not accounted for by sst
plot(mod1_resid ~ drop_na(size_regression_df, sst)$sst, main = "GSI ~ SST Residuals")
abline(0,0, col = "darkred")
plot(mod1_resid ~ drop_na(size_regression_df, sst)$year, main = "GSI ~ SST Residuals Timeline")



# Include these in the model instead of gsi
mod1_resid <- data.frame(
  year = drop_na(size_regression_df, sst)$year,
  survey_area = drop_na(size_regression_df, sst)$survey_area,
  gsi_oisst_resid = mod1_resid)


# Add them back
model_data <- left_join(size_regression_df, mod1_resid)




# Do for ERSST?
plot(size_regression_df$ersst_anom, size_regression_df$gsi, main = "GSI ~ ERSST Anomalies")
abline(0,1, col = "darkred")
cor.test(size_regression_df$ersst_anom, size_regression_df$gsi)
mod2 <- lm(gsi ~ ersst_anom, data = drop_na(size_regression_df, ersst))
mod2_resid <- resid(mod2) # represents variance in GSI not accounted for by sst
plot(mod2_resid ~ drop_na(size_regression_df, ersst)$year, main = "GSI ~ ERSST Residuals Timeline")
mod2_resid <- data.frame(
  year = drop_na(size_regression_df, ersst)$year,
  survey_area = drop_na(size_regression_df, ersst)$survey_area,
  gsi_ersst_resid = mod2_resid)


# Add them back
model_data <- left_join(size_regression_df, mod1_resid)
model_data <- left_join(model_data, mod2_resid)





#####  Rolling Averages  ####

# hypothesis, sustained landings, or sustqained SST matter more than yearly calues
model_data <- model_data %>% 
  group_by(survey_area) %>% 
  arrange(survey_area, year) %>% 
  mutate(
    land_5 = zoo::rollapply(landings, 5, mean, na.rm = T, align = "right",  fill = NA),
    sst_5  = zoo::rollapply(sst, 5, mean, na.rm = T, align = "right",  fill = NA),
    ersst_5  = zoo::rollapply(ersst_anom, 5, mean, na.rm = T, align = "right",  fill = NA)
  )





####  Export  ####
write_csv(model_data, here::here("data/size_and_spectra_model_data.R"))
write_csv()







####_____________________####
#####____  Old EDA  _____####  
####___________________####

####  Supplementals Figures  ####

#### Fig 1. Community Indices Figure  ####

# Put Spectra slope, overall abundance, overall biomass, and body size changes on same figure
# Add some decadal flare

# Reshape body sizes
sizes_long <- bind_rows(shelf_sizes, region_sizes) %>% 
  left_join(area_df) %>% 
  select(year = Year, 
         area_titles,
         "Average Weight (kg)" = mean_wt_kg,
         "Average Length (cm)" = mean_len_cm) %>% 
  pivot_longer(
    cols = starts_with("Average"), 
    names_to = "metric", values_to = "value") %>% 
  mutate(year = as.numeric(year))

# Reshape abundance and Biomss
totals_long <- all_totals %>% 
  select(year = Year,
         area_titles,
         "Total Abundance" = abund_mill,
         "Total Biomass"   = lwbio_mill) %>% 
  pivot_longer(cols = starts_with("total"), names_to = "metric", values_to = "value") %>% 
  mutate(year = as.numeric(year),
         value = value*1e6)



# Put the indicators together
indicators_all <- bind_rows(list(spectra_long, sizes_long, totals_long)) %>% 
  mutate(area_titles = factor(area_titles, area_levels_long)) 


indicators_decades <- indicators_all %>% 
  mutate(decade = floor_decade(year)) %>% 
  group_by(area_titles, metric,decade) %>% 
  summarise(value = mean(value),.groups = "drop") %>% 
  mutate(year = as.numeric(as.character(decade))+5)

indicators_rolling <- indicators_all %>% 
  group_by(area_titles, metric) %>% 
  arrange(year) %>% 
  mutate(
    value = zoo::rollapply(value, 10, mean, na.rm = T, align = "center",  fill = NA)
  )


# Plot
ggplot() +
  geom_point(
    data = filter(indicators_all, area_titles != "Northeast Shelf"), 
    aes(year, value), size = 1, alpha = 0.35) +
  # geom_line(
  #   data = filter(indicators_decades, area_titles != "Northeast Shelf"),
  #   aes(year, value), linewidth = 1) +
  geom_line(
    data = filter(indicators_rolling, area_titles != "Northeast Shelf"),
    aes(year, value), linewidth = 1.25) +
  facet_grid(metric~area_titles, 
             scales = "free_y",
             labeller = labeller(metric = label_wrap_gen(width = 10)))  +
  scale_y_continuous(labels = label_comma()) + 
  theme_gmri() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype = 3, colour = "#AEAEAE")) +
  labs(x = "Year", y = "Community Metric")








####  Fig 2. Major Driver Timeseries  ####

# Do one figure where SSt/Landings/Zooplankton are included

# B. All Drivers and Slope

all_slope <- region_indices %>% 
  select(survey_area,
         year = yr, 
         value = b) %>% 
  mutate(metric = "Size Spectra Slope",
         driver_area = case_when(
           survey_area == "GoM" ~ "Gulf of Maine",
           survey_area == "GB" ~ "Georges Bank",
           survey_area == "SNE" ~ "Southern New England",
           survey_area == "MAB" ~ "Mid-Atlantic Bight"
         ))

all_drivers <- unscaled_drivers %>% 
  filter(metric != "Spectra Slope") %>% 
  mutate(
    survey_area = case_when(
      driver_area == "Gulf of Maine" ~"GoM",
      driver_area == "Georges Bank" ~ "GB",
      driver_area == "Southern New England" ~ "SNE",
      driver_area == "Mid-Atlantic Bight" ~ "MAB"
    )
  )

# Join for plot
area_levels_long <- c("Northeast Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")
all_slope <- all_slope %>% 
  bind_rows(all_drivers) %>% 
  mutate(
    metric = ifelse(metric == "zp_large", "Large Zooplankton Index", metric),
    metric = ifelse(metric == "zp_small", "Small Zooplankton Index", metric),
    metric = factor(
      metric,
      levels = c(
        "Size Spectra Slope", "Commercial Landings",
        "GSI", "SST", "Large Zooplankton Index", "Small Zooplankton Index")),
    driver_area = factor(driver_area, levels = area_levels_long))





# Driver Figure
driver_ts_all <- all_slope %>% 
  filter(metric != "Size Spectra Slope") %>% 
  ggplot() +
  geom_line(aes(year, value), linewidth = 1.25) +
  facet_grid(metric~driver_area, scales = "free_y",
             labeller = labeller(metric = label_wrap_gen(width = 10))) +
  labs(x = "Year", y = "Regional Metric") +
  scale_y_continuous(labels = label_comma()) +
  scale_x_continuous(breaks = seq(160,2020,10)) +
  theme_gmri(axis.text.x = element_text(size = 12)) +
  theme(panel.border = element_rect(color = "black", fill = NA))

# and plot
driver_ts_all





####  Fig 3. Driver CCF Figure  ####



# do a bunch of reshaping
driver_ccf_prep <- drivers_scaled %>% 
  select(-All_sst_anom) %>% 
  #rownames_to_column(var = "year") %>% 
  pivot_longer(names_to = "spectra_param", 
               values_to = "spectra_values", 
               cols = ends_with("slope") | ends_with("int")) %>% 
  pivot_longer(names_to = "driver_var", 
               values_to = "driver_values", 
               cols = -matches("spectra|year")) %>% 
  mutate(
    # C. These are the areas associated with the Spectra Features
    spectra_area = case_when(
      str_detect(spectra_param, "All")          ~ "All",
      str_detect(spectra_param, "Georges")      ~ "Georges Bank",
      str_detect(spectra_param, "Gulf")         ~ "Gulf of Maine",
      str_detect(spectra_param, "Southern")     ~ "Southern New England",
      str_detect(spectra_param, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    spectra_area = factor(spectra_area, 
                          levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")),
    
    # D. These are the areas of the drivers
    driver_area = case_when(
      str_detect(driver_var, "All")          ~ "All",
      str_detect(driver_var, "Georges")      ~ "Georges Bank",
      str_detect(driver_var, "Gulf")         ~ "Gulf of Maine",
      str_detect(driver_var, "Southern")     ~ "Southern New England",
      str_detect(driver_var, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    driver_area = factor(driver_area, 
                         levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight"))) 



# Filter it out so there is only cases where the driver area matches
driver_ccf_prep <- driver_ccf_prep %>% 
  filter(driver_area == spectra_area | driver_area == "All") %>% 
  arrange(year, driver_var)



# Walk through each xvariable, and see how it correlates with each yvar
# split on xvar
# then split on yvars
# make sure order is good
ccf_relationships <- driver_ccf_prep %>% 
  drop_na() %>% 
  split(.$spectra_param) %>% 
  map_dfr(function(x_param){
    x_param %>% 
      split(.$driver_var) %>% 
      map_dfr(function(driver_y_data){
        
        # Get the CCF Information
        ccf_df <- get_ccf_vector(
          x = driver_y_data$spectra_values,
          y = driver_y_data$driver_values,
          lagmax = 8
        )
      }, .id = "driver_var")
  }, .id = "spectra_param")




# Build back out the labels for plotting
ccf_plot_data <- ccf_relationships %>% 
  mutate(
    # Flag what the driver type was
    driver_type = case_when(
      str_detect(driver_var, "landings") ~ "Commercial Landings",
      str_detect(driver_var, "sst")      ~ "SST",
      str_detect(driver_var, "index")    ~ "GSI",
      str_detect(driver_var, "small")    ~ "ZP-S",
      str_detect(driver_var, "large")    ~ "ZP-L",
      TRUE ~ "Missed Something"),
    param_feature = case_when(
      # Flag what the Size Distribution Parameter was
      str_detect(spectra_param, "int")   ~ "Spectra Intercept",
      str_detect(spectra_param, "isd")   ~ "ISD Exponent",
      str_detect(spectra_param, "slope") ~ "Spectra Slope",
      TRUE ~ "Missed Something"),
    spectra_region = case_when(
      # Flag what region the driver was coming from
      str_detect(spectra_param, "All")          ~ "All",
      str_detect(spectra_param, "Georges")      ~ "Georges Bank",
      str_detect(spectra_param, "Gulf")         ~ "Gulf of Maine",
      str_detect(spectra_param, "Southern")     ~ "Southern New England",
      str_detect(spectra_param, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    # Flag when it crosses threshold
    sig_flag = ifelse(acf < signeg | acf > sigpos, T, F),
    # Set Factor Levels
    spectra_region = factor(spectra_region, 
                            levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")))


# Limit to one Response:
ccf_plot_data <- ccf_plot_data %>% 
  filter(param_feature == "ISD Exponent") %>% 
  filter(spectra_region != "All")


# Plot them
ccf_figure <- ccf_plot_data %>% 
  mutate(
    driver_type = ifelse(driver_type == "ZP-L", "Large Zooplankton Index", driver_type),
    driver_type = ifelse(driver_type == "ZP-S", "Small Zooplankton Index", driver_type),
    driver_type = factor(
      driver_type,
      levels = c(
        "Size Spectra Slope", "Commercial Landings",
        "GSI", "SST", "Large Zooplankton Index", "Small Zooplankton Index")),
    vjust =  ifelse(acf < 0, -0.5,1.5),
    sig_alpha = ifelse(sig_flag, 1, 0.5),
    sig_color = ifelse(sig_flag, "black", "transparent")
  ) %>% 
  #filter(spectra_region == "Gulf of Maine") %>% 
  ggplot() +
  geom_vline(xintercept = 0, color = "gray25", linewidth = 1) +
  geom_col(aes(lag, acf, fill = driver_type, color = I(sig_color), alpha = I(sig_alpha))) +
  geom_text(aes(lag, y = 0, label = lag, vjust = I(vjust), color = I(sig_color)), size = 3) +
  geom_line(aes(x = lag, y = sigpos), linetype = 2, color = "gray25") +
  geom_line(aes(x = lag, y = signeg), linetype = 2, color = "gray25") +
  geom_hline(yintercept = 0, color = "gray25", linewidth = 1) +
  #scale_color_gmri() +
  scale_fill_gmri() +
  facet_grid(
    driver_type~spectra_region,
    labeller = labeller(
      driver_type = label_wrap_gen(width = 10),
      spectra_region = label_wrap_gen(width = 10))) +
  scale_x_continuous(
    limits = c(-8, 8), 
    breaks = seq(from = -8, to = 8, by = 2),
    expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme_gmri(
    legend.position = "bottom",
    axis.text = element_text(size = 10)) +
  theme(legend.position = "none") +
  labs(x = "Driver Lag", y = "Correlation", 
       fill = "Driver Type:", color = "Driver Type:",
       title = "Cross-Correlation Function:\nSpectrum Slope ~ Regional Drivers",
       caption = "All predictors scaled over year range of 1982-2019. Lags exceeding significance threshold outlined and numbered.")




ccf_figure




####  Fig 4. Size Bin Heatmaps  ####

# This figure will track biomass across size bins through time
# Might be interesting to do spectra and the heatmap as pairs


# Take all the data:
# Get the average to normalize against
spectra_heatmap_data <- catch_size_bins %>% 
  mutate(area_titles = "Northeast Shelf") %>% 
  bind_rows(catch_size_bins) %>% 
  split(.$area_titles) %>% 
  map_dfr(function(x){
    
    # Take the data for a region
    # Group on the size bins - and get total abundance/biomass
    
    # what should we use for mean/sd??
    
    annual_bin_totals <- x %>% 
      group_by(area_titles, Year, left_lim, bin_label, bin_width) %>% 
      dplyr::summarize(
        single_year_total_abund = sum(strat_total_abund_s, na.rm = T),
        single_year_total_bio   = sum(strat_total_lwbio_s, na.rm = T),
        .groups = "drop")
    
    
    #### Standarding Abundance/Biomass Across All Size Bins
    # They Should in theory have the same biomass in each....
    # So we should be able to center on that and see everything
    overall_bin_avg <- annual_bin_totals %>% 
      group_by(area_titles) %>% 
      summarise(
        overall_avg_abund = mean(single_year_total_abund, na.rm = T),
        overall_sd_abund = sd(single_year_total_abund, na.rm = T),
        overall_avg_bio = mean(single_year_total_bio, na.rm = T),
        overall_sd_bio = sd(single_year_total_bio, na.rm = T),
      ) %>% 
      left_join(annual_bin_totals) %>% 
      mutate(
        single_year_abund_z = (single_year_total_abund - overall_avg_abund) / overall_sd_abund,
        single_year_bio_z = (single_year_total_bio - overall_avg_bio) / overall_sd_bio
      )
    
    #### Standarding Abundance/Biomass Within Size Bins
    size_bin_overall <- overall_bin_avg %>% 
      group_by(area_titles, left_lim, bin_label, bin_width) %>% 
      mutate(
        overall_sbin_abund = sum(single_year_total_abund),
        overall_sbin_bio   = sum(single_year_total_bio),
        overall_sbin_avg_abund  = mean(single_year_total_abund, na.rm = T),
        overall_sbin_sd_abund   = sd(single_year_total_abund, na.rm = T),
        overall_sbin_avg_bio    = mean(single_year_total_bio, na.rm = T),
        overall_sbin_sd_bio     = sd(single_year_total_bio, na.rm = T)) %>% 
      ungroup() %>% 
      mutate(
        sbin_abund_z = (single_year_total_abund - overall_sbin_avg_abund) / overall_sbin_sd_abund,
        sbin_bio_z   = (single_year_total_bio - overall_sbin_avg_bio) / overall_sbin_sd_bio)
    
    #return(overall_bin_avg)
    return(size_bin_overall)
  }, .id = "area_titles") %>% 
  mutate(area_titles = factor(area_titles, levels = area_levels_long))




# Plot the heatmap: Biomass - Centered within regions
spectra_heatmap_data %>% 
  ggplot(aes(
    x = Year, 
    y = 2^left_lim, 
    fill = single_year_bio_z)) +
  geom_tile(color = "white") +
  scale_x_continuous(expand = expansion(add = c(0.5,0.5))) +
  scale_y_continuous(
    trans = "log2",
    labels = label_log(base = 2),
    breaks = 2^seq(0,12, by = 2),
    expand = expansion(add = c(0.5, 0.5))) +
  scale_fill_distiller(
    limits = c(-2,2), 
    labels = c("<-4", c(-2, 0, 2), ">4"),
    oob = oob_squish, 
    palette = "RdBu") +
  facet_wrap(~area_titles)  +
  theme_gmri(
    legend.title = element_text(angle= 0, size = 12),
    panel.grid.major = element_blank(),
    legend.position = "bottom", 
    legend.direction = "horizontal") +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black", 
      ticks.colour = "black",
      barwidth = unit(10, "cm"),
      barheight = unit(1, "cm"))) +
  labs(title = "Departures from All-Size Equal Biomass Presumption",
       y = "Body Weight Bin (g)",
       fill = "Equilibrium Spectra Biomass Anomaly (z)")




# Plot the heatmap: Biomass - Centered within regions & Size Bins
heatmap_plots <- spectra_heatmap_data %>% 
  #filter(area_titles == "Gulf of Maine") %>% 
  ggplot(aes(x = Year, y = 2^left_lim, fill = sbin_bio_z)) +
  geom_tile(color = "white") +
  scale_x_continuous(expand = expansion(add = c(0.5,0.5))) +
  scale_y_continuous(
    trans = "log2",
    labels = label_log(base = 2),
    breaks = 2^seq(0,12, by = 2),
    expand = expansion(add = c(0.5, 0.5))) +
  scale_fill_distiller(
    limits = c(-2,2), 
    labels = c("<-2", c(-1, 0, 1), ">2"),
    oob = oob_squish, 
    palette = "RdBu") +
  facet_wrap(~area_titles)  +
  theme_gmri(
    legend.title = element_text(angle= 0, size = 12),
    panel.grid.major = element_blank(),
    legend.position = "bottom", 
    legend.direction = "horizontal") +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black", 
      ticks.colour = "black",
      barwidth = unit(5, "cm"),
      barheight = unit(1, "cm"))) +
  labs(title = "Changes in Size-Class Biomass",
       y = "Body Weight Bin (g)",
       fill = "Relative Biomass\nWithin Each Size-Bin\n(z-score)")


# Plot just the heatmaps
heatmap_plots



#### Pair Size Bin Change with Spectrum changes
pair_plots <- area_levels_long %>% 
  map(function(name_x){
    
    # Plot the Spectra Changes
    spectra_plot <- spectra_long %>% 
      filter(area_titles == name_x) %>% 
      ggplot(aes(year, value)) +
      geom_line(linewidth = 1) +
      scale_x_continuous(
        limits = c(1970, 2020),
        expand = expansion(add = c(0.75,0.25))) +
      #geom_point(size = 2) + 
      theme(plot.margin = margin(t = 0, b = 0, r = 0, l = 0)) +
      facet_wrap(~area_titles) +
      theme_gmri() +
      labs(x = NULL, y = "Spectra Slope")
    
    
    # Plot the heatmap
    sizes_plot <- spectra_heatmap_data %>% 
      filter(area_titles == name_x) %>% 
      ggplot(aes(x = Year, y = 2^left_lim, fill = sbin_bio_z)) +
      geom_tile(color = "white", show.legend = F) +
      scale_x_continuous(expand = expansion(add = c(0.5,0.5))) +
      scale_y_continuous(
        trans = "log2",
        labels = label_log(base = 2),
        breaks = 2^seq(0,12, by = 2),
        expand = expansion(add = c(0.5, 0.5))) +
      scale_fill_distiller(
        limits = c(-2,2), 
        labels = c("<-2", c(-1, 0, 1), ">2"),
        oob = oob_squish, 
        palette = "RdBu") +
      theme_gmri(
        axis.ticks.y = element_line(),
        panel.grid.major = element_blank()) +
      theme(plot.margin = margin(t = 0, b = 0, r = 0, l = 0)) +
      labs(y = "Body-Size Bin (g)")
    
    # Return the combo
    return(spectra_plot / sizes_plot)
  })

# Set up the legend as a plot
leg <- as_ggplot(get_legend(heatmap_plots)) +
  theme(plot.margin = margin(t = 0, b = 0, r = 0, l = 0))

# Combine Everything
spectra_pair_figure <- 
  ((pair_plots[[1]] + labs(x= NULL)) | 
     (pair_plots[[2]] + labs(x= NULL) & labs(y = NULL) ) | 
     (pair_plots[[3]] & labs(y = NULL) )) / 
  (pair_plots[[4]] | 
     (pair_plots[[5]] & labs(y = NULL)) | 
     leg) + plot_layout(widths = c(1,1,1))

# This looks much better, you can see where a bump in single tiles messes the whole thing up
spectra_pair_figure



####  Figure S1. Inter-annual Variability in Abundance  ####

# super one-off, but worth checking
all_totals %>% 
  select(
    year = Year,
    area_titles,
    abund_z,
    lwbio_z
  ) %>% 
  pivot_longer(cols = ends_with("z"),
               names_to = "metric",
               values_to = "value") %>% 
  ggplot(aes(area_titles, value, color = area_titles)) +
  geom_boxplot() +
  facet_grid(metric~.)

all_totals %>% 
  select(
    year = Year,
    area_titles,
    abund_z,
    lwbio_z) %>% 
  pivot_longer(cols = ends_with("z"),
               names_to = "metric",
               values_to = "value") %>% 
  ggplot(aes(year, value, color = area_titles)) +
  geom_line() +
  facet_grid(metric~.)







####  Stream Graph Mirages  ####
# I wonder what a streamplot would look like
library(ggstream)


# # Pretty Cool Actually
# 
# # Proportional Biomass
# spectra_heatmap_data %>% 
#   mutate(
#     bin_label = fct_reorder(bin_label, left_lim, median),
#     single_year_total_bio = (single_year_total_bio/1000)/1e6) %>% 
#   ggplot(aes(x = Year, y = single_year_total_bio, 
#              fill = fct_rev(bin_label))) +
#   geom_stream(
#     type = "proportional",
#     true_range = "none",
#     color = "white") +
#   scale_x_continuous(expand = expansion(add = c(1,1))) +
#   scale_y_continuous(labels = label_percent()) +
#   facet_wrap(~area_titles, ncol = 1, scale = "free_y") +
#   theme_gmri() +
#   labs(y = "Total Biomass",
#        fill = "Size Class")
# 
# # Absolute Units - Ridge Version
# # Confirm lwbio units
# spectra_heatmap_data %>% 
#   mutate(
#     bin_label = fct_reorder(bin_label, left_lim, median),
#     single_year_total_bio = (single_year_total_bio/1000)/1e6) %>% 
#   ggplot(aes(x = Year, y = single_year_total_bio, 
#              fill = fct_rev(bin_label))) +
#   geom_stream(
#     type = "ridge",
#     #extra_span = 0.2
#     true_range = "none",
#     color = "white"
#     ) +
#   geom_stream(
#     type = "ridge",
#     #extra_span = 0.2,
#     true_range = "none",
#     alpha = 0.3) +
#   scale_y_continuous(labels = label_comma(suffix = "M")) +
#   scale_x_continuous(expand = expansion(add = c(1,1))) +
#   facet_wrap(~area_titles, ncol = 1, scale = "free_y") +
#   theme_gmri() +
#   #scale_fill_brewer() +
#   labs(y = "Total Biomass",
#        fill = "Size Class")





####  Functional Group Stream Graph  ####

# Just need total biomass for the area-fgroup-year
# Then I can make plot, this could get at the functional group changes

fgroup_summary <- catch_size_bins %>% 
  mutate(area_titles = "Northeast Shelf") %>% 
  bind_rows(catch_size_bins) %>% 
  group_by(area_titles, year = Year, hare_group) %>% 
  summarise(
    total_abund = sum(strat_total_abund_s),
    total_bio = sum(strat_total_lwbio_s),
    .groups = "drop") %>% 
  mutate(area_titles = factor(area_titles, levels = area_levels_long))


# Abundance
abundance_stream <- fgroup_summary %>% 
  ggplot(aes(year, total_abund, fill = hare_group)) +
  geom_stream(
    type = "proportional",
    true_range = "none",
    color = "white",
    show.legend = F) +
  scale_x_continuous(expand = expansion(add = c(0,0))) +
  scale_y_continuous(labels = label_percent(),
                     expand = expansion(add = c(0,0))) +
  facet_wrap(~area_titles, nrow = 2, scale = "free_y") +
  theme_gmri() +
  scale_fill_brewer() +
  labs(y = "Total Abundance",
       x = NULL,
       title = "Abundance Composition",
       fill = "Functional Group")

# Biomass
bio_stream <- fgroup_summary %>% 
  ggplot(aes(year, total_bio, fill = hare_group)) +
  geom_stream(
    type = "proportional",
    true_range = "none",
    color = "white", show.legend = T) +
  scale_x_continuous(expand = expansion(add = c(0,0))) +
  scale_y_continuous(labels = label_percent(),
                     expand = expansion(add = c(0,0))) +
  facet_wrap(~area_titles, nrow = 2, scale = "free_y") +
  theme_gmri(legend.position = c(0.775, 0.2)) +
  scale_fill_brewer() +
  labs(y = "Total Biomass",
       x = "Year",
       title = "Biomass Composition",
       fill = "Functional Group")



# These look cool, maybe do both but reorient:
abundance_stream / bio_stream






#####__________________####
####  Dropping - Exploratory Data Analyses  ####
####  ENVCpt Breakpoint Structures  ####

# Does it make sense to do this?
# Gulf of Maine only region with breakpoints*
# But only with AR2 term, not without
# all other regions, no

# Pull Regions to test independently:
regression_df <- regression_df %>% 
  rename(
    landings = `Commercial Landings`,
    b = `ISD Exponent`)

# Gulf of Maine
gom_df <- filter(regression_df, survey_area == "Gulf of Maine") %>% 
  drop_na()



# Check each for different changepoint structures
gom_cpt <- envcpt(gom_df$b, verbose = F)



# Pulling out details of best supported model
cpt_res <- list(
  "gom" =  list("cpt_res" = gom_cpt)#,
  # "gb"  =  list("cpt_res" = gb_cpt),
  # "sne" =  list("cpt_res" = sne_cpt),
  # "mab" =  list("cpt_res" = mab_cpt)
) %>% map(function(x){
  best <- names(which.min(AIC(x$cpt_res)))
  best_summ <- x$cpt_res[[best]]
  x$best <- best
  x$best_summ <- best_summ
  return(x)
})


# # Plot GOM, the only one with changepoints
plot(gom_cpt, main = "Gulf of Maine, Break Point Analysis")



# Put the rest of the results into a delta AIC table:
bind_rows(map(cpt_res, ~tibble("Most Supported Model" = .x[["best"]])), .id = "Region") %>% knitr::kable()






####  Gulf of Maine Breakpoint Model Selection  ####


gom_df <- mutate(gom_df, 
                 regime = case_when(
                   year < 1999 ~ "yrs_1982-1998",
                   year < 2007 ~ "yrs_1999_2006",
                   year >= 2007 ~ "yrs_2007_2019"))


# Add the autoregressive predictors neatly/manually
gom_df_lag <- gom_df %>% 
  mutate(
    blag1   = dplyr::lag(b, 1),
    blag2   = dplyr::lag(b, 2),
    zpslag1 = dplyr::lag(zp_small, 1),
    zpslag2 = dplyr::lag(zp_small, 2)) %>% 
  drop_na(blag2) %>% 
  mutate(regime = factor(regime, levels = c("yrs_1982-1998", "yrs_1999_2006", "yrs_2007_2019")))


# kitchen sink
gom_ar_models = lm(
  b ~ 
    regime * SST + 
    regime * GSI + 
    regime * landings + 
    regime * zp_small + 
    regime * zp_large + 
    regime * blag1 + 
    regime * blag2, 
  data = gom_df_lag)



# Model selection on models with lags at 1+2 years
library(MuMIn)
library(broom)

# change na. action for dredge
options(na.action = "na.fail")
results_gom <- dredge(gom_ar_models)


# Best autoregressive model
best_gom <- get.models(results_gom, subset = delta == 0)[[1]]



####  "Best" Model Exploration ####

# What do the partial dependence plots look like:
# they look like trash because the interactions are NA


# Create plot-data data frame with mean values of the other predictors
# And the different factor levels too
# langings ranges from the observed minimum to maximum


# Set up the daata to plot the range of predicted values for landings across regimes
plotdata <- expand.grid(
  landings = seq(min(gom_df_lag$landings), max(gom_df_lag$landings), by = .01),
  regime = unique(gom_df_lag$regime),
  blag2 = mean(gom_df_lag$blag2, na.rm = T),
  zp_large = mean(gom_df_lag$zp_large, na.rm = T),
  zp_small = mean(gom_df_lag$zp_small, na.rm = T),
  SST = mean(gom_df_lag$SST, na.rm = T)
)



# Restrict the ranges to the ranges seen within each regime?
# Makes it easier to see the actual ranges we observed
plotdata <- plotdata %>% 
  split(.$regime) %>% 
  imap_dfr(function(x,y){
    observed_range <- filter(gom_df_lag, regime == y) %>% pull(landings)
    filter(x, landings >= min(observed_range, na.rm = T), landings <= max(observed_range, na.rm = T))
  })


# Add Best model predictions
bmod_se <- predict(best_gom, plotdata, se.fit = T)[["se.fit"]]
plotdata <- plotdata %>% 
  mutate(bmod_preds = predict(best_gom, plotdata),
         landings_se = bmod_se,
         landings_hi = bmod_preds + 1.96*landings_se,
         landings_lo = bmod_preds - 1.96*landings_se)

# Plotting Differential Impact of landings
ggplot() +
  geom_ribbon(
    data = plotdata,
    aes(ymin = landings_lo, ymax = landings_hi, x = landings, fill = regime), alpha = 0.5)+
  geom_line(
    data = plotdata,
    aes(landings, bmod_preds, color = regime), linewidth = 1) +
  facet_wrap(~regime, scales = "free") +
  geom_point(data = gom_df_lag, aes(landings, b, color = regime), size = 2) +
  scale_color_gmri() +
  scale_fill_gmri() +
  labs(x = "Landings (z)", y = "Multi-Species Spectrum Slope", caption = "Points = Observed Values")








# Do it again for large zooplankton


# Set up the daata to plot the range of predicted values for landings across regimes
plotdata <- expand.grid(
  landings = mean(gom_df_lag$landings, na.rm = T),
  regime = unique(gom_df_lag$regime, na.rm = T),
  blag2 = mean(gom_df_lag$blag2, na.rm = T),
  zp_large = seq(min(gom_df_lag$zp_large), max(gom_df_lag$zp_large), by = .01),
  zp_small = mean(gom_df_lag$zp_small, na.rm = T),
  SST = mean(gom_df_lag$SST, na.rm = T)
)



# Restrict the ranges to the ranges seen within each regime?
# Makes it easier to see the actual ranges we observed
plotdata <- plotdata %>% 
  split(.$regime) %>% 
  imap_dfr(function(x,y){
    observed_range <- filter(gom_df_lag, regime == y) %>% pull(zp_large)
    filter(x, zp_large >= min(observed_range, na.rm = T), zp_large <= max(observed_range, na.rm = T))
  })


# Add Best model predictions
bmod_se <- predict(best_gom, plotdata, se.fit = T)[["se.fit"]]
plotdata <- plotdata %>% 
  mutate(bmod_preds = predict(best_gom, plotdata),
         zp_large_se = bmod_se,
         zp_large_hi = bmod_preds + 1.96*zp_large_se,
         zp_large_lo = bmod_preds - 1.96*zp_large_se)


# Plotting Differential Impact of zp_large
ggplot() +
  geom_ribbon(
    data = plotdata,
    aes(ymin = zp_large_lo, ymax = zp_large_hi, x = zp_large, fill = regime), alpha = 0.5)+
  geom_line(
    data = plotdata,
    aes(zp_large, bmod_preds, color = regime), linewidth = 1) +
  facet_wrap(~regime, scales = "free") +
  geom_point(data = gom_df_lag, aes(zp_large, b, color = regime), size = 2) +
  scale_color_gmri() +
  scale_fill_gmri() +
  labs(x = "Large Zooplankton Index (z)", y = "Multi-Species Spectrum Slope", caption = "Points = Observed Values",
       title = "Really Terrible Fit for Large Zooplankton and Regime Interaction",
       subtitle = "*At Mean Values for Other Predictors")



#### EMMEANS Interactions Check  ####
# Use emmeans to speed life up:

# Partial Dependence plots, way less work...
emmip(best_gom, regime ~ landings, 
      at = list(
        landings = gom_df_lag$landings))

# Zooplankton Large
emmip(best_gom, regime ~ zp_large, 
      at = list(
        zp_large = gom_df_lag$zp_large
      ))




####  Top Performing Model Weights/PRedictions  ####

# Top performers
# Pull details for models with delta AIC < 2
top_mods_gom <- get.models(results_gom, subset = delta < 2)

# Peformance, Gets the AIC/AICc order 
top_perf <- map_dfr(top_mods_gom, function(x){
  broom::glance(x) %>% 
    mutate(AICc = AICc(x))}, 
  .id = "id") %>% 
  arrange(AICc)

# Calculate delta AICc
top_perf$Delta <- top_perf$AICc - min(top_perf$AICc)

# Give the ranking based on that order
top_perf <- top_perf %>% 
  arrange(Delta) %>% 
  mutate(`Model ID` = row_number(),
         `Model ID` = str_c("Model ", `Model ID`))



# Retrieve Model Coefficients from "top" models
# join in the model number based on AIC order
top_coef <- map_dfr(top_mods_gom, tidy, .id = "id") %>% 
  left_join(select(top_perf, id, `Model ID`)) %>% 
  select(-id)


# Reshape Performance Specs to bind_row later
top_perf <- select(top_perf, `Model ID`, r2 = r.squared, Delta) %>% 
  pivot_longer(cols = -1, names_to = "term", values_to = "estimate")




###  Model Weights  ####

# Sum of weights - predictor importance
sw(top_mods_gom)

# Model averaging
gom.ests <- model.avg(top_mods_gom, revised.var = TRUE)
gom.ests


# Model preds for each
model.preds <- sapply(setNames(top_mods_gom, c("mod 1", "mod 2", "mod 3")), predict, newdata = gom_df_lag)

# Weighted average prediction
out.put <- model.sel(top_mods_gom$`10459`, top_mods_gom$`2073`, top_mods_gom$`10491`)
mod.ave.preds <- model.preds %*% Weights(out.put)

gom_df_lag %>% cbind(mod.ave.preds) %>% 
  ggplot() +
  geom_line(aes(year, b, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, mod.ave.preds, color = "Model-Averaged Predicition"), linewidth = 1) +
  labs(title = "Model Averaged Prediction")




#### Model Averaged Partial Dependence Plots  ####

# More interesting application: everything but one predictor set to mean value

# langings ranges from the observed minimum to maximum
landings <- seq(min(gom_df_lag$landings), max(gom_df_lag$landings), by = .01)

# Create plot-data data frame with mean values of the other predictors
# And the different factor levels too
# Set up the daata to plot the range of predicted values for landings across regimes
plotdata <- expand.grid(
  landings = seq(min(gom_df_lag$landings), max(gom_df_lag$landings), by = .01),
  regime = unique(gom_df_lag$regime),
  blag2 = mean(gom_df_lag$blag2, na.rm = T),
  zp_large = mean(gom_df_lag$zp_large, na.rm = T),
  zp_small = mean(gom_df_lag$zp_small, na.rm = T),
  SST = mean(gom_df_lag$SST, na.rm = T)
)



# Restrict the ranges to the ranges seen within each regime?
# Makes it easier to see the actual ranges we observed
plotdata <- plotdata %>% 
  split(.$regime) %>% 
  imap_dfr(function(x,y){
    observed_range <- filter(gom_df_lag, regime == y) %>% pull(landings)
    filter(x, landings >= min(observed_range, na.rm = T), landings <= max(observed_range, na.rm = T))
  })


# Predict response for the plot data with each model
model.preds = sapply(top_mods_gom, predict, newdata = plotdata)

# Weight the prediction from each model by its AIC weight
# and sum (matrix multiplication)
mod.ave4plot <- model.preds %*% Weights(out.put)

# plot the model averaged predicted densities vs elevation

#plot(mod.ave4plot ~ plotdata$landings, type = 'l', xlab="Elevation (m)", ylab="Model averaged predicted density")

plotdata %>% 
  mutate(avg_pred = mod.ave4plot[,1]) %>% 
  ggplot() +
  geom_line(aes(landings, avg_pred, color = regime, linetype = "Model Averaged Prediction")) +
  geom_point(data = gom_df_lag, aes(landings, b, color = regime), size = 2) +
  facet_wrap(~regime, nrow = 1, scales = "free") +
  scale_color_gmri() +
  scale_fill_gmri() +
  labs(
    y = "Multi-Species Spectrum Slope", x = "Landings Index (z)", color = "Regime")




# What about with confidence intervals?
# documentation: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf
avg_preds <- predict(model.avg(top_mods_gom), plotdata, se.fit = T)
plotdata %>% 
  mutate(avg_pred = avg_preds$fit,
         pred_se = avg_preds$se.fit,
         pred_hi = avg_pred + 1.96*pred_se,
         pred_lo = avg_pred - 1.96*pred_se) %>% 
  ggplot() +
  geom_ribbon(aes(landings, ymin = pred_lo,  ymax = pred_hi, fill = regime), alpha = 0.3) +
  geom_line(aes(landings, avg_pred, color = regime, linetype = "Model Averaged Prediction")) +
  geom_point(data = gom_df_lag, aes(landings, b, color = regime), size = 2) +
  scale_color_gmri() +
  scale_fill_gmri() +
  facet_wrap(~regime, scales = "free") +
  labs(
    y = "Multi-Species Spectrum Slope", x = "Landings Index (z)", color = "Regime", fill = "Regime",
    title = "AIC Top Model Ensemble Prediction: Model-Averaging Prediction")




# # year and area emmip
# emmeans::emmip(
#   object = length_lm, 
#   formula = ~ year * survey_area, 
#   at = list(
#     "survey_area" = area_levels_long[c(2:5)], 
#     "year" = c(1982:2019)), 
#   type = "response", 
#   plotit = F) %>% 
#   ggplot(aes(year, yvar, color = survey_area)) +
#   geom_line()

