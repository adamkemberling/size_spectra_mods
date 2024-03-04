# Conceptual Figures
library(tidyverse)
library(gmRi)
library(scales)




####  Foundation - Some Size Spectra  ####

# Support functions
source(here::here("R/support/sizeSpectra_support.R"))
set.seed(123)

# Define the structure
fake_spectra <- define_log2_bins(log2_min = 0, log2_max = 13, log2_increment = 1)

# Set some biomass level
fake_intercept <- 100000

# Add some fake removals
# Make some hypothetical fishing removals
fake_f <- c(rep(0,8), c(.2, .275, .375, .5, .65))
fake_release <- rev(fake_f)
fake_f_both <- c(1+fake_release[1:7], 1-fake_f[8:13])

# Add some fake prey release
fake_spectra <- fake_spectra %>% mutate(
  eq_abund = (fake_intercept*2^-left_lim),
  abund_f = eq_abund*(1-fake_f),
  abund_f_release = eq_abund*(1+fake_release),
  abund_f_both = eq_abund*fake_f_both)



# Show the idealized equilibrium spectra
ggplot(fake_spectra) +
  geom_col(aes(left_lim, eq_abund), fill = "lightgray") +
  scale_y_log10(labels = label_log(), 
                limits = c(1, 10^6),
                expand = expansion(add = c(0,0))) +
  scale_x_continuous(labels = math_format(2^.x)) +
  theme_gmri() +
  labs(
    title = "Idealized Equilibrium Abundance Size Spectra",
    x = "Body Weight (g)",
    y = "Abundance")



# Show the fishing removals explicitly
fake_spectra %>% 
  mutate(
    log2_bins = as.numeric(as.character(log2_bins)),
    is_f = ifelse(log2_bins > 7, "darkred", "lightgrey")) %>% 
ggplot() +
  geom_col(aes(left_lim, eq_abund), fill = "lightgray") +
  geom_col(aes(left_lim, abund_f, fill = I(is_f))) +
  scale_y_log10(labels = label_log(), 
                limits = c(1, 10^9),
                expand = expansion(add = c(0,0))) +
  scale_x_continuous(labels = math_format(2^.x)) +
  theme_gmri() +
  labs(
    title = "Direct Impact on Slope Via Removals",
    x = "Body Weight (g)",
    y = "Abundance")

# Show the fishing removals
fake_spectra %>% 
  mutate(
    log2_bins = as.numeric(as.character(log2_bins)),
    is_f = ifelse(log2_bins > 7, "darkred", "lightgrey"),
    is_released = ifelse(log2_bins < 6, "royalblue", "lightgray")) %>% 
  ggplot() +
  geom_col(aes(left_lim, abund_f_release, fill = I(is_released)), alpha = 0.8) +
  geom_col(aes(left_lim, abund_f), fill = "lightgray") +
  geom_col(aes(left_lim, abund_f, fill = I(is_f)), alpha = 0.6) +
  scale_y_log10(labels = label_log(), 
                limits = c(1, 10^9),
                expand = expansion(add = c(0,0))) +
  scale_x_continuous(labels = math_format(2^.x)) +
  theme_gmri() +
  labs(
    title = "Direct Impact on Slope Via Removals,\nIndirect Effect of Prey Release",
    x = "Body Weight (g)",
    y = "Abundance")






# Try some distributions instead, don't need the log axis
fake_histo <- list(
  "x1" = data.frame("x" = 1, "size" = rnorm(100000, 1, sd = .2)),
  "x2" = data.frame("x" = 2, "size" = rnorm(10000, 2, sd = .2)),
  "x3" = data.frame("x" = 3, "size" = rnorm(1000, 3, sd = .2)),
  "x4" = data.frame("x" = 4, "size" = rnorm(100, 4, sd = .2))) %>% 
  bind_rows(.id = "group") 
  

#devtools::install_github("nschiett/fishualize", force = TRUE)
library(fishualize)
fake_histo %>% 
  ggplot(aes(x = size, group = group)) +
  geom_histogram(aes(fill = group), alpha = 0.35) +
  scale_fill_gmri() +
  scale_y_log10(labels = label_log(base = 10)) +
  scale_x_continuous(labels = math_format(2^.x)) +
  labs(y = "abundance", x = "size") +
  theme_gmri()




# Panel a. hypothetical size spectrum, with fish overlaid to indicate its a community


# panel b.two slopes, one with big F one with little F
# with panel b., you can do a mechanism figure that represents some
# literature accurate dome or knife edgle selectivity


# Panel c: temperature influence on steepness
# temperature mechanism is TSR based, so metabolism follow a curve, and metabolism plateus with size
# metabolism ~ mass as a figure showing mechanism





####  a. Temperature Impacts ####
# How Temperature impacts Spectra:

library(rTPC)

# Temperature Mechanism:
# thermal performance curve


####  b. Fishing  ####
# How does fishing impact spectra




