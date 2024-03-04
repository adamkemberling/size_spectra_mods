# Building Expectations for Arrhenius equation / TSR
# Credit to Dr. Bart Difiore

# Range of masses: m
m <- seq(2,1000)

# age 0 mass: a0_m
a0_m <- 1

# metabolism B as a function of body size
# 3/4 scaling:
"LINDMARK 2018: Across species, rates are often
assumed and found to scale as power functions 
of mass with exponents of 3/4 for whole organism rates, 
exponentially with temperature, and with 
independent mass and temperature effects (e.g., in
the Arrhenius fractal supply model applied in the metabolic theory
of ecology (MTE) (Brown et al., 2004; Downs et al., 2008; Gillooly
et al., 2001))."


B <- a0_m*m^0.75
plot(B ~ m)

# Arrhenius equation
ar <- function(r0 = 10, E = 0.67, k = 8.62*10^-5, temp){
  r0*exp(-1*E/(k*temp))
}
temp = seq(10,30, length.out = 100)+273.15
R <- ar(temp = temp)


plot(R ~ temp)



# Temperature Size Rule
temp_size <- function(r0 = 10, E = 0.67, k = 8.62*10^-5, temp, mass){
  r0*mass^0.75*exp(-1*E/(k*temp))
}


df <- expand.grid(temp = c(10+273.15, 30+273.15), mass = seq(2,1000))
df$metabolism <- temp_size(temp = df$temp, mass = df$mass)

library(tidyverse)
# Mass
ggplot(df, aes(x = mass, y = metabolism/mass))+
  geom_line(aes(color = as.factor(temp)))+
  scale_x_log10()+
  scale_y_log10()




sharpe_schoolfield <- function(
    c0   = 0.79, 
    E    = 0.73, 
    Eh   = 1.89, 
    k    = 8.62*10^-5, 
    Tc   = 273.15, 
    Th   = 25+273.15, 
    beta = 0.75, 
    temp, 
    mass){
  
  #
  numerator = c0*mass^beta*exp(E*(1/(k*Tc) - 1/(k*temp)))
  denominator = 1 + exp(Eh*(1/(k*Th) - 1/(k*temp)))
  
  numerator/denominator
  
}


df <- expand.grid(temp = c(10+273.15, 30+273.15), mass = seq(2,10))

df$metab <- sharpe_schoolfield(temp = df$temp, mass = df$mass)

library(tidyverse)
p1 <- df %>%
  mutate(temp_cat = ifelse(temp == 10+273.15, "cold", "warm")) %>%
  ggplot(aes(x = mass, y = metab/mass))+
  geom_line(aes(color = temp_cat))+
  scale_color_manual(values = c("blue", "red"))+
  scale_x_log10()+
  scale_y_log10()+
  theme_classic()


df <- expand.grid(temp = seq(10,30)+273.15, mass = c(2,5,10))

df$metab <- sharpe_schoolfield(temp = df$temp, mass = df$mass)

p2 <- df %>%
  mutate(mass_cat = case_when(
    mass == 2 ~ "small", 
    mass == 5 ~ "medium", 
    mass == 10 ~ "large")) %>%
ggplot(aes(x = temp-273.15, y = metab))+
  geom_line(aes(color = mass_cat))+
  scale_color_manual(values = c("darkgreen", "green", "lightgreen"))+
  theme_classic()

df %>%
  mutate(mass_cat = case_when(mass == 2 ~ "small", 
                              mass == 5 ~ "medium", 
                              mass == 10 ~ "large")) %>%
  group_by(mass_cat) %>%
  slice_max(metab) #%>% View()



cowplot::plot_grid(p1, p2, nrow = 1)




p3 <- df %>%
  mutate(mass_cat = case_when(mass == 2 ~ "small", 
                              mass == 5 ~ "medium", 
                              mass == 10 ~ "large")) %>%
  ggplot(aes(x = temp-273.15, y = metab))+
  geom_line(aes(color = mass_cat))+
  scale_color_manual(values = c("darkgreen", "green", "lightgreen"))+
  theme_classic()

p4 <- df %>%
  mutate(mass_cat = case_when(mass == 2 ~ "small", 
                              mass == 5 ~ "medium", 
                              mass == 10 ~ "large")) %>%
  ggplot(aes(x = temp-273.15, y = metab/mass))+
  geom_line(aes(color = mass_cat))+
  scale_color_manual(values = c("darkgreen", "green", "lightgreen"))+
  theme_classic()


cowplot::plot_grid(p3, p4, nrow = 1)




df <- expand.grid(temp = c(19+273.15), mass = seq(2,10, length.out = 100))

df$I_normal <- 1 * df$mass^0.75
df$I_warm <- 1*df$mass^0.5
df$metabolism <- 0.5 *df$mass^0.75


df %>%
  pivot_longer(cols = I_normal:metabolism) %>%
  ggplot(aes(x = mass, y = value/mass))+
    geom_line(aes(color = name), linewidth = 2)+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values = c("blue", "red", "gray50"), labels = c("Ingestion (Normal)", "Ingestion (Warm)", "Standard metabolic demand"))+
  labs(x = "Mass", y = "Energy gains or requirements", color = "")+
  theme_classic()+
  theme(legend.position = c(0.8,0.9))







