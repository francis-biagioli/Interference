#####################################################################################
### Bayesian reanalysis of Hatton et al. 2019. PNAS Mass-Abundance Scaling
#####################################################################################
### CODE WRITTEN BY DR. KYLE COBLENTZ

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load data

scale_data <- read.csv('~/Mass Abundance Scaling/AbundanceScaling.csv')

### plot mass-abundance relationship

ggplot(data = scale_data, aes(x = log(Mass_g), y = log(Density_Nm2), color = Major_taxa)) + geom_point()

### want to estimate the mass-abundance models done by Hatton et al. (but only 
### the ones that we will need for FoRAGE)

### will just be log-log regressions

### Mammals

mammal <- scale_data %>% filter(Major_taxa == 'Mammal')

mammal_fit <- brm(formula = log(Density_Nm2) ~ 0 + Intercept + log(Mass_g), data = mammal, chains = 3, iter = 10000)

### summary

summary(mammal_fit, prob = 0.9)

### Protists

protist <- scale_data %>% filter(Major_taxa == 'Protist')

protist_fit <- brm(formula = log(Density_Nm2) ~ 0 + Intercept + log(Mass_g), data = protist, chains = 3, iter = 10000)

### summary

summary(protist_fit, prob = 0.9)

### Ectotherms -- invertebrate

invertebrate <- scale_data %>% filter(Major_taxa == 'Invertebrate')

invertebrate_fit <- brm(formula = log(Density_Nm2) ~ 0 + Intercept + log(Mass_g), data = invertebrate, chains = 3, iter = 10000)

### summary

summary(invertebrate_fit, prob = 0.9)

### Ectotherms -- vertebrates

ecto_vertebrate <- scale_data %>% filter(Major_taxa == 'Reptilia' | Major_taxa == 'Amphibia' | Major_taxa == 'Fish')

ecto_vertebrate_fit <- brm(formula = log(Density_Nm2) ~ 0 + Intercept + log(Mass_g), data = ecto_vertebrate, chains = 3, iter = 10000)

### summary

summary(ecto_vertebrate_fit, prob = 0.9)

### birds

bird <- scale_data %>% filter(Major_taxa == 'Bird')   

bird_fit <- brm(formula = log(Density_Nm2) ~ 0 + Intercept + log(Mass_g), data = bird, chains = 3, iter = 10000)

### summary

summary(bird_fit, prob = 0.9)

### bacteria -- note that this is really mostly marine phytoplankton <20um from one study

prokaryote <- scale_data %>% filter(Major_taxa == 'Prokaryote')

prokaryote_fit <- brm(formula = log(Density_Nm2) ~ 0 + Intercept + log(Mass_g), data = prokaryote, chains = 3, iter = 10000)

summary(prokaryote_fit, prob = 0.9)

save.image(file = 'MassAbundanceScaling.RData')






