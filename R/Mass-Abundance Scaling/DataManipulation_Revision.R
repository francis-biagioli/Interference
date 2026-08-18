###########################################################################################
### Cleaning the data to predict functional response parameters
###########################################################################################
#CODE ADAPTED FROM DR. KYLE COBLENTZ 


### load libraries

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load bayesian regressions for mass-abundance relationship

load('~/Mass-Abundance Scaling/MassAbundanceScaling.RData')

### need to connect the mass-abundance and metabolic scaling regressions to FoRAGE to get

### load FoRAGE

forage <- read.csv('~/Analysis/FoRAGE_Inter_Revision.csv', stringsAsFactors = TRUE)

### drop rows to not include

#forage <- forage %>% filter(Include. == 1)

### drop levels of Prey Major groups that no longer exist

forage$Prey_group_1 <- droplevels(forage$Prey_group_1)

### drop studies without prey masses

forage <- forage %>% filter(Prey_mass > 0)

### make sure all of the predators have masses

which(is.na(forage$Pred_mass))

# drop studies that do not have predator masses

forage <- forage %>% filter(!is.na(Pred_mass))

### add abundance group column to match studies with the appropriate
### regression for the prey in that study

forage <- forage %>% mutate(PreyAbundanceGroup = ifelse(Prey_Vert_Invert== 'Protozoan' | Prey_Vert_Invert == 'Algae', 'Protist',
                                                        ifelse(Prey_group_1 == 'Mammal', 'Mammal',
                                                               ifelse(Prey_group_1 == 'Bird', 'Bird',
                                                                      ifelse(Prey_Vert_Invert == 'Prokaryote', 'Prokaryote',
                                                                             ifelse(Prey_group_1 == 'Fish'| Prey_group_1 == 'Amphibian', 'ecto_vertebrate', 'invertebrate'))))))  
### add metabolism group column for the predators in each study

forage <- forage %>% mutate(PredAbundanceGroup = ifelse(Pred_Vert_Invert == 'Protozoan' | Pred_Vert_Invert == 'Unicell', 'Protist',
                                                        ifelse(Pred_group_1 == 'Mammal', 'Mammal', 
                                                        ifelse(Pred_group_1 == 'Bird', 'Bird', 
                                                               ifelse(Pred_group_1 == 'Fish' | Pred_group_1 == 'Amphibian' | Pred_group_1 == 'Reptile', 'ecto-vertebrate', 'invertebrate')))))

### add columns for masses of prey and predators in grams

forage$PreyMass_g <- forage$Prey_mass/1000

forage$PredMass_g <- forage$Pred_mass/1000


#Declare num studies and iterations

n_studies <- 40
n_draws <- 30000

### run a loop that uses the appropriate mass-abundance regression to predict
### prey densities given the prey masses

Abundance_prey_full <- matrix(NA, nrow = n_studies, ncol = n_draws)

for (i in 1:n_studies) {
  
  # choose correct mass-abundance model for this prey
  fit_to_use <- switch(forage$PreyAbundanceGroup[i],
                       "Protist"        = protist_fit,
                       "Mammal"         = mammal_fit,
                       "Bird"           = bird_fit,
                       "ecto_vertebrate"= ecto_vertebrate_fit,
                       "invertebrate"   = invertebrate_fit,
                       prokaryote_fit)   # default
  
  # new data must contain Mass_g (model does the log-transform internally)
  newdat <- data.frame(Mass_g = forage$PreyMass_g[i])
  
  # generate 10,000 posterior predictive draws
  post_pred <- posterior_predict(
    object  = fit_to_use,
    newdata = newdat,
    draws   = n_draws
  )
  
  # posterior_predict returns a 1 × draws matrix → convert to vector
  Abundance_prey_full[i, ] <- as.numeric(post_pred)
  
  print(i)
}

### predictions are on the log scale so take the exponential to get back to
### the natural scale

Abundance_prey_full <- exp(Abundance_prey_full)

### convert abundance to a data frame and then collate this data frame 
### with FoRAGE

PreyAbundance_median <- apply(Abundance_prey_full, 1, median)

#colnames(PreyAbundance_median) <- c('PreyAbundance_50')

forage <- forage %>%
  mutate(
    PreyAbundance_50   = PreyAbundance_median
  )

### now need to add predator abundances rates using a similar loop

Abundance_pred_full <- matrix(NA, nrow = n_studies, ncol = n_draws)

for (i in 1:n_studies) {
  
  # choose correct mass-abundance model for this prey
  fit_to_use <- switch(forage$PredAbundanceGroup[i],
                       "Protist"        = protist_fit,
                       "Mammal"         = mammal_fit,
                       "Bird"           = bird_fit,
                       "ecto_vertebrate"= ecto_vertebrate_fit,
                       "invertebrate"   = invertebrate_fit,
                       prokaryote_fit)   # default
  
  # new data must contain Mass_g (model does the log-transform internally)
  newdat <- data.frame(Mass_g = forage$PredMass_g[i])
  
  # generate 4,000 posterior predictive draws
  post_pred <- posterior_predict(
    object  = fit_to_use,
    newdata = newdat,
    draws   = n_draws
  )
  
  # posterior_predict returns a 1 × draws matrix → convert to vector
  Abundance_pred_full[i, ] <- as.numeric(post_pred)
  
  print(i)
}

### predictions are on the log scale so take the exponential to get back to
### the natural scale

Abundance_pred_full <- exp(Abundance_pred_full)

### convert abundance to a data frame and then collate this data frame 
### with FoRAGE

PredAbundance_median <- apply(Abundance_pred_full, 1, median)

forage <- forage %>%
  mutate(
    PredAbundance_50   = PredAbundance_median
  )

# ### save this data set as a .csv

write.csv(Abundance_prey_full, "~/Mass-Abundance Scaling/Abundance_prey_full.csv", row.names = FALSE)

write.csv(Abundance_pred_full, "~/Mass-Abundance Scaling/Abundance_pred_full.csv", row.names = FALSE)

write.csv(forage, file = '~/Analysis/forage_modified_Revision.csv', row.names = FALSE)



