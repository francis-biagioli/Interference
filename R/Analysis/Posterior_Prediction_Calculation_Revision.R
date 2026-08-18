###########################################################################################
### Extract Posterior Predictions of C and R, Matrix Calculations for wC and ahR
###########################################################################################

################################################# load libraries ################################################# 

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load Posteriod distributions of pred and prey abundances

Abundance_prey_full <- as.matrix(read.csv("~/Mass-Abundance Scaling/Abundance_prey_full.csv"))

Abundance_pred_full <- as.matrix(read.csv("~/Mass-Abundance Scaling/Abundance_pred_full.csv"))

### Load in Posterior distributions of parameter estimates from FoRAGE re-fits

w_posterior <- as.matrix(read.csv("~/posterior_w_matrix_Revision.csv"))

h_posterior <- as.matrix(read.csv("~/posterior_h_matrix_Revision.csv"))

a_posterior <- as.matrix(read.csv("~/posterior_scaled_a_matrix_Revision.csv"))

### load FoRAGE

forage <- read.csv('~/Analysis/forage_modified_Revision.csv', stringsAsFactors = TRUE)

#Remove old C/R and wC/ahR columns IF NEEDED
#forage <- forage[,-c(28:39)]

################################################# Calculate w/ah posterior distribution ################################################# 

### Calculate ah

ah_posterior <- a_posterior * h_posterior

## Calculate w/ah

w_ah_posterior <- w_posterior / ah_posterior

w_ah_posterior <- cbind(Inter_ID = seq_len(nrow(w_ah_posterior)),
                        as.data.frame(w_ah_posterior))

write.csv(w_ah_posterior, "~/Analysis/w_ah_post.csv")

################################################# Calculate C/R Matrix ################################################# 

C_R_predict <- Abundance_pred_full / Abundance_prey_full

################################################# Calculate wC/ahR posterior distribution ################################################# 

### Sample posterior parameter draws to make the number of iterations equal to that of the abundance estimate posteriors

#draws <- sample(seq_len(30000), 4000, replace = FALSE)

#a_post_sub <- a_posterior[, draws]
#h_post_sub <- h_posterior[, draws]
#w_post_sub <- w_posterior[, draws]

### Remove rows corresponding to studies omitted in the analysis (14, 25, 38)

#w_post_sub <- w_post_sub[-c(14, 15, 38), ]
w_post_sub <- w_posterior[-c(14, 15, 38), ]

#a_post_sub <- a_post_sub[-c(14, 15, 38), ]
a_post_sub <- a_posterior[-c(14, 15, 38), ]

#h_post_sub <- h_post_sub[-c(14, 15, 38), ]
h_post_sub <- h_posterior[-c(14, 15, 38), ]


### Calculate posteriors of 

wC_predict <- (Abundance_pred_full - (1/1e4)) * w_post_sub

# Set wC values equal to or less than 0 to NA 
wC_predict[wC_predict <= 0] <- NA

ahR_predict <- Abundance_prey_full * (a_post_sub * h_post_sub)

wC_ahR_predict <- wC_predict/ahR_predict

################################################# Calculate Credible Intervals of Each Matrix ################################################# 

# C/R Credible intervals

C_R_ci <- t(apply(C_R_predict, 1, quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE))

colnames(C_R_ci) <- c("C_R_5","C_R_25", "C_R_50", "C_R_75", "C_R_95")

# wC/ahR Credible intervals

wC_ahR_ci <- t(apply(wC_ahR_predict, 1, quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE))

colnames(wC_ahR_ci) <- c("wC_ahR_5", "wC_ahR_25", "wC_ahR_50", "wC_ahR_75", "wC_ahR_95")

# Bind CIs to forage dataset

forage <- cbind(forage, C_R_ci, wC_ahR_ci)

################################################# Save FoRAGE file with CI's of Posteriors ################################################# 

write.csv(wC_ahR_predict, file = '~/Analysis/wC_ahR_predict.csv', row.names = FALSE)

write.csv(forage, file = '~/Analysis/forage_modified_Revision.csv', row.names = FALSE)


