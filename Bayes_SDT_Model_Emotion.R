# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Bayesian signal detection model.
#
# v .90
#
# Andrew L. Cohen
# 6/3/16
#
# Notes:
#   Assumes all rating values used at least once.
#   Assumes equal repetitions across stimuli.
#   To fix the scale across different induced moods:
#     The central criterion for all participants is fixed at 0.
#     The first tested mood group (alphabetical) is assumed to have a SD of 1. 
#   The group graphs use the independent means and SDs for each group. 
#     As the mean and SD are dependent, this is only an approximation.
#
#
# Use of this software comes with ABSOLUTELY NO WARRANTY. 
# Double check everything!
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Startup =================================

# Remove objects --------------
rm(list=ls()) 

# Load required libraries ---------------
require(rjags)
require(coda)
require(runjags)
require(parallel)

# Constants =================================

# Data --------------

# Name of the data file
# File format - each line: 
#   participant induced_mood tested_mood_group stimulus rating
#     participant = subject identifier (integer or factor)
#     induced_mood = induced mood condition (integer or factor)
#     tested_mood_group = mood being queried (integer or factor)
#     stimulus = stimulus within the tested_mood_group (integer or factor)
#     rating = participant rating (integer, higher = more endorsed)
#   Example:
#     1 Anger fear_group fearful 4	
#   The variables names must match these.
#data_file_name <- 'Data/empref_manip_bsdt.csv' 
#data_file_name <- 'Data/text_BayesSDT.csv'
data_file_name <- 'Data/faces_BayesSDT.csv'

# Participants to remove before analysis (NULL=remove none)
# List of (induced_mood, participant)
# Example: list(list('Anger', 1), list('Fear', 5))
remove_participants <- NULL

# Options ------------

# Allow participant x group interactions
pg_interactions <- TRUE
#pg_interactions <- FALSE

# Induced mood conditions to analyze (NULL=analyze all)
# Array of induced moods to analyze
# Example: c('Anger', 'Fear', 'Neutral')
analyze_induced_moods <- NULL
#analyze_induced_moods <- c('calm', 'fear', 'happy')

# Tested mood conditions to analyze (NULL=analyze all)
# Array of tested moods to analyze
# Example: c('angry_group', 'fear_group')
analyze_tested_moods <- NULL
#analyze_tested_moods <- c('angry_group', 'fear_group', 'calm_group')

# Load previously generated mcmc samples
load_samples <- FALSE
load_sample_file_name <- 'SDT_model_mcmc_samples.R'
#load_sample_file_name <- 'SDT_model_mcmc_samples_Andrea_text.R'

# Save mcmc samples
save_samples <- TRUE
#save_sample_file_name <- 'SDT_model_mcmc_samples.R'
#save_sample_file_name <- 'SDT_model_mcmc_samples_Andrea_text.R'
save_sample_file_name <- 'SDT_model_mcmc_samples_Andrea_faces.R'
if (load_samples) save_samples <- FALSE # Can't load and save at same time

# Name of the BUGS code
BUGS_file_name <- 'Bayes_SDT_model.bug'

# Name of results output file
output_file_name <- 'Bayes_SDT_output.txt'

# Which analytics to show ------------

# Global on/off
show_mcmc_analytics <- FALSE

# Individual analytics
print_sample_summary <- FALSE
plot_samples <- FALSE
print_gelman_diag <- FALSE
plot_gelman <- FALSE
autocorr_diag <- FALSE
autocorr_plot <- FALSE
effective_size <- FALSE

# Plots of chains, to check for mixing
#plot_chains <- TRUE
plot_chains <- FALSE

# MCMC constants ------------

n_burnin <- 1000 # Burnin steps
n_adapt <- 1000 # Adaptation steps
n_total_samples <- 10000 # Total number of samples
thin_steps <- 10 # Thinning steps
n_chains <- max(c(1, detectCores() - 2)) # Chains (# processors - 2)
n_samples <- ceiling(n_total_samples/n_chains) # Samples per chain

# Model parameters  ------------

# Shape and rate for SD variability.
# Not a direct parameter in the model.
# mode=1, sd=3
sd_uninformed_shape <- 1.393 
sd_uninformed_rate <- .393

# Determines variability of deltas at mood level 
d_shape_shape <- 1.588 # mode=1.4, sd=3
d_shape_rate <- .42
d_rate_shape <- 1.142 # mode=.4, sd=3
d_rate_rate <- .356

# Determines variability of mood group sigma
g_sigma_lambda_shape <- sd_uninformed_shape 
g_sigma_lambda_rate <- sd_uninformed_rate

p_mu <- 0 # Center of participant effects
p_lambda_shape <- sd_uninformed_shape # Determines SD of center of participant effects
p_lambda_rate <- sd_uninformed_rate

g_mu_min <- -10 # Min/max tested mood group mean
g_mu_max <- 10 

# Only used if pg_interactions is TRUE
if(pg_interactions) {
  pg_mu <- 0 # Center of participant x group effects
  pg_lambda_shape <- sd_uninformed_shape # Determines SD of center of participant x group effects
  pg_lambda_rate <- sd_uninformed_rate
}

s_mu <- 0 # Center of stimulus effects
s_lambda_shape <- sd_uninformed_shape # Determines SD of center of stimulus effects
s_lambda_rate <- sd_uninformed_rate

# Plotting constants ------------------

# HDI 
HDI_prop <- .95 # HDI proportion (e.g., .95 = 95%)
HDI_mult <- c(.95, .90) # For staggering height of criterion error bars

# Da contrasts
# Defines which Da contrasts to perform
#   group_1 = A single tested mood group for each induced mood
#   group_2 = A set (1 or more) of tested mood groups to compare against group_1 for each induced mood
# Default (NULL) is to compare the top two tested mood groups.
# Example: 
# da_contrasts <- list(group_1=c('angry_group', 'fear_group', 'calm_group'), 
#                      group_2=list(c('fear_group'), c('angry_group', 'calm_group'), c('fear_group')))
da_contrasts <- NULL

# Plot buffers
sd_buffer <- 2.5 # Determines x-axis for distribution plots
y_max_buffer <- .1 # How much to add (proportion) to y-axis in distribution plots

# Colors/lines
distr_lwds <- 2 # One for all or one for each tested mood group
distr_ltys <- 1 # One for all or one for each tested mood group
distr_cols <- NULL # Null cycles 2, 3, etc or one for each tested mood group
da_col <- 'gray'

# Save the distribution plots
# Doesn't save the sample plots
save_plots <- FALSE
plot_file_name <- 'Bayes_SDT_plot.png'

# Get data ===========================

# Read raw data
d <- read.csv(data_file_name, header=TRUE, sep=',')

# Remove unwanted participants
if (!is.null(remove_participants)) {
  for (part in remove_participants) {
    d <- d[!(d$induced_mood==part[[1]] & d$participant==part[[2]]), ]
  }
}

# Remove unwanted induced moods
if (!is.null(analyze_induced_moods)) {
  d <- d[is.element(d$induced_mood, analyze_induced_moods), ]
}

# Remove unwanted tested moods
if (!is.null(analyze_tested_moods)) {
  d <- d[is.element(d$tested_mood_group, analyze_tested_moods), ]
}

# Store level names (for plotting)
induced_mood_names <- levels(droplevels(d$induced_mood))
tested_mood_group_names <- levels(droplevels(d$tested_mood_group)) 

# Recode data
# Start all counts from 1
# Necessary for use with table & to get correct counts
d_order <- d[order(d$induced_mood, d$participant, d$tested_mood_group, d$stimulus),]

induced_mood_new <- as.numeric(d_order$induced_mood)
participant_new <- numeric()
participant_old <- list()
k <- 1
for (m in unique(d_order$induced_mood)) {
  participant_new <- c(participant_new, 
                       as.numeric(as.factor(droplevels(d_order[d_order$induced_mood==m, 'participant']))))
  participant_old[[k]] <- as.character(unique(d_order[d_order$induced_mood==m, 'participant']))
  k <- k + 1
}
tested_mood_group_new <- as.numeric(d_order$tested_mood_group)
stimulus_new <- numeric()
for (g in unique(d_order$tested_mood_group)) {
  stimulus_new <- c(stimulus_new,
                    as.numeric(as.factor(as.character(d_order[d_order$tested_mood_group==g, 'stimulus']))))
}
rating_new <- d_order$rating

# Store recoded data
d_new <- data.frame(induced_mood=induced_mood_new,
                    participant=participant_new,
                    stimulus=stimulus_new,
                    tested_mood_group=tested_mood_group_new,
                    rating=rating_new,
                    stringsAsFactors=FALSE)

# Put into the right format (array/table) for analysis
# Assumes all rating values used at least once
ratings <- with(d_new, table(induced_mood,participant,tested_mood_group,stimulus,rating))

# Counts
n_m <- length(unique(d$induced_mood))
d_tmp <- d[!duplicated(d[, c('induced_mood', 'participant')]), c('induced_mood', 'participant')]
n_p <- as.numeric(table(d_tmp$induced_mood))
n_g <- length(unique(d$tested_mood_group))
d_tmp <- d[!duplicated(d[, c('tested_mood_group', 'stimulus')]), c('tested_mood_group', 'stimulus')]
n_s <- as.numeric(table(d_tmp$tested_mood_group))
n_s <- n_s[n_s > 0]
n_r <- length(unique(d$rating))
n_s_reps <- max(rowSums(ratings[1,,1,1,])) # Assumes equal number of repetitions across stimuli

# Store the stimuli names
stimulus_names <- matrix(NA, nrow=n_g, ncol=max(n_s))
k <- 1
for (g in unique(d_order$tested_mood_group)) {
  stimulus_names[k,] <- as.vector(unique(droplevels(d[d$tested_mood_group==g, 'stimulus'])))
  k <- k + 1
}

# The Bayesian model ===============================

if (!load_samples) {
  
  # The data list ----------------------
  
  # Passed into runjags
  data_list <- list(ratings=ratings, # Data
                    n_m=n_m, # Data constants
                    n_p=n_p,
                    n_g=n_g,
                    n_s=n_s,
                    n_r=n_r,
                    n_s_reps=n_s_reps,
                    d_shape_shape=d_shape_shape, # Prior constants
                    d_shape_rate=d_shape_rate,
                    d_rate_shape=d_rate_shape,
                    d_rate_rate=d_rate_rate,
                    g_sigma_lambda_shape=g_sigma_lambda_shape,
                    g_sigma_lambda_rate=g_sigma_lambda_rate,
                    p_mu=p_mu,
                    p_lambda_shape=p_lambda_shape,
                    p_lambda_rate=p_lambda_rate,
                    g_mu_min=g_mu_min,
                    g_mu_max=g_mu_max,
                    s_mu=s_mu,
                    s_lambda_shape=s_lambda_shape,
                    s_lambda_rate=s_lambda_rate)
  
  # Allow pxg interactions
  if (pg_interactions) {
    data_list <- c(data_list, 
                   pg_mu=pg_mu,
                   pg_lambda_shape=pg_lambda_shape,
                   pg_lambda_rate=pg_lambda_rate)
  }
  
}

# The model ----------------------

if (!load_samples) {
  
  model_str_opening <- "
model {

# Indices throughout:
#   m = induced Mood
#   p = Participant
#   g = tested mood Group
#   s = Stimulus
#   r = Rating
"
  model_str_ratings <- "
# Ratings distributed as a multinomial
# ratings: n_m x n_p x n_g x n_s x n_r
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {
    for (g in 1:n_g) {
      for (s in 1:n_s[g]) {

        ratings[m,p,g,s,] ~ dmulti(probs[m,p,g,s,1:n_r], n_s_reps)

      }
    }
  }
}
"
model_str_probs <- "
# Multinomial probabilities
# probs: m_n x n_p x n_g x n_s x n_r
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {
    for (g in 1:n_g) {
      for (s in 1:n_s[g]) {
        
        # First bin
        probs[m,p,g,s,1] <- phi((k[m,p,1]-mus[m,p,g,s])/g_sigmas[m,g])

        # Last bin
        probs[m,p,g,s,n_r] <- 1-phi((k[m,p,n_r-1]-mus[m,p,g,s])/g_sigmas[m,g])

        # Other bins
        for (r in 2:(n_r-1)) {
    
          probs[m,p,g,s,r] <- phi((k[m,p,r]-mus[m,p,g,s])/g_sigmas[m,g]) - 
                              phi((k[m,p,r-1]-mus[m,p,g,s])/g_sigmas[m,g])

        }

      }
    }
  }
}
"

if (n_r%%2==0) {
  # Even ratings
  
  model_str_k <- "
# Criteria for each participant 
# Assumes even number of rating categories
# k: n_m x n_p x n_r-1
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {

    k[m,p,round(n_r/2)] <- 0 # Middle criterion is always 0 for every participant

    for (r in 1:round(n_r/2-1)) {
      k[m,p,r] <- k[m,p,round(n_r/2)] - sum(deltas[m,p,1:round(n_r/2-r)])
    }
  
    for (r in round(n_r/2+1):(n_r-1)) {
      k[m,p,r] <- k[m,p,round(n_r/2)] + sum(deltas[m,p,round(n_r/2):(r-1)])
    }
      
  }
}" } else {
  # Odd ratings
  
  model_str_k <- "
# Criteria for each participant 
# Assumes odd number of rating categories
# k: n_m x n_p x n_r-1
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {
  
    # Centered on 0
    # Middle response region defined by delta 1 (1/2 on each side)
    k[m,p,(n_r-1)/2] <- -deltas[m,p,1]/2
    k[m,p,(n_r-1)/2+1] <- deltas[m,p,1]/2

    for (r in 1:(((n_r-1)/2)-1)) {
      k[m,p,r] <- -deltas[m,p,1]/2 - sum(deltas[m,p,2:((n_r-1)/2-r+1)])
    }

    for (r in (((n_r-1)/2)+2):(n_r-1)) {
      k[m,p,r] <- deltas[m,p,1]/2 + sum(deltas[m,p,(((n_r-1)/2)+1):(r-1)])
    }
  
  }
}
"
}

model_str_deltas <- "
# Delta for each participant
# Deltas define the distances between criteria
# deltas: n_m x n_p x n_r-2
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {
    for (r in 1:(n_r-2)) {

      deltas[m,p,r] ~ dgamma(d_shape[m,r], d_rate[m,r])

    }
  }
}
"

model_str_d_shape_rate <- "
# Shape and rate for deltas
# shape, rate: n_m x n_r
for (m in 1:n_m) {
  for (r in 1:(n_r-2)) {

    d_shape[m,r] ~ dgamma(d_shape_shape, d_shape_rate)
    d_rate[m,r] ~ dgamma(d_rate_shape, d_rate_rate)

  }
}
" 

if (!pg_interactions) {
  # Don't allow pxg interactions
  
  model_str_mus <- "
# Overall means 
# No mood effects because coufounded with participant effects
#   That is, assume different participants in each mood group
# mus: n_m x n_p x n_g x n_s
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {
    for (g in 1:n_g) {
      for (s in 1:n_s[g]) {

        mus[m,p,g,s] <- p_effects[m,p] + g_effects[m,g] + s_effects[g,s]

      }
    }
  }
}
" } else {
  # Allow pxg interactions
  
  model_str_mus <- "
# Overall means 
# No mood effects because coufounded with participant effects
#   That is, assume different participants in each mood group
# mus: n_m x n_p x n_g x n_s
for (m in 1:n_m) {
  for (p in 1:max(n_p)) {
    for (g in 1:n_g) {
      for (s in 1:n_s[g]) {

        mus[m,p,g,s] <- p_effects[m,p] + g_effects[m,g] + pg_effects[m,p,g] + s_effects[g,s]

      }

      # Monitor for later participant x group da tests
      mpg_mus[m,p,g] <- p_effects[m,p] + g_effects[m,g] + pg_effects[m,p,g]
    }
  }
}
"  
}

model_str_p_effects <- "
# Partipant effects 
# p_effects: n_m x n_p
for (m in 1:n_m) {
  for (p in 1:max(n_p)) { ###### CHECK THIS!!!!!!!!!!!!!!!!!!!!!! why not n_p[m]

    p_effects[m,p] ~ dnorm(p_mu, p_lambda) 

  }
}

# SD for participant effects
p_lambda ~ dgamma(p_lambda_shape, p_lambda_rate)
p_sigma <- 1/sqrt(p_lambda)
"

model_str_g_effects <- "
# Test mood group effects
# Varies by induced mood
# Uniform b/c assume mood effects are not normally distributed
# g_effects: n_m x n_g
for (m in 1:n_m) {
  for (g in 1:n_g) {
    
    g_effects[m,g] ~ dunif(g_mu_min, g_mu_max)

  }
}
"

if (!pg_interactions) {
  # Don't allow pxg interactions
  
  model_str_pg_effects <- ""
  
} else {
  # Allow pxg interactions
  
  model_str_pg_effects <- "
# Partipant x group interaction effects 
# pg_effects: n_m x n_p x n_g
for (m in 1:n_m) {
  for (p in 1:max(n_p)) { ###### CHECK THIS!!!!!!!!!!!!!!!!!!!!!! why not n_p[m]
    for (g in 1:n_g) {

      pg_effects[m,p,g] ~ dnorm(pg_mu, pg_lambda) 

    }
  }
}
 
# SD for participant effects
pg_lambda ~ dgamma(pg_lambda_shape, pg_lambda_rate) 
pg_sigma <- 1/sqrt(pg_lambda) 
"  
}

model_str_s_effects <- "
# Stimulus effects
# Assumed independent of induced mood 
# s_effects: n_g x n_s
for (g in 1:n_g) {
  for (s in 1:n_s[g]) {
  
    s_effects[g,s] ~ dnorm(s_mu, s_lambda)

  }
}

# SD for stimulus effects
s_lambda ~ dgamma(s_lambda_shape, s_lambda_rate) 
s_sigma <- 1/sqrt(s_lambda)
"

model_str_g_lambdas_sigmas <- "
# SD of moods overall
# Assumes same sd for all participants for a given mood
#   Model was unstable with individual participant mood SDs
for (m in 1:n_m) {
  g_sigmas[m,1] <- 1 # Fix for scaling
  for (g in 2:n_g) {

    g_lambdas[m,g] ~ dgamma(g_sigma_lambda_shape, g_sigma_lambda_rate)
    g_sigmas[m,g] <- 1/sqrt(g_lambdas[m,g])

  }
}
"

model_str_closing <- "
}
"

# Write the BUGS model to a file
model_str <- paste(model_str_opening, model_str_ratings, model_str_probs, model_str_k, model_str_deltas, model_str_d_shape_rate, model_str_mus, model_str_p_effects, model_str_g_effects, model_str_pg_effects, model_str_s_effects, model_str_g_lambdas_sigmas, model_str_closing)
writeLines(model_str, con=BUGS_file_name)

}

# Parameters to watch -----------------

if (!load_samples) {
  parameters = c("p_effects", "g_effects", "s_effects", "g_sigmas", "p_sigma", "s_sigma", "d_shape", "d_rate")
  if (pg_interactions) parameters <- c(parameters, "mpg_mus")
}

# Run the model  -----------------

if (!load_samples) {
  
  # Turn off some warnings
  runjags.options(inits.warning=FALSE, rng.warning=FALSE)
  
  # The runjags model
  jags.model <- run.jags(method='parallel',
                         model=BUGS_file_name,
                         monitor=parameters,
                         data=data_list,
                         n.chains=n_chains,
                         adapt=n_adapt,
                         burnin=n_burnin,
                         sample=n_samples,
                         thin=thin_steps,
                         summarise=FALSE,
                         plots=FALSE)
  samps <- as.mcmc.list(jags.model)
  
  # Save the samples, if requested
  if (save_samples) {
    date <- date()
    save(samps, 
         date,
         data_file_name, 
         pg_interactions,
         remove_participants, analyze_induced_moods, analyze_tested_moods,
         n_m, n_p, n_g, n_s, n_r, n_s_reps, 
         n_burnin, n_adapt, n_total_samples, thin_steps, n_chains, n_samples,
         sd_uninformed_shape, sd_uninformed_rate, d_shape_shape, d_shape_rate, d_rate_shape, d_rate_rate, g_sigma_lambda_shape, g_sigma_lambda_rate, p_mu, p_lambda_shape, p_lambda_rate, pg_mu, pg_lambda_shape, pg_lambda_rate,
         g_mu_min, g_mu_max, s_mu, s_lambda_shape, s_lambda_rate,
         HDI_prop, da_contrasts,
         file=save_sample_file_name)
  }
  
} else {
  
  # Load the samples, if requested
  load(load_sample_file_name)
  
}

# MCMC analytics ==========================

if (show_mcmc_analytics) {
  
  if (print_sample_summary) print(summary(samps))
  if (plot_samples) plot(samps)
  if (print_gelman_diag) print(gelman.diag(samps))
  if (plot_gelman) gelman.plot(samps)
  if (autocorr_diag) print(autocorr.plot(samps))
  if (autocorr_plot) autocorr.plot(samps)
  if (effective_size) effectiveSize(samps)
  
}

# Analyze and plot output ==========================

# Put samples into matrix form
samps_matrix <- as.matrix(samps)
n_total_samps <- nrow(samps_matrix)

# Plotting constants ------------------

g_lwd <- numeric(n_g)
g_lty <- numeric(n_g)
g_cols <- numeric(n_g)
g_names <- character(n_g)
for (g in 1:n_g) {
  
  if (length(distr_lwds)==1) g_lwd[g] <- distr_lwds
  else g_lwd[g] <- distr_lwds[g]
  
  if (length(distr_ltys)==1) g_lty[g] <- distr_ltys
  else g_lty[g] <- distr_ltys[g]
  
  if (is.null(distr_cols)) g_cols[g] <- g+1
  else g_cols[g] <- distr_cols[g]
  
  g_names[g] <- tested_mood_group_names[g]
  
}

# Graphs -----------------------

if (save_plots) png(filename=plot_file_name, width=15, height=10, units='in', res=300)

# Induced mood graphs and information ----------

# General graphing parameters
n_fig_rows <- length(induced_mood_names)
n_fig_cols <- 4
layout_mat <- matrix(c(1:n_fig_rows, 1:n_fig_rows, 
                       (n_fig_rows+1):(2*n_fig_rows),
                       (2*n_fig_rows+1):(3*n_fig_rows),
                       (3*n_fig_rows+1):(4*n_fig_rows)),
                     nrow=n_fig_rows, byrow=FALSE)
layout(layout_mat)

# Distribution graphs ---------------
# NOTE: USING INDEPENDENT MEANS FOR MEANS AND SDS.

# Graphing parameters
par(mar=c(5, 5, 3, 1))
par(pty='m')

# Determine distribution parameters from samples

# Group effects
g_effect_means <- matrix(NA, nrow=n_m, ncol=n_g)
g_effect_HDI <- array(NA, c(n_m, n_g, 2))
g_sigmas <- matrix(NA, nrow=n_m, ncol=n_g)
g_sigmas_HDI <- array(NA, c(n_m, n_g, 2))
x_min <- 0; x_max <- 0; y_max <- 0
for(m in 1:n_m) {
  for (g in 1:n_g) {
    
    g_effect_means[m, g] <- mean(samps_matrix[, paste("g_effects[", m, ",", g, "]", sep='')])
    g_effect_HDI[m, g, ] <- HPDinterval(as.mcmc(samps_matrix[, paste("g_effects[", m, ",", g, "]", sep='')]), prob=HDI_prop)
    g_sigmas[m, g] <- mean(samps_matrix[, paste("g_sigmas[", m, ",", g, "]", sep='')])
    g_sigmas_HDI[m, g, ] <- HPDinterval(as.mcmc(samps_matrix[, paste("g_sigmas[", m, ",", g, "]", sep='')]), prob=HDI_prop)
    
    x_min_tmp <- g_effect_means[m, g] - sd_buffer*g_sigmas[m, g]
    if (x_min_tmp < x_min) x_min <- x_min_tmp
    x_max_tmp <- g_effect_means[m, g] + sd_buffer*g_sigmas[m, g]
    if (x_max_tmp > x_max) x_max <- x_max_tmp
    y_max_tmp <- dnorm(0, 0, g_sigmas[m, g])
    if (y_max_tmp > y_max) y_max <- y_max_tmp
    
  }
  
}

# Criteria
m_criteria_mean <- matrix(NA, nrow=n_m, ncol=n_r-1)
m_criteria_HDI <- array(NA, c(n_m, n_r-1, 2))
delta_samples <- matrix(NA, nrow=dim(samps_matrix)[1], ncol=n_r-2)
for(m in 1:n_m) {
  
  # Delta samples
  # Compute from rate and shape
  for (r in 1:(n_r-2)) {
    delta_samples[,r] <- samps_matrix[, paste("d_shape[", m, ",", r, "]", sep='')] /
      samps_matrix[, paste("d_rate[", m, ",", r, "]", sep='')]
  }
  
  # The criteria
  if (n_r%%2 == 0){
    # Even rating levels
    
    for (r in 1:round(n_r/2-1)) {
      tmp_samples <- -apply(as.matrix(delta_samples[,1:round(n_r/2-r)]),1,sum)
      m_criteria_mean[m,r] <- mean(tmp_samples)
      m_criteria_HDI[m,r,] <- HPDinterval(as.mcmc(tmp_samples), prob=HDI_prop)
    }
    
    for (r in round(n_r/2+1):(n_r-1)) {
      tmp_samples <- apply(as.matrix(delta_samples[,round(n_r/2):(r-1)]),1,sum)
      m_criteria_mean[m,r] <- mean(tmp_samples)
      m_criteria_HDI[m,r,] <- HPDinterval(as.mcmc(tmp_samples), prob=HDI_prop)      
    }
    
    m_criteria_mean[m,round(n_r/2)] <- 0
    m_criteria_HDI[m,round(n_r/2),1:2] <- c(0, 0)
    
  } else {
    # Odd rating levels
    
    r <- (n_r - 1)/2
    tmp_samples <- -delta_samples[,1]
    m_criteria_mean[m,r] <- mean(tmp_samples/2)
    m_criteria_HDI[m,r,] <- HPDinterval(as.mcmc(tmp_samples/2), prob=HDI_prop)
    
    r <- (n_r - 1)/2 + 1
    tmp_samples <- delta_samples[,1]
    m_criteria_mean[m,r] <- mean(tmp_samples/2)
    m_criteria_HDI[m,r,] <- HPDinterval(as.mcmc(tmp_samples/2), prob=HDI_prop)
    
    for (r in 1:((n_r-1)/2-1)) {
      tmp_samples <- apply(as.matrix(-delta_samples[,2:((n_r-1)/2-r+1)]),1,sum)
      m_criteria_mean[m,r] <- mean(tmp_samples)
      m_criteria_HDI[m,r,] <- HPDinterval(as.mcmc(tmp_samples), prob=HDI_prop)      
    }
    
    for (r in ((n_r-1)/2+2):(n_r-1)) {
      tmp_samples <- apply(as.matrix(delta_samples[,(((n_r-1)/2)+1):(r-1)]),1,sum)
      m_criteria_mean[m,r] <- mean(tmp_samples)
      m_criteria_HDI[m,r,] <- HPDinterval(as.mcmc(tmp_samples), prob=HDI_prop)      
    }
    
  }    
  
}

# Draw the distribution graphs
for (m in 1:n_m) {
  
  range <- seq(x_min, x_max, .01)
  y_min <- 0
  y_max <- y_max+y_max_buffer*y_max
  plot(NA, 
       xlim=c(x_min, x_max), ylim=c(y_min, y_max),
       xlab='', ylab='density', main=paste('induced mood: ', induced_mood_names[m], sep=''),
       xaxt='n',
       cex.lab=2, cex.axis=1.5, cex.main=2)
  points(c(0, 0), c(0, y_max), type='l', lwd=2) # Show 0 even when odd number of ratings
  for (g in 1:n_g) {
    points(range, dnorm(range, g_effect_means[m, g], g_sigmas[m, g]), 
           type='l', col=g_cols[g], lwd=g_lwd[g], lty=g_lty[g])
  }
  if (n_r%%2==0) {
    # Even rating levels
    for (r in c(1:(n_r/2-1), (n_r/2+1):(n_r-1))) {
      points(c(m_criteria_mean[m,r], m_criteria_mean[m,r]), c(0, y_max), type='l', lty=2)
      arrows(m_criteria_HDI[m,r,1], y_max*HDI_mult[abs(r%%2-2)], 
             m_criteria_HDI[m,r,2], y_max*HDI_mult[abs(r%%2-2)], 
             code=3, length=.1, angle=90, col='gray')
    }
  } else {  
    # Odd rating levels
    for (r in 1:(n_r-1)) {
      points(c(m_criteria_mean[m,r], m_criteria_mean[m,r]), c(0, y_max), type='l', lty=2)
      arrows(m_criteria_HDI[m,r,1], y_max*HDI_mult[abs(r%%2-2)], 
             m_criteria_HDI[m,r,2], y_max*HDI_mult[abs(r%%2-2)], 
             code=3, length=.1, angle=90, col='gray')
    }
  }
  if (m==1) {
    legend(x='topleft', legend=g_names[g_names!=''], col=g_cols[g_names!=''], lty=1, bty='n', cex=1.5)
  }
}
mtext('mood endorsement strength (dimensionless)', side=1, line=2.5, cex=1.5)

# Da graphs -----------

da_total <- matrix(0, nrow=n_m, ncol=n_total_samps)
da_mean <- numeric(n_m)
da_HDI <- matrix(NA, nrow=n_m, ncol=2)
p_da_total_gt0 <- numeric(n_m)
da_titles <- character(n_m)
for (m in 1:n_m) {
  
  # Da contrasts
  if (!is.null(da_contrasts)) {
    # User defined
    group_1 <- which(tested_mood_group_names==da_contrasts$group_1[m])
    group_2 <- as.numeric(lapply(da_contrasts$group_2[[m]], 
                                 function(x) {which(tested_mood_group_names==x)}))
  } else {
    # Default is to compare the groups with the top two means
    group_1 <- which.max(g_effect_means[m,])
    group_2 <- order(g_effect_means[m,])[length(g_effect_means[m,]) - 1]
  }
  
  da_titles[m] <- paste(g_names[group_1], 'v', paste(g_names[group_2], collapse=', '), sep=' ')
  
  for (g2 in group_2) {
    
    g1_m <- samps_matrix[, paste('g_effects[', m, ',', group_1, ']', sep='')]
    g1_sd <- samps_matrix[, paste('g_sigmas[', m, ',', group_1, ']', sep='')]
    
    g2_m <- samps_matrix[, paste('g_effects[', m, ',', g2, ']', sep='')]
    g2_sd <- samps_matrix[, paste('g_sigmas[', m, ',', g2, ']', sep='')]
    
    da_total[m,] <- da_total[m,] + (g1_m - g2_m)/sqrt((g1_sd^2 + g2_sd^2)/2)
    
  }
  
  da_mean[m] <- mean(da_total[m,])
  da_HDI[m,] <- as.numeric(HPDinterval(as.mcmc(da_total[m,]), prob=HDI_prop))
  
  p_da_total_gt0[m] <- round(sum(da_total[m,]>0)/length(da_total[m,]), 2)
  
}

y_max_hist <- 0
for (m in 1:n_m) {
  h_tmp <- hist(da_total[m,], plot=FALSE)
  if (max(h_tmp$density) > y_max_hist) y_max_hist <- max(h_tmp$density)
}
y_max_hist <- y_max_hist + y_max_hist*y_max_buffer

par(pty='m')
par(mar=c(5, 5, 3, .5))

x_min <- min(da_total)
x_max <- max(da_total)
for (m in 1:n_m) {
  if (m!=n_m) xlab=''  
  if (m==n_m) xlab=expression('d'[a])
  h_tmp <- hist(da_total[m,], 
                main=da_titles[m], xlab=xlab, ylab='density',
                freq=FALSE,
                xlim=c(x_min, x_max),
                ylim=c(0, y_max_hist),
                cex.axis=1.5, cex.lab=2, cex.main=2,
                xaxt='n',
                col=da_col, border='gray')
  if (m!=n_m) axis(side=1, labels=FALSE, cex.axis=1.5, cex.lab=1.5)
  if (m==n_m) axis(side=1, cex.axis=1.5, cex.lab=1.5)
  
  points(c(da_HDI[m,1], da_HDI[m,1]), c(0, y_max_hist), ty='l', lty=2, lwd=2)
  points(c(da_HDI[m,2], da_HDI[m,2]), c(0, y_max_hist), ty='l', lty=2, lwd=2)
  
  p_da_total_gt0_formatted <- formatC(signif(p_da_total_gt0[m], digits=3), 
                                      digits=3, format="fg", flag="#")
  if(da_mean[m] > 0) prob_text <- 'P(da>0)\n'
  else prob_text <- 'P(da<0)\n'
  text(da_mean[m], .075*y_max_hist,
       paste(prob_text, p_da_total_gt0_formatted, sep=''))
  
  points(c(0, 0), c(0, y_max_hist), ty='l', lwd=2)
  
}

# Participant effects
p_effect_means <- matrix(NA, nrow=n_m, ncol=max(n_p))
for (m in 1:n_m) {
  for (p in 1:n_p[m]) {
    
    p_effect_means[m,p] <- mean(samps_matrix[, paste('p_effects[', m, ',', p, ']', sep='')])
    
  }
}

if (pg_interactions) {
  
  da_mp_total <- array(0, c(n_m, max(n_p), n_total_samps))
  da_mp_mean <- matrix(NA, nrow=n_m, ncol=max(n_p))
  da_mp_HDI <- array(NA, c(n_m, max(n_p), 2))
  da_titles <- character(n_m)
  da_mp_in <- list()
  da_mp_out <- list()
  for (m in 1:n_m) {
    
    # Da contrasts
    if (!is.null(da_contrasts)) {
      # User defined
      group_1 <- which(tested_mood_group_names==da_contrasts$group_1[m])
      group_2 <- as.numeric(lapply(da_contrasts$group_2[[m]], 
                                   function(x) {which(tested_mood_group_names==x)}))
    } else {
      # Default is to compare the groups with the top two means
      group_1 <- which.max(g_effect_means[m,])
      group_2 <- order(g_effect_means[m,])[length(g_effect_means[m,]) - 1]
    }
    
    da_titles[m] <- paste(g_names[group_1], 'v', paste(g_names[group_2], collapse=', '), sep=' ')
    
    for (p in 1:n_p[m]) {
      
      for (g2 in group_2) {
        
        g1_m <- samps_matrix[, paste('mpg_mus[', m, ',', p, ',', group_1, ']', sep='')]
        g1_sd <- samps_matrix[, paste('g_sigmas[', m, ',', group_1, ']', sep='')]
        
        g2_m <- samps_matrix[, paste('mpg_mus[', m, ',', p, ',', g2, ']', sep='')]
        g2_sd <- samps_matrix[, paste('g_sigmas[', m, ',', g2, ']', sep='')]
        
        da_mp_total[m,p,] <- da_mp_total[m,p,] + (g1_m - g2_m)/sqrt((g1_sd^2 + g2_sd^2)/2)
        
      }
      
      da_mp_mean[m,p] <- mean(da_mp_total[m,p,])
      da_mp_HDI[m,p,] <- as.numeric(HPDinterval(as.mcmc(da_mp_total[m,p,]), prob=HDI_prop))
      
    }
   
    da_mp_in[[m]] <- as.vector(na.omit(da_mp_HDI[m,,1] > 0))
    da_mp_out[[m]] <- as.vector(na.omit(da_mp_HDI[m,,1] < 0))
     
  }
  
}

par(pty='s')
par(mar=c(5, 4.5, 3, 5.5))
par(new=FALSE)
p_effect_outliers <- list()
for (m in 1:n_m) {
  
  # Participants
  y_min <- min(p_effect_means, na.rm=T); y_min <- y_min*1.1
  y_max <- max(p_effect_means, na.rm=T); y_max <- y_max*1.1
  p_effects_HDI <- as.numeric(HPDinterval(as.mcmc(p_effect_means[1,])))
  p_effects_in_HDI <- na.omit(p_effect_means[m,p_effect_means[m,]>=p_effects_HDI[1] & p_effect_means[m,]<=p_effects_HDI[2]])
  p_effects_out_HDI <- na.omit(p_effect_means[m,p_effect_means[m,]<p_effects_HDI[1] | p_effect_means[m,]>p_effects_HDI[2]])
  p_effect_outliers[[m]] <- which(p_effect_means[m,]<p_effects_HDI[1] | p_effect_means[m,]>p_effects_HDI[2])
  if (pg_interactions) x_pos <- 1/3
  else x_pos <- 1/2
  plot(rep(x_pos, length(p_effects_in_HDI)), 
       p_effects_in_HDI,  
       xlim=c(0, 1), 
       ylim=c(y_min, y_max),
       xlab='', ylab='participant effects',
       cex.lab=2, cex.axis=1.5, cex=1, 
       xaxt='n', col='black')
  points(rep(x_pos, length(p_effects_out_HDI)), 
         p_effects_out_HDI,
         pch=4, col='gray')
  
  points(c(x_pos-.05, x_pos+.05), c(0, 0), ty='l', lwd=2)
  
  # Participant x group das
  if (pg_interactions) {
    
    par(new=TRUE)
    plot(1, 
         ty='n',
         xlim=c(0, 1),
         ylim=c(min(0, min(da_mp_mean, na.rm=TRUE)), max(da_mp_mean, na.rm=TRUE)),
         axes=FALSE, 
         xlab='', ylab='',
         cex.axis=1.5, cex.lab=2)
    
    points(rep(2/3, length(da_mp_mean[m,da_mp_in[[m]]])), 
           da_mp_mean[m,da_mp_in[[m]]],
           pch=1, cex=1, col='black')
    points(rep(2/3, length(da_mp_mean[m,da_mp_out[[m]]])), 
           da_mp_mean[m,da_mp_out[[m]]],
           pch=4, cex=1, col='gray')
    
    points(c(2/3-.05, 2/3+.05), c(0, 0), ty='l', lwd=2)
    
    axis(side=4, cex.axis=1.5)
    mtext('participant da', side=4, line=3.5, cex=1.33) # Not sure why cex isn't 1.5
    
  }
  
}

if (pg_interactions) {
  axis(side=1, at=c(1/3, 2/3), labels=c('effect', 'da'), cex.axis=1.5)
} else {
  axis(side=1, at=c(1/2), labels=c('effect'), cex.axis=1.5)
}
par(new=FALSE)

# Stimulus effects

par(pty='m')
par(mar=c(5, 4.5, 3, 5.5))

s_effect_means <- matrix(NA, nrow=n_g, ncol=max(n_s))
s_effect_HDI <- array(NA, c(n_g, max(n_s), 2))
for (g in 1:n_g) {
  for (s in 1:n_s[g]) {
    tmp_samps <- samps_matrix[, paste('s_effects[', g, ',', s, ']', sep='')]    
    s_effect_means[g, s] <- mean(tmp_samps)
    s_effect_HDI[g, s, ] <- as.numeric(HPDinterval(as.mcmc(tmp_samps)))    
  }
}

for (m in 1:n_m) {
  
  y_min <- min(s_effect_means, na.rm=T); y_min <- y_min - abs(y_min)*.1
  y_max <- max(s_effect_means, na.rm=T); y_max <- y_max + abs(y_max)*.1
  plot(1, 
       type="n", 
       xaxt='n',
       xlab="", ylab="stimulus effects",
       xlim=c(0, 1),
       ylim=c(y_min, y_max),
       cex.lab=2, cex.axis=1.5)
  
  for (g in 1:n_g) {
    for (s in 1:n_s[g]) {
      
      text(1/(n_g+1)*g, na.omit(s_effect_means[g, s]), s, cex=1.5, col=g_cols[g])
      points(c(1/(n_g+1)*g, 1/(n_g+1)*g)-.001, c(min(s_effect_means, na.rm=T), max(s_effect_means, na.rm=T)), ty='l', col=g_cols[g])
      
    }
  }
  
}
mtext('tested mood', side=1, line=2.5, cex=1.5)

if (save_plots) dev.off()

# MCMC graphs
if (plot_chains) {
  
  par(mar=c(2,3,2,1))
  
  # g_effects
  par(pty='s')
  par(mfrow=c(n_m, n_g))
  par(mar=c(3, 4, 4, 1))
  for (m in 1:n_m) {
    for (g in 1:n_g) {
      effect <- paste("g_effects[", m, ",", g, "]", sep='')
      effect_mean <- round(mean(samps_matrix[, effect]), 2)
      y_min <- min(samps_matrix[, effect])
      y_max <- max(samps_matrix[, effect])
      plot(NA, 
           xlim=c(0,n_samples), ylim=c(y_min, y_max), 
           xlab="", ylab='')
      mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
      for (i in 1:n_chains) {
        samp_min <- (i-1)*n_samples + 1
        samp_max <- i*n_samples
        points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
      }
      if (g == 1) mtext(induced_mood_names[m], side=2, line=2.25)
      if (m == 1) mtext(tested_mood_group_names[g], side=3, line=2)
    }
  }
  
  # s_effects
  par(pty='s')
  par(mfcol=c(max(n_s), n_g))
  par(mar=c(3, 4, 3, 1))
  for (g in 1:n_g) {
    for (s in 1:n_s[g]) {
      effect <- paste("s_effects[", g, ",", s, "]", sep='')
      effect_mean <- round(mean(samps_matrix[, effect]), 2)
      y_min <- min(samps_matrix[, effect])
      y_max <- max(samps_matrix[, effect])
      plot(NA, 
           xlim=c(0,n_samples), ylim=c(y_min, y_max), 
           xlab="", ylab='')
      mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
      for (i in 1:n_chains) {
        samp_min <- (i-1)*n_samples + 1
        samp_max <- i*n_samples
        points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
      }
      if (s == 1) mtext(tested_mood_group_names[g], side=3, line=2)
    }
  }
  
  # p_effects
  # 7 random participants in each mood group
  par(pty='s')
  par(mfrow=c(n_m, 7))
  par(mar=c(3, 4, 3, 1))
  for (m in 1:n_m) {
    rand_ps <- sort(sample(1:n_p[m], 7, replace=FALSE))
    for (p in 1:7) {
      effect <- paste("p_effects[", m, ",", rand_ps[p], "]", sep='')
      effect_mean <- round(mean(samps_matrix[, effect]), 2)
      y_min <- min(samps_matrix[, effect])
      y_max <- max(samps_matrix[, effect])
      plot(NA, 
           xlim=c(0,n_samples), ylim=c(y_min, y_max), 
           xlab="", ylab='')
      mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
      for (i in 1:n_chains) {
        samp_min <- (i-1)*n_samples + 1
        samp_max <- i*n_samples
        points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
      }
      if (p == 1) mtext(text=induced_mood_names[m], side=2, line=2)
    }
  }
  
  # sigmas
  # First group is always 1, so skip
  par(pty='s')
  par(mfrow=c(n_m, n_g-1))
  par(mar=c(3, 4, 4, 1))
  for (m in 1:n_m) {
    for (g in 2:n_g) {
      effect <- paste("g_sigmas[", m, ",", g, "]", sep='')
      effect_mean <- round(mean(samps_matrix[, effect]), 2)
      y_min <- min(samps_matrix[, effect])
      y_max <- max(samps_matrix[, effect])
      plot(NA, 
           xlim=c(0,n_samples), ylim=c(y_min, y_max), 
           xlab="", ylab='')
      mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
      for (i in 1:n_chains) {
        samp_min <- (i-1)*n_samples + 1
        samp_max <- i*n_samples
        points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
      }
      if (g == 2) mtext(text=induced_mood_names[m], side=2, line=2.25)
      if (m == 1) mtext(text=tested_mood_group_names[g], side=3, line=2)
    }
  }
  
  # shape
  par(pty='s')
  par(mfrow=c(n_m, n_r-2))
  par(mar=c(3, 4, 4, 1))
  for (m in 1:n_m) {
    for (r in 1:(n_r-2)) {
      effect <- paste("d_shape[", m, ",", r, "]", sep='')
      effect_mean <- round(mean(samps_matrix[, effect]), 2)
      y_min <- min(samps_matrix[, effect])
      y_max <- max(samps_matrix[, effect])
      plot(NA, 
           xlim=c(0,n_samples), ylim=c(y_min, y_max), 
           xlab="", ylab='')
      mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
      for (i in 1:n_chains) {
        samp_min <- (i-1)*n_samples + 1
        samp_max <- i*n_samples
        points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
      }
      if (r == 1) mtext(text=induced_mood_names[m], side=2, line=2)
      if (m == 1) mtext(text=paste('delta ', r), side=3, line=2)
    }
  }
  
  # rate
  par(pty='s')
  par(mfrow=c(n_m, n_r-2))
  par(mar=c(3, 4, 4, 1))
  for (m in 1:n_m) {
    for (r in 1:(n_r-2)) {
      effect <- paste("d_rate[", m, ",", r, "]", sep='')
      effect_mean <- round(mean(samps_matrix[, effect]), 2)
      y_min <- min(samps_matrix[, effect])
      y_max <- max(samps_matrix[, effect])
      plot(NA, 
           xlim=c(0,n_samples), ylim=c(y_min, y_max), 
           xlab="", ylab='')
      mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
      for (i in 1:n_chains) {
        samp_min <- (i-1)*n_samples + 1
        samp_max <- i*n_samples
        points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
      }
      if (r == 1) mtext(text=induced_mood_names[m], side=2, line=2)
      if (m == 1) mtext(text=paste('delta ', r), side=3, line=2)
    }
  }
  
  # other sigmas
  par(pty='s')
  par(mfrow=c(1, 2))
  par(mar=c(3, 4, 4, 1))
  for (i in 1:2) {
    if (i == 1) effect <- 'p_sigma'
    if (i == 2) effect <- 's_sigma'
    effect_mean <- round(mean(samps_matrix[, effect]), 2)
    y_min <- min(samps_matrix[, effect])
    y_max <- max(samps_matrix[, effect])
    plot(NA, 
         xlim=c(0,n_samples), ylim=c(y_min, y_max), 
         xlab="", ylab='')
    mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
    for (i in 1:n_chains) {
      samp_min <- (i-1)*n_samples + 1
      samp_max <- i*n_samples
      points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
    }
  }
  
  # Too many to show all.
  # 1 random participant from each induced mood x group
  if (pg_interactions) {
    
    par(pty='s')
    par(mfrow=c(n_m, n_g))
    par(mar=c(3, 4, 4, 1))
    for (m in 1:n_m) {
      for (g in 1:n_g) {
        rand_ps <- sort(sample(1:n_p[m], 1, replace=FALSE))
        effect <- paste("mpg_mus[", m, ",", rand_ps, ",", g, "]", sep='')
        effect_mean <- round(mean(samps_matrix[, effect]), 2)
        y_min <- min(samps_matrix[, effect])
        y_max <- max(samps_matrix[, effect])
        plot(NA, 
             xlim=c(0,n_samples), ylim=c(y_min, y_max), 
             xlab="", ylab='')
        mtext(text=paste(effect, "=", effect_mean, sep=''), side=3, line=.5)
        for (i in 1:n_chains) {
          samp_min <- (i-1)*n_samples + 1
          samp_max <- i*n_samples
          points(samps_matrix[samp_min:samp_max, effect], ty='l', col=i)
        }
        if (g == 1) mtext(induced_mood_names[m], 2, 2)
        if (m == 1) mtext(tested_mood_group_names[g], 3, 2)
      }
    }
     
  }
  
}

# Text output ------------------------

output_txt <- character()
k <- 1
output_txt[k] <- '==================================='
k <- k+1
output_txt[k] <- 'Data info -------------------------'
k <- k+1
output_txt[k] <- paste('Date model run: ', date())
k <- k+1
output_txt[k] <- paste('data_file_name: ', data_file_name, sep='')
k <- k+1
if(is.null(remove_participants)){
  output_txt[k] <- paste('Removed participants: None')
} else {
  output_txt[k] <- paste('Removed participants: ', paste(unlist(remove_participants), collapse=' '), sep='')
}
k <- k+1
if(pg_interactions) {
  output_txt[k] <- paste('pxg interactions allowed')
} else {
  output_txt[k] <- paste('pxg interactions not allowed')
}
k <- k+1
output_txt[k] <- 'MCMC info -------------------------'
k <- k+1
output_txt[k] <- paste('n_burnin: ', n_burnin, sep='')
k <- k+1
output_txt[k] <- paste('n_adapt: ', n_adapt, sep='')
k <- k+1
output_txt[k] <- paste('n_total_samples: ', n_total_samples, sep='')
k <- k+1
output_txt[k] <- paste('thin_steps: ', thin_steps, sep='')
k <- k+1
output_txt[k] <- paste('n_chains: ', n_chains, sep='')
k <- k+1
output_txt[k] <- paste('n_samples: ', n_samples, sep='')
k <- k+1
output_txt[k] <- 'Model parameter info --------------'
k <- k+1
output_txt[k] <- paste('sd_uninformed_shape: ', sd_uninformed_shape, sep='')
k <- k+1
output_txt[k] <- paste('sd_uninformed_rate: ', sd_uninformed_rate, sep='')
k <- k+1
output_txt[k] <- paste('d_shape_shape: ', d_shape_shape, sep='')
k <- k+1
output_txt[k] <- paste('d_shape_rate: ', d_shape_rate, sep='')
k <- k+1
output_txt[k] <- paste('d_rate_shape: ', d_rate_shape, sep='')
k <- k+1
output_txt[k] <- paste('d_rate_rate: ', d_rate_rate, sep='')
k <- k+1
output_txt[k] <- paste('g_sigma_lambda_shape: ', g_sigma_lambda_shape, sep='')
k <- k+1
output_txt[k] <- paste('g_sigma_lambda_rate: ', g_sigma_lambda_rate, sep='')
k <- k+1
output_txt[k] <- paste('p_mu: ', p_mu, sep='')
k <- k+1
output_txt[k] <- paste('p_lambda_shape: ', p_lambda_shape, sep='')
k <- k+1
output_txt[k] <- paste('p_lambda_rate: ', p_lambda_rate, sep='')
k <- k+1
output_txt[k] <- paste('g_mu_min: ', g_mu_min, sep='')
k <- k+1
output_txt[k] <- paste('g_mu_max: ', g_mu_max, sep='')
k <- k+1
if (pg_interactions) {
  output_txt[k] <- paste('pg_mu: ', pg_mu, sep='')
  k <- k+1
  output_txt[k] <- paste('pg_lambda_shape: ', pg_lambda_shape, sep='')
  k <- k+1
  output_txt[k] <- paste('pg_lambda_rate: ', pg_lambda_rate, sep='')
  k <- k+1
}
output_txt[k] <- paste('s_mu: ', s_mu, sep='')
k <- k+1
output_txt[k] <- paste('s_lambda_shape: ', s_lambda_shape, sep='')
k <- k+1
output_txt[k] <- paste('s_lambda_rate: ', s_lambda_rate, sep='')
k <- k+1
output_txt[k] <- paste('g_mu_max: ', g_mu_max, sep='')
k <- k+1
output_txt[k] <- paste('g_mu_max: ', g_mu_max, sep='')
k <- k+1
output_txt[k] <- paste('g_mu_max: ', g_mu_max, sep='')
k <- k+1
output_txt[k] <- 'Test info -------------------------'
k <- k+1
if (is.null(analyze_induced_moods)) {
  output_txt[k] <- paste('Analyzed induced moods: Default')  
} else {
  output_txt[k] <- paste('Analyzed induced moods: ', paste(unlist(analyze_induced_moods), collapse=' '), sep='')
}
k <- k+1
if (is.null(analyze_tested_moods)) {
  output_txt[k] <- paste('Analyzed tested moods: Default')
} else {
  output_txt[k] <- paste('Analyzed tested moods: ', paste(unlist(analyze_tested_moods), collapse=' '), sep='')
}
k <- k+1
output_txt[k] <- paste('HDI_prop: ', HDI_prop, sep='')
k <- k+1
if (is.null(da_contrasts)) {
  output_txt[k] <- paste('da_contrasts: Default')
} else {
  output_txt[k] <- paste('da_contrasts: ', paste(unlist(da_contrasts), collapse=' '), sep='')
}
k <- k+1
output_txt[k] <- 'Stimulus legend -------------------'
k <- k+1
for (g in 1:n_g) {
  group <- tested_mood_group_names[g]
  stimuli <- paste(1:n_s[g], stimulus_names[g,], collapse=' ')
  output_txt[k] <- paste(group, ': ', stimuli, sep='')
  k <- k + 1
}
output_txt[k] <- 'MCMC summary ----------------------'
k <- k+1
output_txt[k] <- '--Group means & HDIs'
k <- k+1
for (m in 1:n_m) {
  for (g in 1:n_g) {
    output_txt[k] <- paste(induced_mood_names[m], ', ', tested_mood_group_names[g], ': (', g_effect_HDI[m,g,1], ', ', g_effect_means[m,g], ', ', g_effect_HDI[m,g,2], ')', sep='')    
    k <- k+1
  }
}
output_txt[k] <- '--Group SDs & HDIs'
k <- k+1
for (m in 1:n_m) {
  for (g in 1:n_g) {
    output_txt[k] <- paste(induced_mood_names[m], ', ', tested_mood_group_names[g], ': (', g_sigmas_HDI[m,g,1], ', ', g_sigmas[m,g], ', ', g_sigmas_HDI[m,g,2], ')', sep='')    
    k <- k+1
  }
}
output_txt[k] <- '--Criterion means & HDIs'
k <- k+1
for (m in 1:n_m) {
  for (r in 1:(n_r-1)) {
    output_txt[k] <- paste(induced_mood_names[m], ', criterion ', r, ': (', m_criteria_HDI[m,r,1], ', ', m_criteria_mean[m,r], ', ', m_criteria_HDI[m,r,2], ')', sep='')    
    k <- k+1
  }
}
output_txt[k] <- 'Da tests --------------------------'
k <- k+1
output_txt[k] <- '--Da means, HDIs, & p>0'
k <- k+1
for (m in 1:n_m) {
  output_txt[k] <- paste(induced_mood_names[m], ', ', da_titles[m], ': (', da_HDI[m,1], ', ', da_mean[m], ', ', da_HDI[m,2], '), P(da>0)=', p_da_total_gt0[m], sep='')    
  k <- k+1
}
output_txt[k] <- 'Stimulus effects ------------------'
k <- k+1
output_txt[k] <- '--Stimulus effect means, HDIs'
k <- k+1
for (g in 1:n_g) {
  for (s in 1:n_s[g]) {
    output_txt[k] <- paste(tested_mood_group_names[g], ', ', stimulus_names[g,s], ': (', s_effect_HDI[g,s,1], ', ', s_effect_means[g,s], ', ', s_effect_HDI[g,s,2], ')', sep='')    
    k <- k+1
  }
}
output_txt[k] <- 'Outliers --------------------------'
k <- k+1
output_txt[k] <- '--p_effect outliers (Bayes_SDT IDs) (original IDs)'
k <- k+1
# DOUBLE CHECK ORGINAL NUMBERING!!!!!
for (m in 1:n_m) {
  output_txt[k] <- paste(induced_mood_names[m], 
                         ': (', 
                         paste(p_effect_outliers[[m]], collapse=' '), 
                         ') (', 
                         paste(participant_old[[m]][p_effect_outliers[[m]]], collapse=' '), 
                         ')', 
                         sep='')    
  k <- k+1
}
output_txt[k] <- '--participant x group outliers (Bayes_SDT IDs) (original IDs)'
k <- k+1
# DOUBLE CHECK ORGINAL NUMBERING!!!!!
for (m in 1:n_m) {
  da_mp_out_nums <- which(da_mp_out[[m]])
  output_txt[k] <- paste(induced_mood_names[m], 
                         ': (', 
                         paste(da_mp_out_nums, collapse=' '), 
                         ') (', 
                         paste(participant_old[[m]][da_mp_out_nums], collapse=' '), 
                         ')', 
                         sep='')    
  k <- k+1
}
output_txt[k] <- '==================================='

# Show on screen
print(output_txt)

# Save in file
writeLines(output_txt, con=output_file_name)



