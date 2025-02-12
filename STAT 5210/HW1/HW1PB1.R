
# Given values

y <- c(16, 21, 18) # study values (y1, y2, y3)
pr=c(0.5,0.25,0.25) # Sample probability, P(Ai)
pi <- c(0.75, 0.75, 0.50) # first-order inclusion probabilities (pi1, pi2, pi3)
pi_ij <- matrix(c(NA, 0.5, 0.25,  # pairwise inclusion probabilities
                  0.5, NA, 0.25,
                  0.25, 0.25, NA), 
                nrow=3, byrow=TRUE)

# Define samples and their probabilities

samples <- list(c(1,2), c(1,3), c(2,3))
sample_probs <- c(0.5, 0.25, 0.25)

# Initialize results dataframe

results <- data.frame(Sample = c("{1,2}", "{1,3}", "{2,3}"),
                      HT_Estimator = NA,
                      HT_Variance = NA,
                      SYG_Variance = NA)

# Function to compute HT estimator

compute_HT <- function(sample) {
  return(sum(y[sample] / pi[sample]))
}

# Function to compute HT variance estimator

compute_HT_variance <- function(sample) {
  i <- sample[1]
  j <- sample[2]
  
  # Pairwise term (times 2 for ordering)
  pairwise_term <- 2 * ((pi_ij[i,j] - pi[i] * pi[j]) / pi_ij[i,j]) * (y[i] / pi[i]) * (y[j] / pi[j])
  
  # Single terms
  single_term <- ((1 - pi[i]) / pi[i]^2) * y[i]^2 + ((1 - pi[j]) / pi[j]^2) * y[j]^2
  
  return(pairwise_term + single_term)
}

# Function to compute SYG variance estimator
compute_SYG_variance <- function(sample) {
  i <- sample[1]
  j <- sample[2]
  
  return( (-1) * ((pi_ij[i,j] - pi[i] * pi[j]) / pi_ij[i,j]) * ((y[i] / pi[i]) - (y[j] / pi[j]))^2)
} # Here, single term is zero so only cross term is enough.


# Compute values for each sample

for (k in 1:length(samples)) {
  results$HT_Estimator[k] <- compute_HT(samples[[k]])
  results$HT_Variance[k] <- compute_HT_variance(samples[[k]])
  results$SYG_Variance[k] <- compute_SYG_variance(samples[[k]])
}

results

# Compute true variance of HT estimator
HT_estimator_mean=sum(pr*results$HT_Estimator)
variance_HT=sum((results$HT_Estimator-HT_estimator_mean)^2*pr)
variance_HT

# Check Unbiasedness
round(sum(pr*results$HT_Variance),3)==round(variance_HT,3)
round(sum(pr*results$SYG_Variance),3)==round(variance_HT,3)

