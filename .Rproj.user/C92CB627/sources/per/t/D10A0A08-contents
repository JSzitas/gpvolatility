


# previous parameters is some matrix of N samples of k parameters 
# previous density is just the vector of previous importance weights
# lambda is the shrinkage parameter
rapf_point_estimate <- function( previous_parameters, previous_density, lambda = 0.5 ) {
  
  # n_particles <- nrow(previous_parameters)
  # n_params <- ncol(previous_parameters)
  
  parameter_means <- colMeans( exp(previous_density) * previous_parameters )
  # return(parameter_means)
  # shrink parameters towards their means 
  lambda * previous_parameters + (1-lambda) * parameter_means
}


