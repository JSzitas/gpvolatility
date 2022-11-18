
rapcf <- function( y,
                   n_state_lags = 1,
                   n_observed_lags = 1,
                   n_covar_pars = 3,
                   mean_fun = mean_rapcf,
                   covariance_fun = covariance_rapcf,
                   n_particles = 100,
                   lambda = 0.95,
                   burn_in = max( ceiling(0.1* length(y)), 10 )#,
                   # add initial parameters if not random
                   # observation_likelihood_fun
                   ) {

  max_lags = max( n_observed_lags, n_state_lags )
  # burn_in -
  # initialize states
  states <- matrix( 0, nrow = length(y), ncol = n_particles )
  states[seq_len(burn_in+max_lags),] <- log( rnorm( (burn_in + max_lags) * n_particles,
                                                    mean = 0,
                                                    sd = stats::sd(y[seq_len(burn_in)]) )^2
  )
  # initialize particle hyper-parameters
  # fix hyperparameters for now - hence the 5
  # since this is easy to overlook, we store steps as the first array dimension
  # and indexing into it yields a matrix where particles are rows and parameters are
  # columns

  n_pars = n_state_lags + n_observed_lags + n_covar_pars
  particle_parameters = array( NA, dim = c( length(y), n_particles, n_pars ) )
  particle_parameters[seq_len(burn_in+max_lags+1),
                      seq_len(n_particles),] = runif( (burn_in+max_lags+1) * n_particles * n_pars, 0,5)

  weights <- matrix( 1/n_particles, nrow = length(y) - 1, ncol = n_particles )
  pred_likelihood <- matrix(NA, nrow = length(y), ncol = n_particles)
    # particle filter loop
  plot.ts(y[(burn_in+1):length(y)])


    for( step in (burn_in+1):(length(y) - 1) ) {
      cat(step, '\n')

      previous_density = weights[ step - 1,]
      # get previous parameters from particle array
      previous_parameters = particle_parameters[step - 1,,]

      # shrink towards the mean using a discounting factor lambda - kind of resembles ad-hoc
      # exponential smoothing, but if it works, it aint broken :)
      parameter_estimates = rapf_point_estimate( previous_parameters,
                                                 previous_density,
                                                 lambda = lambda )
      # turn states and x into a particle data structure so its easier to pass along?
      particles = particle_data( matrix( states[ 1:(step-1), ],
                                         ncol = n_particles),
                                 # observed volatility
                                 y[ 1:(step-1) ],
                                 state_lags = n_state_lags,
                                 observed_lags = n_observed_lags,
                                 parameter_estimates )
      # turn states and x into a particle data structure so its easier to pass along?
      # samples of state from posterior gp for each particle
      states_samples = gp_predict_particles( particles, mean_fun, covariance_fun )
      # pdfs on expected
      #   # here we need to get model densities - ie. using either mnormt::dmnorm for normal or
      #   # mnormt::dmt for multivariate T dist
      m_density_mu = rapcf_density_fun( y[step], states_samples, transform = exp )
      # if( any(m_density_mu < 0) ) {
      #   cat("Negative densities encountered - returning for debug.")
      #   return(list( y[step], states_samples ))
      # }
      # return(list(m_density_mu, previous_density))
      # g = m_density_mu + previous_density
      g = m_density_mu + previous_density
      # resample using importance weights
      cat(paste0("Unique densities: ", length(unique(g)), "\n"))
      # if(length(unique(g))< 2) {
      #   cat("Total particle sampler degradation - returning for debug.")
      #   return( list( densities = weights[(burn_in+1):step,],
      #                 states = states_samples, particles = particles, g = g ) )
      # }
      # return(g)
      # first normalize the weights
      new_weights = normalize( g )
      # cat( new_weights, "\n" )
      # get indices for new weights
      # return(new_weights)
      weight_indices = resample_alright( new_weights )
      # return(list(new_weights, weight_indices))
      cat(paste0("Unique particles sampled: ", length(unique(weight_indices)), "\n"))
      # weight_indices = resample_direct( new_weights )
      # cat(weight_indices, "\n")
      # Propage importance sampled particles
      states[ 1:(step-1), ] = states[ 1:(step-1), weight_indices ]
      # find sd of hyperparameters for jitter
      parameter_sd <- apply( previous_parameters, 2, sd)
      # generate new hyperparameters
      new_parameters <- mnormt::rmnorm( n = n_particles,
                                        varcov = (1-(lambda^2)) * diag( parameter_sd )) +
        parameter_estimates[weight_indices]

      # propose new log variances
      particles = particle_data( matrix( states[ 1:(step-1), ],
                                         ncol = n_particles),
                                 # same as a bit higher
                                 y[ 1:(step-1) ],
                                 state_lags = n_state_lags,
                                 observed_lags = n_observed_lags,
                                 new_parameters )
      new_states = gp_predict_particles( particles, mean_fun, covariance_fun )

      # calculate observation probability on new states
      m_density = rapcf_density_fun( y[step], new_states, transform = exp )
      # adjust for proposal density
      # this sometimes leads to negative weights - I think this is a bug
      # the original code has -, but the paper says / - and this holds for
      # logs, but not for transformed densities (imo) - so perhaps?
      weights[ step,] = m_density - m_density_mu[weight_indices]
      # return(weights[step,])
      # weights[step,] = normalize(weights[step,])
      # cat(weights[step,],"\n")
      # store parameters
      particle_parameters[step,,] = new_parameters
      # store states
      states[step,] = new_states
      # plot current steps
      plot.ts(y[(burn_in):length(y)])
      lines(apply(states[(burn_in):step,], 1, mean)+exp(1), col = "red" )
      # # one step forward prediction
      # particles = particle_data( matrix( states[ 1:step, ],
      #                                    ncol = n_particles),
      #                            # observed volatility
      #                            y[ 1:step ],
      #                            state_lags = n_state_lags,
      #                            observed_lags = n_observed_lags,
      #                            new_parameters )
      # predicted = gp_predict_particles( particles, mean_fun, covariance_fun )
      # # predictive lkelihood
      # pred_likelihood[step,] = rapcf_density_fun( y[step+1], predicted, transform = exp )
    }
  # )
  return(list( states = states[(burn_in+1):nrow(states),],
               weights = weights[(burn_in+1):nrow(weights),],#[nrow(weights),],
               particle_parameters = particle_parameters[(burn_in+1):nrow(particle_parameters),,],
               predictive_likelihood = pred_likelihood[(burn_in+1):nrow(pred_likelihood),],
               y = y[(burn_in+1):length(y)]))
}


#   % predict 1- step forward
#   predParts = gpVolPredFun(parts(1:t,:),data(1:t,:),gpStruct,newParams,isT);
#   % prediction of t+1 at t
#   predLikAll(t,:)=modelStruct.obsfunc(data(t+1,:),predParts,newParams,isT,modelStruct.transfunc);
