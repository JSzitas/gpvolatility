make_particle_data <- function(state,
                               observed,
                               state_lags = 1,
                               observed_lags = 1) {
  # calculate first "valid" (not NA) observation
  first_non_na_index <- max(state_lags, observed_lags) + 1
  state_df <- purrr::map(
    seq_len(state_lags),
    ~ dplyr::lag(state, .x)
  )
  observed_df <- purrr::map(
    seq_len(observed_lags),
    ~ dplyr::lag(observed, .x)
  )
  state_lags <- matrix(do.call(cbind, state_df), ncol = state_lags)
  observed_lags <- matrix(do.call(cbind, observed_df), ncol = observed_lags)

  x <- cbind(state_lags, observed_lags)
  # we need this for casting back to a matrix
  # x_cols = ncol(x)

  total_length <- nrow(x)
  x_test <- x[total_length, , drop = FALSE]
  # use the valid index from earlier
  x <- x[first_non_na_index:(total_length - 1), , drop = FALSE]
  y <- state[first_non_na_index:(total_length - 1)]

  list(y = y, x = x, x_test = x_test)
}

particle_data <- function( x_data,
                           states,
                           state_lags = 1,
                           observed_lags = 1,
                           parameter_estimates) {
  particles <- list()
  for (particle in seq_len(dim(parameter_estimates)[1])) {
    particles[[particle]] <- make_particle_data(x_data[, particle],
      states,
      state_lags = 1,
      observed_lags = 1
    )
    particles[[particle]][["parameters"]] <- parameter_estimates[particle, ]
  }
  return(particles)
}


normalize <- function(x) {

  n <- length(x)
  # this should be an adjustment for numerical stability, apparently
  weights <- exp(x - max(x))
  weights <- weights/sum(weights)
  weights[is.nan(weights)] <- 0
  weights
}

effective_sample_size <- function(x) {
  1 / sum(c(x^2, 1))
}

resample <- function(x) {
  len_x = length(x)
  sample( which( x > (runif(len_x, 0, 1)/len_x) ),
          size = len_x,
          replace = TRUE )
}

resample_alright <- function(x) {
  sample( length(x),
          size = length(x),
          replace = TRUE,
          prob = x
          )
}


resample_direct <- function(x) {
  N = length(x);
  Q = cumsum(x);

  T_ = rep(NA, N)
  T_[1:N] = seq(from = 0, to = 1-(1/N), length.out = N) + runif(1)/N
  T_[N] <- 1
  # T_[N+1] = 1;

  i=1;
  j=1;

  indx <- rep(NA, N)
  while (i<N) {
    tryCatch(
      {
        if (T_[i] < Q[j] ) {
          indx[i] = j;
          i=i+1;
        }
        else {
          j=j+1;
        }
      },
    error = function(e) {
      dbg <<- list( T_, Q, i, j )
      stop()
    })
  }
  return(indx)
}



