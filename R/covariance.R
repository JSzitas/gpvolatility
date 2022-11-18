
single_mat_squared_abs_distance <- function( X ) {
  D <- matrix( 0, nrow = nrow(X), ncol = nrow(X) )
  for( i in seq_len(nrow(X)) ) {
    for(j in seq_len(nrow(X))) {
      for( k in seq_len(ncol(X)) ) {
        D[i,j] <- D[i,j] + (X[i,k] - X[j,k])^2
      }
    }
  }
  return(D)
}

two_mat_squared_abs_distance <- function( X1, X2 ) {
    # unsigned int i,j,k;
  D <- matrix( 0, nrow = nrow(X1), ncol = nrow(X2) )

    # /* for each row of X1 and X2 */
      for( i in seq_len(nrow(X1))) {
        for( j in seq_len(nrow(X2)) ) {
          # /* sum the squared entries */
            for( k in seq_len(ncol(X1)) ) {
              D[i,j] <- D[i,j] + (X1[i,k] - X2[j,k])^2
            }
        }
      }
  return(D)
}

squared_abs_distance <- function( X1, X2 = NULL ) {
  if(is.null(X2)) {
    return( single_mat_squared_abs_distance(X1) )
  }
  return( two_mat_squared_abs_distance(X1, X2) )
}

# X is a matrix
# length_scale and gamma are just random BS scaling parameters
# covariance_rapcf <- function( X1, X2= NULL, parameters = list( length_scale = 10, gamma = 0.5 ) ) {
#   length_scale = parameters$length_scale
#   gamma = parameters$gamma
#   gamma * exp( -1/2 * plgp::distance(X1, X2)/length_scale^2  )
# }

covariance_rapcf <- function( X1, X2= NULL, parameters = list( length_scale = 10, gamma = 0.5 ) ) {
  length_scale = parameters$length_scale
  gamma = parameters$gamma
  gamma * exp( -1/2 * plgp::distance(X1, X2)/length_scale^2)
                 #squared_abs_distance(X1, X2)/length_scale^2  )
}

mean_rapcf <- function( X, parameters ) {
  c(X %*% unlist(parameters))
}

# x_test is a single test point used to generate the prediction -
# for gpVol it is the last row of x_test
gp <- function( y,
                x,
                x_test,
                sigma_n = 0.5,
                mean_fun = mean_rapcf,
                mean_pars = list( a = 1, b = 1 ),
                covariance_fun = covariance_rapcf,
                covariance_pars = list( length_scale = 10, gamma = 0.5 ) ) {
  # in case x_test was cast to vector from a matrix, this should fix it
  x_test <- matrix(x_test, ncol = ncol(x))
  # compute kernel, m and m_star (means for the gaussian process both for training
  # and test data )
  K <- covariance_fun(x, NULL, covariance_pars)

  m <- mean_fun(x, mean_pars)
  m_star <- mean_fun( x_test, mean_pars )
  # precompute K(X, X_star) - which is the same as transpose of K( X_star, X)
  K_x_star <- covariance_fun( x, x_test, covariance_pars )
  K_star_x <- t(K_x_star)
  # precompute K( X_star, X_star)
  K_star_star <- covariance_fun( x_test, x_test, covariance_pars )

  # precompute the inverse of the kernel matrix so we dont have to recompute it
  # many times - here is where the sigma_n enters in
  K_inv <- solve( K + diag( ncol(K) ) * sigma_n^2 )

  # mean equation becomes
  # mu = m + K( sigma2 * I + K-1 )-1 (y - m)
  # see Rasmussen, Williams, chapter 2, eq. 2.38
  # http://www.gaussianprocess.org/gpml/chapters/RW2.pdf
  mu_star <- m_star + K_star_x %*% K_inv %*% ( y - m )
  # sigma is unchanged by inclusion of explicit mean
  sigma_star <- K_star_star - (K_star_x %*% K_inv) %*% K_x_star
  return( c(mu = mu_star, sigma = sigma_star) )
}
