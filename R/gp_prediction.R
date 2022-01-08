
# gpVolPredFun(parts(1:t-1,:),data(1:t-1,:),gpStruct,paramEst,isT);

# this only takes one set of parameters - ie only one particle 
gp_predict_fun <- function( y, x, x_test, mean_fun, covariance_fun, parameters ) {
  # y = latent states
  # x = observed states 
  
  # covariance_pars = parameters$covariance_pars
  # mean_pars = parameters$mean_pars
  # sigma_n = parameters$sigma_n
  # hardcode for now :/ 
  # a, b (the mean parameters)
  mean_pars = parameters[ 1:2 ]
  mean_pars = as.list(mean_pars)
  names(mean_pars) <- c("a","b")
  # the lengthscale and gamma 
  covariance_pars = parameters[ 3:4 ]
  covariance_pars = as.list(covariance_pars)
  names(covariance_pars) <- c("length_scale","gamma")
  # the sigma_n 
  sigma_n = parameters[ 5 ]
  # gp prediction
  pars <- gp( y,
              x, 
              x_test,
              sigma_n = sigma_n,
              mean_fun = mean_fun,
              mean_pars = mean_pars,
              covariance_fun = covariance_fun,
              covariance_pars = covariance_pars )
  # sample state from the posterior
  rnorm( 1, mean = pars[1], sd = pars[2] )
}

gp_predict_particles <- function( particles,
                                  mean_fun,
                                  covariance_fun) {
  # particles is a list whose elements are particle "properties" - ie the state, 
  # the covariates, and the corresponding parameters 
  purrr::map_dbl( particles, function(particle) {
    y = particle[["y"]]
    x = particle[["x"]]
    x_test = particle[["x_test"]]
    pars = particle[["parameters"]]
    gp_predict_fun( y,
                    x,
                    x_test,
                    mean_fun = mean_fun,
                    covariance_fun = covariance_fun,
                    pars )
  })
}



