#' @title Log-Likelihood function SEIFR
#' @description
#' Gaussian likelihood function for SEIFR model
#' @param parameters vector of parameters, with last parameters the standard deviation for each species
#' @param simulation list like object that contains output_data (output of simulation)
#' @param args a list like object of arguments that includes a boolean include_data_constraints, that indicates if data constraints should be included, and output_data (output of observed data)
#' @return log likelihood of simulation outputs against empirical data if args$include_data_constraints set to TRUE, 0 otherwise.
#' @export
#'
log_likelihood_Gaussian_seifr <- function (parameters, simulation, args) {
  n_data <- nrow(args$output_data)
  sigma_log_likelihood <- parameters[seq(12,16)] #each species has its own measurement noise

  if (args$include_data_constraints) {
    log_likelihood <- c()
    for (sp in seq_along(sigma_log_likelihood)){
      log_likelihood[sp] <- -n_data * log(sigma_log_likelihood[sp]) - sum(sum((simulation$y_data[,sp] -
                                                                                 args$output_data[,sp])^2/(2 * sigma_log_likelihood[sp]^2)))
    }
    return(sum(log_likelihood))
  }
  else {
    return(0)
  }
}
