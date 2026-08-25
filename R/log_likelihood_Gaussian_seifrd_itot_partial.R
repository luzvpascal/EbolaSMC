#' @title Log-Likelihood function for the logistic growth case study
#' @description
#' Gaussian likelihood function for SEIR model, with limited data (e.g. only I)
#' @param parameters vector of parameters, with last parameters the standard deviation for each species
#' @param simulation list like object that contains output_data (output of simulation)
#' @param args a list like object of arguments that includes a boolean include_data_constraints, that indicates if data constraints should be included, and output_data (output of observed data)
#' @return log likelihood of simulation outputs against empirical data if args$include_data_constraints set to TRUE, 0 otherwise.
#' @export
#'
log_likelihood_Gaussian_seifrd_itot_partial <- function (parameters, simulation, args) {
  n_data <- nrow(args$output_data)
  names_obs <- names(args$output_data)
  names <- c("S","E", "I", "Fu", "R", "D")
  all_names <- c("S","E", "I", "Fu", "R", "D","Itot","Dinst")
  sigma_log_likelihood <- parameters[seq_along(all_names)+
                                       6+
                                       length(names)] #each species has its own measurement noise

  if (args$include_data_constraints) {
    log_likelihood <- c()
    selected_sp <- which(all_names %in% names_obs)
    for (id in seq_along(selected_sp)){
      sp <- selected_sp[id]
      log_likelihood[id] <- -n_data * log(sigma_log_likelihood[sp]) - sum(sum((simulation$y_data[,sp] -
                                                                                 args$output_data[,id])^2/(2 * sigma_log_likelihood[sp]^2)))
    }
    return(sum(log_likelihood))
  } else {
    return(0)
  }
}
