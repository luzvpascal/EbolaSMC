#' @title Log-Likelihood function for the logistic growth case study
#' @description
#' Negative likelihood function for SEIR model, with limited data (e.g. only I)
#' @param parameters vector of parameters, with last parameters the standard deviation for each species
#' @param simulation list like object that contains output_data (output of simulation)
#' @param args a list like object of arguments that includes a boolean include_data_constraints, that indicates if data constraints should be included, and output_data (output of observed data)
#' @return log likelihood of simulation outputs against empirical data if args$include_data_constraints set to TRUE, 0 otherwise.
#' @export
#'
log_likelihood_neg_binomial_seird_itot_partial <- function (parameters, simulation, args) {
  n_data <- nrow(args$output_data)
  names_obs <- names(args$output_data) #delta_Itot and delta_D

  names <- c("S","E","I","R","D")
  all_names <- c("S","E","I","R","D","Itot")

  parms_rep_rate <- parameters[4+#4 parameters
                                 length(names)+ #initial conditions
                                 c(1,3) #reporting rate delta Itot and delta D
                                 ]
  parms_dispersion <- (parameters[4+#4 parameters
                             length(names)+ #initial conditions
                             c(2,4) #dispersion parameter delta Itot and delta D
                              ])

  #compute delta_Itot and delta_D for the simulation
  diff_time <- diff(args$input_data)
  y_data <- data.frame(delta_Itot=pmax(1e-6,diff(simulation$y_data$Itot)/diff_time))
  y_data$delta_D <- pmax(1e-6,diff(simulation$y_data$D)/diff_time)

  if (args$include_data_constraints) {
    log_likelihood <- c()
    selected_sp <- c("delta_Itot","delta_D")
    for (id in seq_along(selected_sp)){
      log_likelihood[id] <- sum(dnbinom(args$output_data[,id],
                                    size=parms_dispersion[id],
                                    mu=parms_rep_rate[id]*y_data[,id],
                                    log = TRUE))
    }
    return(sum(log_likelihood))
  } else {
    return(0)
  }
}
