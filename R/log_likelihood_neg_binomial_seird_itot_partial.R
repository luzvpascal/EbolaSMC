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
  names_obs <- names(args$output_data) #delta_R, delta_D, delta_Itot

  names <- c("S","E","I","R","D")
  all_names <- c("S","E","I","R","D","Itot")

  parms_rep_rate <- parameters[4+#4 parameters
                                 length(names)+ #initial conditions
                                 seq(1,2*length(names_obs),2) #reporting rate delta_H, delta_RH, delta_Itot and delta_Dtot
  ]

  parms_dispersion <- (parameters[4+#4 parameters
                                    length(names)+ #initial conditions
                                    seq(2,2*length(names_obs),2) #dispersion parameter delta_H, delta_RH, delta_Itot and delta_Dtot
  ])


  #compute delta_Itot and delta_Dtot for the simulation
  diff_time <- diff(args$input_data)
  y_data <- data.frame(delta_R=pmax(1e-6,diff(simulation$y_data$R)/diff_time))
  y_data$delta_D <- pmax(1e-6,diff(simulation$y_data$D)/diff_time)
  y_data$delta_Itot <- pmax(1e-6,diff(simulation$y_data$Itot)/diff_time)

  if (args$include_data_constraints) {
    log_likelihood <- c()
    names_negbin <- c("delta_R", "delta_D", "delta_Itot")
    selected_sp_negbin <- which(names_negbin %in% names_obs)
    for (id in seq_along(selected_sp_negbin)){
      sp <- selected_sp_negbin[id]
      log_likelihood[id] <- sum(dnbinom(args$output_data[,id],
                                    size=parms_dispersion[id],
                                    mu=parms_rep_rate[id]*y_data[,sp],
                                    log = TRUE))
    }
    return(sum(log_likelihood))
  } else {
    return(0)
  }
}
