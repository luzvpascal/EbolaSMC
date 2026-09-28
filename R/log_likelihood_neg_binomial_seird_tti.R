#' @title Log-Likelihood function for the logistic growth case study
#' @description
#' Negative likelihood function for SEIR model, with limited data (e.g. only I)
#' @param parameters vector of parameters, with last parameters the standard deviation for each species
#' @param simulation list like object that contains output_data (output of simulation)
#' @param args a list like object of arguments that includes a boolean include_data_constraints, that indicates if data constraints should be included, and output_data (output of observed data)
#' @return log likelihood of simulation outputs against empirical data if args$include_data_constraints set to TRUE, 0 otherwise.
#' @export
#'
log_likelihood_neg_binomial_seird_tti <- function (parameters, simulation, args) {
  n_data <- nrow(args$output_data)
  names_obs <- names(args$output_data) #delta_Du, delta_ItotD, delta_Dtot

  all_names_states <- names_states("SEIRD")
  names <- all_names_states$names_tti
  all_names <- all_names_states$all_names_tti

  names_pars <- names_parameters("SEIRD")$names_tti

  parms_rep_rate <- parameters[length(names_pars)+#4 parameters
                                 length(names)+ #initial conditions
                                 seq(1,2*length(names_obs),2) #reporting rate
  ]

  parms_dispersion <- (parameters[length(names_pars)+#4 parameters
                                    length(names)+ #initial conditions
                                    seq(2,2*length(names_obs),2) #dispersion parameter
  ])


  #compute delta_Itot and delta_Dtot for the simulation
  diff_time <- diff(args$input_data)
  y_data <- data.frame(delta_Du=pmax(1e-6,diff(simulation$y_data$Du)/diff_time))
  y_data$delta_Idtot <- pmax(1e-6,diff(simulation$y_data$Idtot)/diff_time)
  y_data$delta_Dtot <- pmax(1e-6,diff(simulation$y_data$Dtot)/diff_time)

  if (args$include_data_constraints) {
    log_likelihood <- c()
    names_negbin <- c("delta_Du", "delta_Idtot","delta_Dtot")
    selected_sp_negbin <- which(names_negbin %in% names_obs)
    for (id in seq_along(selected_sp_negbin)){
      sp <- selected_sp_negbin[id]
      indexes_not_NA <- which(!is.na(args$output_data$neg_bin[,id]))

      log_likelihood[id] <- sum(dnbinom(round(args$output_data$neg_bin[indexes_not_NA,id]),
                                                size=parms_dispersion[id],
                                                mu=parms_rep_rate[id]*y_data[indexes_not_NA,sp],
                                                log = TRUE))
    }
    return(sum(log_likelihood))
  } else {
    return(0)
  }
}
