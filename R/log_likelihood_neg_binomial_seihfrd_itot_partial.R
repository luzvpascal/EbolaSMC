#' @title Log-Likelihood function for the logistic growth case study
#' @description
#' Negative likelihood function for SEIR model, with limited data (e.g. only I)
#' @param parameters vector of parameters, with last parameters the standard deviation for each species
#' @param simulation list like object that contains output_data (output of simulation)
#' @param args a list like object of arguments that includes a boolean include_data_constraints, that indicates if data constraints should be included, and output_data (output of observed data)
#' @return log likelihood of simulation outputs against empirical data if args$include_data_constraints set to TRUE, 0 otherwise.
#' @export
#'
log_likelihood_neg_binomial_seihfrd_itot_partial <- function (parameters, simulation, args) {
  names_obs_neg_bin <- names(args$output_data$neg_bin) #delta_Itot and delta_Dtot
  names_obs_gaussian <- names(args$output_data$gaussian) #H RH

  names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D")
  all_names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D", "Itot", "Dinst")

  # Negative binomial
  parms_rep_rate <- parameters[11+#11 parameters
                                 length(names)+ #initial conditions
                                 seq(1,2*length(names_obs_neg_bin),2) #reporting rate delta_H, delta_RH, delta_Itot and delta_Dtot
  ]

  parms_dispersion <- (parameters[11+#11 parameters
                                    length(names)+ #initial conditions
                                    seq(2,2*length(names_obs_neg_bin),2) #dispersion parameter delta_H, delta_RH, delta_Itot and delta_Dtot
  ])


  #compute delta_Itot and delta_Dtot for the simulation
  diff_time <- diff(args$input_data)
  y_data <- data.frame(delta_H=pmax(1e-6,diff(simulation$y_data$H)/diff_time))
  y_data$delta_RH <- pmax(1e-6,diff(simulation$y_data$RH)/diff_time)
  y_data$delta_Itot <- pmax(1e-6,diff(simulation$y_data$Itot)/diff_time)
  y_data$delta_Dinst <- pmax(1e-6,diff(simulation$y_data$Dinst)/diff_time)

  # Gaussian
  sigma_log_likelihood <- parameters[seq_along(names_obs_gaussian)+#each species has its own measurement noise
                                       11+#9 parameters
                                       2*length(names_obs_neg_bin)+#4 reporting rate and dispersion parameters
                                       length(names)] #initial conditions

  if (args$include_data_constraints) {
    log_likelihood_neg_bin <- c()
    names_negbin <- c("delta_H","delta_RH","delta_Itot","delta_Dinst")
    selected_sp_negbin <- which(names_negbin %in% names_obs_neg_bin)
    for (id in seq_along(selected_sp_negbin)){
      sp <- selected_sp_negbin[id]
      log_likelihood_neg_bin[id] <- sum(dnbinom(round(args$output_data$neg_bin[,id]),
                                                size=parms_dispersion[id],
                                                mu=parms_rep_rate[id]*y_data[,sp],
                                                log = TRUE))
    }

    log_likelihood_gaussian <- c()
    selected_sp_gaussian <- which(all_names %in% names_obs_gaussian)
    for (id in seq_along(selected_sp_gaussian)){
      sp <- selected_sp_gaussian[id]
      log_likelihood_gaussian[id] <- sum(dnorm(simulation$y_data[,sp],
                                               mean=args$output_data$gaussian[,id],
                                               sd = sigma_log_likelihood[id],
                                               log=TRUE))
    }
    return(sum(log_likelihood_neg_bin)+sum(log_likelihood_gaussian))
  } else {
    return(0)
  }
}

# log_likelihood_neg_binomial_seihfrd_itot_partial <- function (parameters, simulation, args) {
#   names_obs <- names(args$output_data) #delta_H, delta_RH, delta_Itot and delta_Dtot
#
#   names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D")
#   all_names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D", "Itot", "Dinst")
#
#   # Negative binomial
#   parms_rep_rate <- parameters[11+#11 parameters
#                                  length(names)+ #initial conditions
#                                  seq(1,2*length(names_obs),2) #reporting rate delta_H, delta_RH, delta_Itot and delta_Dtot
#   ]
#
#   parms_dispersion <- (parameters[11+ #11 parameters
#                                     length(names)+ #initial conditions
#                                     seq(2,2*length(names_obs),2) #dispersion parameter delta_H, delta_RH, delta_Itot and delta_Dtot
#   ])
#
#
#   #compute delta_Itot and delta_Dtot for the simulation
#   diff_time <- diff(args$input_data)
#   y_data <- data.frame(delta_H=pmax(1e-6,diff(simulation$y_data$H)/diff_time))
#   y_data$delta_RH <- pmax(1e-6,diff(simulation$y_data$RH)/diff_time)
#   y_data$delta_Itot <- pmax(1e-6,diff(simulation$y_data$Itot)/diff_time)
#   y_data$delta_Dinst <- pmax(1e-6,diff(simulation$y_data$Dinst)/diff_time)
#
#   if (args$include_data_constraints) {
#     log_likelihood_neg_bin <- c()
#     names_negbin <-  c("delta_H","delta_RH","delta_Itot","delta_Dinst")
#     selected_sp_negbin <- which(names_negbin %in% names_obs)
#     for (id in seq_along(selected_sp_negbin)){
#       sp <- selected_sp_negbin[id]
#       log_likelihood_neg_bin[id] <- sum(dnbinom(round(args$output_data[,id]),
#                                                 size=parms_dispersion[id],
#                                                 mu=parms_rep_rate[id]*y_data[,sp],
#                                                 log = TRUE))
#     }
#     return(sum(log_likelihood_neg_bin))
#   } else {
#     return(0)
#   }
# }

