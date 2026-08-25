#' @title Model simulation on input data for SEIR model
#' @description
#' Simulate outputs of model for plotting purposes
#'
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0
#' @param simulation list of simulation arguments, as returned by \link[SMCfeatures]{model_seir}
#' @param args a list of arguments as defined in \link[SMCfeatures]{SMCfeatures}
#' @return list like object representing simulation outputs
#' t_data: time steps for data (args$input_data)
#' y_data: output of simulation on input data as dataframe of SEIR
#' @export

model_for_plots_seir <- function (parameters, simulation, args) {
  simulation$t_plot <- seq(0, 500, by = 1)
  simulation$y_plot <- args$simulate_model(parameters, simulation$t_plot)
  return(simulation)
}
