#' @title Differential model for SEIRV
#' @description
#' Differential model for SEIRV
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIRV projection over time giving by input
#' @export

model_seir <- function(parameters, input){

  beta <- parameters[1]
  sigma <- parameters[2]
  gamma <- parameters[3]

  # Initial state
  params <- c(beta=beta,
              sigma=sigma,
              gamma=gamma)

  initial_state <- parameters[c(4,5,6,7)]
  names(initial_state) <- c("S","E", "I", "R")
  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = SMCfeatures::diff_model_seir,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
