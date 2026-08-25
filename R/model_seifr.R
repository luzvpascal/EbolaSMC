#' @title Differential model for SEIR
#' @description
#' Differential model for SEIR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seifr <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaFu <- parameters[2] #transmission rate from funerals
  sigma <- parameters[3] # 1/incubation time
  delta <- parameters[4] #case fatality ratio
  gammaD <- parameters[5] #death rate if infected
  gammaI <- gammaD#parameters[6] #recovery rate if survivor
  gammaFu <- parameters[6] #1/funeral time

  # Initial state
  params <- c(beta=beta,
              betaFu=betaFu,
              sigma=sigma,
              delta=delta,
              gammaD=gammaD,
              gammaI=gammaI,
              gammaFu=gammaFu)

  initial_state <- parameters[seq(7,11)]
  names(initial_state) <- c("S","E", "I", "Fu", "R")
  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = SMCfeatures::diff_model_seifr,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
