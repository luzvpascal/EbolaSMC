#' @title Differential model for SEIR
#' @description
#' Differential model for SEIR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihr_identifiable <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaH <- parameters[2] #transmission rate from funerals
  sigma <- parameters[3] # 1/incubation time
  phi_gammaH <- parameters[4] #case fatality ratio*death rate if infected
  one_minus_phi_gammaI <- parameters[5] #
  gammaIH <- parameters[6] #1/time in hospital

  # Initial state
  params <- c(beta=beta,
              betaH=betaH,
              sigma=sigma,
              phi_gammaH=phi_gammaH,
              one_minus_phi_gammaI=one_minus_phi_gammaI,
              gammaIH=gammaIH)

  initial_state <- parameters[seq(7,11)]
  names(initial_state) <- c("S","E", "I", "H", "R")
  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = SMCfeatures::diff_model_seihr_identifiable,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
