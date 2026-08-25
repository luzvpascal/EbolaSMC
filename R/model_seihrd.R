#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta, betaH, betaFu, sigma, phi, delta, delta2, gammaH, gammaD, gammaI, gammaIH, gammaDH,gammaFu, S0, E0, I0, H0, F0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_H0,sigma_likelihood_F0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihrd <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaH <- parameters[2] #transmission rate from funerals

  sigma <- parameters[3] # 1/incubation time

  phi <- parameters[4] #probability of going to the hospital

  gammaH <- parameters[5] #rate of hospital admissions

  #parameters with D compartement
  delta <- parameters[6] #case fatality ratio
  deltaH <- parameters[7] #case fatality ratio in hospital

  gammaDI <- parameters[8] #death rate if not in hospital
  gammaRI <- gammaDI #recovery rate if survivor

  gammaDH <- gammaDI #death rate if in hospital
  gammaRH <- gammaDI #recovery rate if in hospital

  # Initial state
  params <- c(beta=beta,
              betaH=betaH,

              sigma=sigma,

              phi=phi,

              gammaH=gammaH,

              delta=delta,
              gammaDI=gammaDI,
              gammaRI=gammaRI,

              deltaH=deltaH,
              gammaDH=gammaDH,
              gammaRH=gammaRH
              )

  all_names <- c("S","E", "I", "H", "R", "RH", "D", "DH")
  initial_state <- parameters[seq_along(all_names)+8]
  names(initial_state) <- all_names
  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seihrd,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
