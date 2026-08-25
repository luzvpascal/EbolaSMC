#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta, betaH, betaFu, sigma, phi, delta, delta2, gammaH, gammaD, gammaI, gammaIH, gammaDH,gammaFu, S0, E0, I0, H0, F0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_H0,sigma_likelihood_F0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihr <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaH <- parameters[2] #transmission rate from funerals

  sigma <- parameters[3] # 1/incubation time

  phi <- parameters[4] #probability of going to the hospital

  gammaH <- parameters[5] #rate of hospital admissions
  gammaI <- parameters[6] #recovery rate if survivor
  gammaIH <- parameters[7] #recovery rate if survivor from the hospital

  # Initial state
  params <- c(beta=beta,
              betaH=betaH,

              sigma=sigma,

              phi=phi,

              gammaH=gammaH,
              gammaI=gammaI,
              gammaIH=gammaIH)

  initial_state <- parameters[seq(8,12)]
  names(initial_state) <- c("S","E", "I", "H", "R")
  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seihr,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
