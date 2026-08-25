#' @title Differential model for SEIHFR
#' @description
#' Differential model for SEIHFR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta, betaH, betaFu, sigma, phi, delta, delta2, gammaH, gammaD, gammaI, gammaIH, gammaDH,gammaFu, S0, E0, I0, H0, F0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_H0,sigma_likelihood_F0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihfr <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaH <- parameters[2] #transmission rate from funerals
  betaFu <- parameters[3] #transmission rate from funerals

  sigma <- parameters[4] # 1/incubation time

  phi <- parameters[5] #probability of going to the hospital
  delta <- parameters[6] #case fatality ratio
  delta2 <- parameters[7] #case fatality ratio in hospital

  gammaH <- parameters[8] #rate of hospital admissions
  gammaD <- parameters[9] #fatality rate if infected not going to the hospital
  gammaI <- parameters[10] #recovery rate if survivor
  gammaIH <- parameters[11] #recovery rate if survivor from the hospital
  gammaDH <- parameters[12] #fatality rate in the hospital
  gammaFu <- parameters[13] #1/funeral time

  # Initial state
  params <- c(beta=beta,
              betaH=betaH,
              betaFu=betaFu,

              sigma=sigma,

              phi=phi,
              delta=delta,
              delta2=delta2,

              gammaH=gammaH,
              gammaD=gammaD,
              gammaI=gammaI,
              gammaIH=gammaIH,
              gammaDH=gammaDH,
              gammaFu=gammaFu)

  initial_state <- parameters[seq(14,19)]
  names(initial_state) <- c("S","E", "I", "H", "Fu", "R")
  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seihfr,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
