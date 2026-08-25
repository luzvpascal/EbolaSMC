#' @title Differential model for SEIRD
#' @description
#' Differential model for SEIRD
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, D0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0,sigma_likelihood_D0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIRD projection over time giving by input
#' @export

model_seifrd_vax_full <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaFu <- parameters[2] #transmission rate from funerals
  sigma <- parameters[3] # 1/incubation time
  delta <- parameters[4] #case fatality ratio
  gammaR <- parameters[5] #death rate if infected
  gammaD <- gammaR  #recovery rate if survivor
  gammaFu <- parameters[6] #1/funeral time

  #vax parameters
  mu <- parameters[7] #vaccination rate

  betaV <- parameters[8] #infections from vaccinated people to non-vax
  betaI <- parameters[9] #infections from non-vax people to vaccinated
  betaVV <- parameters[10] #infections from vaccinated people to vaccinated

  betaFuV <- parameters[11] # infections from vaccindated funeral to non-vax person
  betaVFu <- parameters[12] # infections from unvaccindated funeral to vax person
  betaVFuVFu <- parameters[13] # infections from vaccindated funeral to vax person

  deltaV <- parameters[14] #case fatality ratio vaccinated

  gammaRV <- gammaD #recovery rate
  gammaDV <- gammaD #death rate


  # Initial state
  params <- c(beta=beta,
              betaFu=betaFu,
              sigma=sigma,
              delta=delta,

              gammaR=gammaR,
              gammaD=gammaD,
              gammaFu=gammaFu,

              mu=mu,

              betaV=betaV,
              betaI=betaI,
              betaVV=betaVV,

              betaFuV=betaFuV,# infections from vaccindated funeral to non-vax person
              betaVFu=betaVFu, # infections from unvaccindated funeral to vax person
              betaVFuVFu=betaVFuVFu, # infections from vaccindated funeral to vax person

              deltaV=deltaV,
              gammaRV=gammaRV,
              gammaDV=gammaDV
              )

  names <- c("S","E", "I", "Fu", "R", "D")
  initial_state <- parameters[seq_along(names)+14]
  names(initial_state) <- names

  initial_state[length(names)+1] = initial_state["I"]+initial_state["R"]+initial_state["D"]
  names(initial_state)[length(names)+1] <- "Itot"

  #add 0 vaccinated people at the start
  initial_state[length(initial_state)+seq(4)] <- 0
  names(initial_state)[length(initial_state)-seq(3,0)] <- c("V", "EV", "IV", "FuV")

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seifrd_vax_full,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
