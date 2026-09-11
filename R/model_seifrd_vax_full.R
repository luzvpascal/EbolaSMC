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
  betaFuVV <- parameters[13] # infections from vaccindated funeral to vax person

  deltaV <- parameters[14] #case fatality ratio vaccinated

  gammaRV <- gammaR #recovery rate vaccinated
  gammaDV <- gammaR #death rate vaccinated


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
              betaFuVV=betaFuVV, # infections from vaccindated funeral to vax person

              deltaV=deltaV,
              gammaRV=gammaRV,
              gammaDV=gammaDV
              )

  all_names_states <- names_states("SEIFRD")
  names <- all_names_states$names
  names_vax <- all_names_states$names_vax

  all_names_pars <- names_parameters("SEIFRD")
  names_pars <- all_names_pars$names
  names_pars_vax <- all_names_pars$names_vax

  initial_state <- parameters[seq_along(names)+length(names_pars)+length(names_pars_vax)]
  names(initial_state) <- names

  initial_state[length(names)+1] = initial_state["I"]+initial_state["R"]+initial_state["D"]
  names(initial_state)[length(names)+1] <- "Itot"

  initial_state[length(initial_state)+1] <- initial_state["D"] + initial_state["Fu"]
  names(initial_state)[length(initial_state)] <- "Dinst"

  #vaccinated people at the start
  initial_state[length(initial_state)+seq_along(names_vax)] <- parameters[length(names)
                                                                      +length(names_pars)
                                                                      +length(names_pars_vax)
                                                                      +seq_along(names_vax)]
  names(initial_state)[length(names)+2+seq_along(names_vax)] <- names_vax

  #total vaccinated
  initial_state[length(initial_state)+1] = sum(parameters[length(names)
                                                          +length(names_pars)
                                                          +length(names_pars_vax)
                                                          +seq_along(names_vax)])

  names(initial_state)[length(initial_state)] <- "Vtot"

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seifrd_vax_full,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
