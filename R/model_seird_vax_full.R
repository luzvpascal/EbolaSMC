#' @title Differential model for SEIRD
#' @description
#' Differential model for SEIRD
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, D0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0,sigma_likelihood_D0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIRD projection over time giving by input
#' @export

model_seird_vax_full <- function(parameters, input){

  beta <- parameters[1]
  sigma <- parameters[2] #same for vaccinated and non-vaccinated
  delta <- parameters[3] #case fatality ratio non-vaccinated
  gammaR <- parameters[4] #recovery rate
  gammaD <- gammaR #death rate

  #vax parameters
  mu <- parameters[5] #vaccination coverage

  betaV <- parameters[6] #infections from vaccinated people to non-vax
  betaI <- parameters[7] #infections from non-vax people to vaccinated
  betaVV <- parameters[8] #infections from vaccinated people to vaccinated

  deltaV <- parameters[9] #case fatality ratio vaccinated

  gammaRV <- gammaR #recovery rate
  gammaDV <- gammaR #death rate

  # Initial state
  params <- c(beta=beta,
              sigma=sigma,
              delta=delta,
              gammaR=gammaR,
              gammaD=gammaD,

              mu=mu,

              betaV=betaV,
              betaI=betaI,
              betaVV=betaVV,
              deltaV=deltaV,
              gammaRV=gammaRV,
              gammaDV=gammaDV)

  all_names_states <- names_states("SEIRD")
  names <- all_names_states$names
  names_vax <- all_names_states$names_vax

  all_names_pars <- names_parameters("SEIRD")
  names_pars <- all_names_pars$names
  names_pars_vax <- all_names_pars$names_vax

  initial_state <- parameters[seq_along(names)+length(names_pars)+length(names_pars_vax)]
  names(initial_state) <- names

  initial_state[length(names)+1] = initial_state["I"]+initial_state["R"]+initial_state["D"]
  names(initial_state)[length(names)+1] <- "Itot"

  #vaccinated people at the start
  initial_state[length(initial_state)+seq_along(names_vax)] <- parameters[length(names)
                                                                      +length(names_pars)
                                                                      +length(names_pars_vax)
                                                                      +seq_along(names_vax)]
  names(initial_state)[length(names)+1+seq_along(names_vax)] <- names_vax

  #total vaccinated
  initial_state[length(initial_state)+1] = sum(parameters[length(names)
                                                      +length(names_pars)
                                                      +length(names_pars_vax)
                                                      +seq_along(names_vax)])

  names(initial_state)[length(initial_state)] <- "Vtot"

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seird_vax_full,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
