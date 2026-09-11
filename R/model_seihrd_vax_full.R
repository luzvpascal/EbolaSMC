#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta, betaH, betaFu, sigma, phi, delta, delta2, gammaH, gammaD, gammaI, gammaIH, gammaDH,gammaFu, S0, E0, I0, H0, F0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_H0,sigma_likelihood_F0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihrd_vax_full <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaH <- parameters[2] #transmission rate from hospital

  sigma <- parameters[3]# 1/incubation time

  phi <- parameters[4] #probability of going to the hospital

  gammaH <- parameters[5] #rate of hospital admissions

  #parameters with D compartement
  delta <- parameters[6] #case fatality ratio
  deltaH <- parameters[7] #case fatality ratio in hospital

  gammaR <- parameters[8] #death rate if not in hospital
  gammaD <- gammaR #recovery rate if survivor not in hospital

  gammaDH <- parameters[9]#death rate if in hospital
  gammaRH <- gammaDH #recovery rate if in hospital

  ##################################################
  # vaccination parameters
  ##################################################

  mu <- parameters[10] #vaccination rate

  betaV  <- parameters[11]#infections from vaccinated people to non-vax
  betaI  <- parameters[12]#infections from non-vax people to vaccinated
  betaVV <- parameters[13]#infections from vaccinated people to vaccinated

  betaHV  <- parameters[14] # infections from vaccindated hospitalized to non-vax person
  betaIH  <- parameters[15] # infections from non-vax hospitalized to vax person
  betaHVV <- parameters[16] # infections from vax hospitalized to vax person

  deltaV  <- parameters[17]#case fatality ratio vaccinated
  deltaHV <- parameters[18]#case fatality ratio vaccinated in hospital

  ##################################################
  # vaccinated progression rates
  ##################################################

  gammaRV <- gammaR #recovery rate vaccinated
  gammaDV <- gammaR #death rate vaccinated

  gammaRHV <- gammaRH #recovery rate vaccinated in hospital
  gammaDHV <- gammaRH #death rate vaccinated in hospital

  ##################################################
  # parameter vector
  ##################################################

  params <- c(
    beta=beta,
    betaH=betaH,

    sigma=sigma,

    phi=phi,

    gammaH=gammaH,

    delta=delta,
    deltaH=deltaH,

    gammaR=gammaR,
    gammaD=gammaD,

    gammaDH=gammaDH,
    gammaRH=gammaRH,

    mu=mu,

    betaV=betaV,
    betaI=betaI,
    betaVV=betaVV,

    betaHV=betaHV,
    betaIH=betaIH,
    betaHVV=betaHVV,

    deltaV=deltaV,
    deltaHV=deltaHV,

    gammaDV=gammaDV,
    gammaRV=gammaRV,

    gammaDHV=gammaDHV,
    gammaRHV=gammaRHV
  )

  ##################################################
  # initial conditions
  ##################################################

  all_names_states <- names_states("SEIHRD")
  names <- all_names_states$names
  names_vax <- all_names_states$names_vax

  all_names_pars <- names_parameters("SEIHRD")
  names_pars <- all_names_pars$names
  names_pars_vax <- all_names_pars$names_vax

  initial_state <- parameters[seq_along(names)+length(names_pars)+length(names_pars_vax)]
  names(initial_state) <- names

  ##################################################
  # cumulative incidence
  ##################################################

  initial_state[length(initial_state)+1] <-
    initial_state["I"] +
    initial_state["H"] +
    initial_state["R"] +
    initial_state["RH"] +
    initial_state["D"] +
    initial_state["DH"]

  names(initial_state)[length(initial_state)] <- "Itot"

  ##################################################
  # cumulative deaths
  ##################################################

  initial_state[length(initial_state)+1] <-
    initial_state["D"] +
    initial_state["DH"]

  names(initial_state)[length(initial_state)] <- "Dtot"

  ##################################################
  # vaccinated compartments start at 0
  ##################################################
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

  ##################################################
  # solve ODE
  ##################################################
  output <- deSolve::ode(
    y = initial_state,
    times = input,
    func = EbolaSMC::diff_model_seihrd_vax_full,
    parms = params
  )

  output <- as.data.frame(output)
  output <- output[, -1]

  return(output)
}
