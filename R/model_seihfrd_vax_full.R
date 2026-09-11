#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta, betaH, betaFu, sigma, phi, delta, delta2, gammaH, gammaD, gammaI, gammaIH, gammaDH,gammaFu, S0, E0, I0, H0, F0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_H0,sigma_likelihood_F0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihfrd_vax_full <- function(parameters, input){

  beta <- parameters[1]
  betaH <- parameters[2]
  betaFu <- parameters[3]

  sigma <- parameters[4]

  # hospital admissions
  phi <- parameters[5]
  gammaH <- parameters[6]

  # mortality
  delta <- parameters[7]
  deltaH <- parameters[8]

  gammaD <- parameters[9]
  gammaR <- gammaD

  gammaDH <- parameters[10]
  gammaRH <- gammaDH

  # funeral
  gammaFu <- parameters[11]

  ##################################################
  # vaccination parameters
  ##################################################

  mu <- parameters[12]

  # community transmission
  betaV <- parameters[13]
  betaI <- parameters[14]
  betaVV <- parameters[15]

  # hospital transmission
  betaHV <- parameters[16]
  betaIH <- parameters[17]
  betaHVV <- parameters[18]

  # funeral transmission
  betaFuV <- parameters[19]
  betaVFu <- parameters[20]
  betaFuVV <- parameters[21]

  # mortality
  deltaV <- parameters[22]
  deltaHV <- parameters[23]

  ##################################################
  # vaccinated rates
  ##################################################

  gammaDV <- gammaD
  gammaRV <- gammaR

  gammaDHV <- gammaDH
  gammaRHV <- gammaRH

  ##################################################
  # parameter vector
  ##################################################

  params <- c(
    beta = beta,
    betaH = betaH,
    betaFu = betaFu,

    sigma = sigma,

    phi = phi,
    gammaH = gammaH,

    delta = delta,
    deltaH = deltaH,

    gammaD = gammaD,
    gammaR = gammaR,

    gammaDH = gammaDH,
    gammaRH = gammaRH,

    gammaFu = gammaFu,

    mu = mu,

    betaV = betaV,
    betaI = betaI,
    betaVV = betaVV,

    betaHV = betaHV,
    betaIH = betaIH,
    betaHVV = betaHVV,

    betaFuV = betaFuV,
    betaVFu = betaVFu,
    betaFuVV = betaFuVV,

    deltaV = deltaV,
    deltaHV = deltaHV,

    gammaDV = gammaDV,
    gammaRV = gammaRV,

    gammaDHV = gammaDHV,
    gammaRHV = gammaRHV
  )

  ##################################################
  # initial conditions
  ##################################################

  all_names_states <- names_states("SEIHFRD")
  names <- all_names_states$names
  names_vax <- all_names_states$names_vax

  all_names_pars <- names_parameters("SEIHFRD")
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
    initial_state["Fu"] +
    initial_state["R"] +
    initial_state["RH"] +
    initial_state["D"]

  names(initial_state)[length(initial_state)] <- "Itot"

  ##################################################
  # instantaneous deaths
  ##################################################

  initial_state[length(initial_state)+1] <-
    initial_state["D"] +
    initial_state["Fu"]

  names(initial_state)[length(initial_state)] <- "Dinst"

  ##################################################
  # vaccinated compartments start at zero
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
    func = EbolaSMC::diff_model_seihfrd_vax_full,
    parms = params
  )

  output <- as.data.frame(output)
  output <- output[, -1]

  return(output)

}
