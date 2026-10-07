#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param parameters a vector of parameters. Here parameters must be input in this order beta, betaH, betaFu, sigma, phi, delta, delta2, gammaH, gammaD, gammaI, gammaIH, gammaDH,gammaFu, S0, E0, I0, H0, F0, R0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_H0,sigma_likelihood_F0,sigma_likelihood_R0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIR projection over time giving by input
#' @export

model_seihrd_tti <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaH <- parameters[2] #transmission rate from funerals

  sigma <- parameters[3] # 1/incubation time

  phi <- parameters[4] #probability of going to the hospital

  gammaH <- parameters[5] #rate of hospital admissions

  #parameters with D compartement
  delta <- parameters[6] #case fatality ratio
  deltaH <- parameters[7] #case fatality ratio in hospital


  gammaR <- parameters[8] #recovery rate if undetected

  gammaDI <- parameters[9] #death rate if detected and not in hospital
  gammaRI <- gammaDI #recovery rate if detected and not in hospital

  gammaDH <- parameters[10] #death rate if detected and in hospital
  gammaRH <- gammaDH #recovery rate if detected and in hospital


  theta = parameters[11]
  eta = parameters[12]

  # Initial state
  params <- c(beta=beta,
              betaH=betaH,

              sigma=sigma,

              phi=phi,

              gammaH=gammaH,

              delta=delta,
              deltaH=deltaH,

              gammaR=gammaR,
              gammaDI=gammaDI,
              gammaRI=gammaRI,

              gammaDH=gammaDH,
              gammaRH=gammaRH,

              theta = theta,
              eta = eta
  )

  names <- names_states("SEIHRD")$names_tti
  pars_names <- names_parameters("SEIHRD")$names_tti

  initial_state <- parameters[seq_along(names)+length(pars_names)]
  names(initial_state) <- names


  initial_state[length(names)+1] = initial_state["Id"]+initial_state["Rc"]+initial_state["Dc"]+initial_state["H"]+initial_state["RH"]+initial_state["DH"]
  names(initial_state)[length(names)+1] <- "Idtot"

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seihrd_tti,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
