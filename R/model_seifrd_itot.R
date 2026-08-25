#' @title Differential model for SEIRD
#' @description
#' Differential model for SEIRD
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, D0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0,sigma_likelihood_D0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIRD projection over time giving by input
#' @export

model_seifrd_itot <- function(parameters, input){

  beta <- parameters[1] #transmission rate from infected
  betaFu <- parameters[2] #transmission rate from funerals
  sigma <- parameters[3] # 1/incubation time
  delta <- parameters[4] #case fatality ratio
  gammaD <- parameters[5] #death rate if infected
  gammaI <- gammaD#parameters[6] #recovery rate if survivor
  gammaFu <- parameters[6] #1/funeral time

  # Initial state
  params <- c(beta=beta,
              betaFu=betaFu,
              sigma=sigma,
              delta=delta,
              gammaD=gammaD,
              gammaI=gammaI,
              gammaFu=gammaFu)

  initial_state <- parameters[seq(7,12)]
  names <- c("S","E", "I", "Fu", "R", "D")
  names(initial_state) <- names

  initial_state[length(names)+1] = initial_state["I"]+initial_state["Fu"]+initial_state["R"]+initial_state["D"]
  names(initial_state)[length(names)+1] <- "Itot"

  initial_state[length(initial_state)+1] = initial_state["Fu"]+initial_state["D"]
  names(initial_state)[length(initial_state)] <- "Dinst"

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seifrd_itot,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
