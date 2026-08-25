#' @title Differential model for SEIRD
#' @description
#' Differential model for SEIRD
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, D0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0,sigma_likelihood_D0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIRD projection over time giving by input
#' @export

model_seird_itot <- function(parameters, input){

  beta <- parameters[1]
  sigma <- parameters[2]
  delta <- parameters[3]
  gammaR <- parameters[4]
  gammaD <- gammaR #same as gamma R

  # Initial state
  params <- c(beta=beta,
              sigma=sigma,
              delta=delta,
              gammaR=gammaR,
              gammaD=gammaD)

  names <- c("S","E", "I", "R","D")
  initial_state <- parameters[seq_along(names)+4]

  names(initial_state) <- names
  initial_state[length(names)+1] = initial_state["I"]+initial_state["R"]+initial_state["D"]
  names(initial_state)[length(names)+1] <- "Itot"

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seird_itot,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
