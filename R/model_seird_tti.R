#' @title Differential model for SEIRD with testing
#' @description
#' Differential model for SEIRD
#' @param parameters a vector of parameters. Here parameters must be input in this order beta,sigma,gamma,S0, E0, I0, R0, D0, sigma_likelihood_S0, sigma_likelihood_E0,sigma_likelihood_I0,sigma_likelihood_R0,sigma_likelihood_D0
#' @param input a vector of inputs (time).
#' @return a data.frame with SEIRD projection over time giving by input
#' @export

model_seird_tti <- function(parameters, input){

  # Parameters
  # params <- c(beta_hat = parameters[1],
  #             c = parameters[2],
  #             sigma = parameters[3],
  #             delta = parameters[4],
  #             gammaR = parameters[5],
  #             gammaD = parameters[6],
  #             theta = parameters[7],
  #             eta = parameters[8])
  params <- c(beta = parameters[1],
              sigma = parameters[2],
              delta = parameters[3],
              gammaR = parameters[4],
              gammaD = parameters[5],
              theta = parameters[6],
              eta = parameters[7])

  names <- names_states("SEIRD")$names_tti
  initial_state <- parameters[seq_along(names)+length(params)]

  #detectable infectious
  names(initial_state) <- names
  initial_state[length(names)+1] = initial_state["Iu"]+initial_state["Ru"]+initial_state["Du"]+
                                    initial_state["Id"]+initial_state["Rd"]+initial_state["Dd"]

  names(initial_state)[length(names)+1] <- "Idtot"

  #total deaths
  initial_state[length(initial_state)+1] = initial_state["Du"]+initial_state["Dd"]
  names(initial_state)[length(initial_state)] <- "Dtot"

  # Solve
  output <- deSolve::ode(y = initial_state, times = input,
                         func = EbolaSMC::diff_model_seird_tti,
                         parms = params)

  output <- as.data.frame(output)
  output <- output[,-1]
  return(output)
}
