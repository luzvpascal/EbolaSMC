#' @title Differential model for SEIR
#' @description
#' Differential model for SEIR
#' @param t time step
#' @param pars vector of parameters (beta, sigma, gamma)
#' @param sigma transmission rate between E and I
#' @param state named vector of S, E, I, R
#' @return a list object with dS, dE, dI, dR
#' @export

diff_model_seir <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + R

    dS <- -beta * S * I / N
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I

    return(list(c(dS, dE, dI, dR)))
  })
}
