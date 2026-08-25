#' @title Differential model for SEIRV
#' @description
#' Differential model for SEIRV
#' @param t time step
#' @param pars vector of parameters (beta, sigma, gamma, mu)
#' @param sigma transmission rate between E and I
#' @param state named vector of S, E, I, R, V
#' @return a list object with dS, dE, dI, dR, dV
#' @export

diff_model_seir_vax <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + R + V

    dS <- -beta * S * I / N - mu*S
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I
    dV <- mu*S

    return(list(c(dS, dE, dI, dR, dV)))
  })
}
