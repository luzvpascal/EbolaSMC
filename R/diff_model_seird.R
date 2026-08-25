#' @title Differential model for SEIR
#' @description
#' Differential model for SEIRD
#' @param t time step
#' @param pars vector of parameters (beta, sigma, gammaR, gammaD, delta)
#' @param sigma transmission rate between E and I
#' @param state named vector of S, E, I, R, D
#' @return a list object with dS, dE, dI, dR, dD
#' @export

diff_model_seird <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + R + D

    dS <- -beta * S * I / N
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - (1-delta)*gammaR * I - delta*gammaD*I
    dR <- (1-delta)*gammaR * I
    dD <- delta*gammaD*I

    return(list(c(dS, dE, dI, dR, dD)))
  })
}
