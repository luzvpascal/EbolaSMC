#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dH, dR
#' @export

diff_model_seihr <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + H + R

    dS <- -beta*S*I/N - betaH*S*H/N

    dE <- beta*S*I/N + betaH*S*H/N - sigma*E

    dI <- sigma*E - (phi*gammaH + (1-phi)*gammaI)*I

    dH <- phi*gammaH*I - gammaIH*H

    dR <- (1-phi)*gammaI*I + gammaIH*H

    return(list(c(dS, dE, dI, dH, dR)))
  })
}
