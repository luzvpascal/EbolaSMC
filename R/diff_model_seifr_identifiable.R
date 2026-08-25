#' @title Differential model for SEIFR
#' @description
#' Differential model for SEIFR
#' @param t time of prediction
#' @param state named vector of S, E, I, Fu, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dFu, dR
#' @export

diff_model_seifr_identifiable <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + Fu + R

    dS <- -beta*S*I/N - betaFu*S*Fu/N

    dE <- beta*S*I/N + betaFu*S*Fu/N - sigma*E

    dI <- sigma*E - (delta_gammaD + one_minus_delta_gammaI)*I

    dFu <- delta_gammaD*I - gammaFu*Fu

    dR <- one_minus_delta_gammaI*I + gammaFu*Fu

    return(list(c(dS, dE, dI, dFu, dR)))
  })
}
