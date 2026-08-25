#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dH, dR
#' @export

diff_model_seifrd_itot <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + Fu + R + D

    dS <- -beta*S*I/N - betaFu*S*Fu/N

    dE <- beta*S*I/N + betaFu*S*Fu/N - sigma*E

    dI <- sigma*E - (delta*gammaD + (1-delta)*gammaI)*I

    dFu <- delta*gammaD*I - gammaFu*Fu

    dR <- (1-delta)*gammaI*I

    dD <- gammaFu*Fu

    dItot <- sigma*E

    dDinst <- delta*gammaD*I

    return(list(c(dS, dE, dI, dFu, dR, dD, dItot, dDinst)))
  })
}
