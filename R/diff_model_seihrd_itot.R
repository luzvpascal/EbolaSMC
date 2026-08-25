#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dH, dR
#' @export

diff_model_seihrd_itot <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + H + R + RH + D + DH

    dS <- -beta*S*I/N - betaH*S*H/N

    dE <- beta*S*I/N + betaH*S*H/N - sigma*E

    dI <- sigma*E - (phi*gammaH + (1-phi)*delta*gammaDI + (1-phi)*(1-delta)*gammaRI)*I

    dH <- phi*gammaH*I - (deltaH*gammaDH + (1-deltaH)*gammaRH)*H

    dR <- (1-phi)*(1-delta)*gammaRI*I

    dRH <- (1-deltaH)*gammaRH*H

    dD <- (1-phi)*delta*gammaDI*I

    dDH <- deltaH*gammaDH*H

    dItot <- sigma*E

    dDtot <- dD + dDH

    return(list(c(dS, dE, dI, dH, dR, dRH, dD, dDH, dItot, dDtot)))
  })
}
