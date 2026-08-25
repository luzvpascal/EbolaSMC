#' @title Differential model for SEIRD with full vaccination chain
#' @description
#' Differential model for SEIRD with full vaccination chain
#' @param t time step
#' @param pars vector of parameters (beta, sigma, gammaR, gammaD, delta)
#' @param sigma transmission rate between E and I
#' @param state named vector of S, E, I, R, D
#' @return a list object with dS, dE, dI, dR, dD
#' @export

diff_model_seird_vax_full <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + R + D + V + EV + IV


    dS <- -beta*S*I/N - betaV*IV*S/N - mu*S

    dE <- beta*S*I/N + betaV*IV*S/N - sigma*E

    dI <- sigma*E - (1-delta)*gammaR*I - delta*gammaD*I


    dV <- mu*S - betaI*I*V/N  - betaVV*IV*V/N

    dEV <- betaI*I*V/N  + betaVV*IV*V/N - sigma*EV

    dIV <- sigma*EV - (1-deltaV)*gammaRV*IV - deltaV*gammaDV*IV



    dR <- (1-delta)*gammaR*I + (1-deltaV)*gammaRV*IV

    dD <- delta*gammaD*I + deltaV*gammaDV*IV

    dItot <- sigma*E + sigma*EV


    return(list(c(dS, dE, dI, dR, dD, dItot, dV, dEV, dIV)))
  })
}
