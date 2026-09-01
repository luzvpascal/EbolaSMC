#' @title Differential model for SEIRD with full vaccination chain
#' @description
#' Differential model for SEIRD with full vaccination chain
#' @param t time step
#' @param pars vector of parameters (beta, sigma, gammaR, gammaD, delta)
#' @param sigma transmission rate between E and I
#' @param state named vector of S, E, I, R, D
#' @return a list object with dS, dE, dI, dR, dD
#' @export

diff_model_seihrd_vax_full <- function(t, state, pars){

  with(as.list(c(state, pars)), {

    N <- S + E + I + H + R + RH + D + DH + V + EV + IV + HV

    ## unvaccinated susceptible
    dS <- -beta*S*I/N - betaH*S*H/N -betaV*S*IV/N -betaHV*S*HV/N-mu*S

    ## unvaccinated exposed
    dE <- beta*S*I/N + betaH*S*H/N + betaV*S*IV/N + betaHV*S*HV/N - sigma*E

    ## unvaccinated infectious
    dI <- sigma*E - (phi*gammaH + (1-phi)*delta*gammaD + (1-phi)*(1-delta)*gammaR)*I

    ## unvaccinated hospitalized
    dH <- phi*gammaH*I - (deltaH*gammaDH + (1-deltaH)*gammaRH)*H

    ####################################################
    ## vaccinated
    ####################################################
    ## vaccinated susceptible
    dV <- mu*S - betaI*I*V/N -betaVV*IV*V/N - betaIH*H*V/N-betaHVV*HV*V/N

    ## vaccinated exposed
    dEV <- betaI*I*V/N + betaVV*IV*V/N + betaIH*H*V/N + betaHVV*HV*V/N - sigma*EV

    ## vaccinated infectious
    dIV <- sigma*EV - phi*gammaH*IV - (1-phi)*deltaV*gammaDV*IV - (1-phi)*(1-deltaV)*gammaRV*IV

    ## vaccinated hospitalized
    dHV <- phi*gammaH*IV - deltaHV*gammaDH*HV - (1-deltaHV)*gammaRH*HV

    ####################################################
    ## outcomes
    ####################################################

    dR <- (1-phi)*(1-delta)*gammaR*I + (1-phi)*(1-deltaV)*gammaRV*IV


    dRH <- (1-deltaH)*gammaRH*H + (1-deltaHV)*gammaRH*HV

    dD <- (1-phi)*delta*gammaD*I + (1-phi)*deltaV*gammaDV*IV

    dDH <- deltaH*gammaDH*H + deltaHV*gammaDH*HV


    ####################################################
    ## cumulative
    ####################################################

    dItot <- sigma*E + sigma*EV

    dDtot <- dD + dDH

    return(list(c(dS, dE, dI, dH, dR, dRH, dD, dDH, dItot, dDtot, dV, dEV, dIV, dHV)))

  })

}
