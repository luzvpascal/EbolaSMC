#' @title Differential model for SEIHR
#' @description
#' Differential model for SEIHR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dH, dR
#' @export

diff_model_seihfrd_vax_full <- function(t, state, pars){

  with(as.list(c(state, pars)), {

    N <- S + E + I + H + Fu + V + EV + IV + HV + FuV + R + RH + D

    ################################################
    # UNVACCINATED
    ################################################

    dS <- -beta*S*I/N - betaV*IV*S/N - betaH*S*H/N - betaHV*S*HV/N - betaFu*S*Fu/N - betaFuV*S*FuV/N - mu*S


    dE <- beta*S*I/N + betaV*IV*S/N + betaFu*S*Fu/N + betaFuV*S*FuV/N + betaH*S*H/N + betaHV*S*HV/N - sigma*E


    dI <- sigma*E - phi*gammaH*I - (1-phi)*delta*gammaD*I - (1-phi)*(1-delta)*gammaR*I

    ## unvaccinated hospitalized
    dH <- phi*gammaH*I - deltaH*gammaDH*H - (1-deltaH)*gammaRH*H

    #funeral unvaccinated
    dFu <- deltaH*gammaDH*H + (1-phi)*delta*gammaD*I - gammaFu*Fu

    ################################################
    # VACCINATED
    ################################################

    dV <- mu*S - betaI*I*V/N - betaVV*IV*V/N - betaIH*H*V/N-betaHVV*HV*V/N - betaVFu*V*Fu/N  - betaFuVV*V*Fu/N


    dEV <- betaI*I*V/N + betaVV*IV*V/N + betaIH*H*V/N-betaHVV*HV*V/N + betaVFu*V*Fu/N + betaFuVV*V*Fu/N - sigma*EV

    dIV <- sigma*EV - phi*gammaH*IV - (1-phi)*deltaV*gammaDV*IV - (1-phi)*(1-deltaV)*gammaRV*IV

    ## vaccinated hospitalized
    dHV <- phi*gammaH*IV - deltaHV*gammaDH*HV - (1-deltaHV)*gammaRH*HV

    dFuV <- (1-phi)*deltaV*gammaDV*IV + deltaHV*gammaDH*HV - gammaFu*FuV


    ################################################
    # OUTCOMES
    ################################################

    dR <- (1-phi)*(1-delta)*gammaR*I + (1-phi)*(1-deltaV)*gammaRV*IV


    dRH <- (1-deltaH)*gammaRH*H +(1-deltaHV)*gammaRH*HV

    dD <- gammaFu*Fu + gammaFu*FuV


    ################################################
    # CUMULATIVE OUTPUTS
    ################################################

    dItot <- sigma*E + sigma*EV

    dDinst <- deltaH*gammaDH*H + (1-phi)*delta*gammaD*I + (1-phi)*deltaV*gammaDV*IV + deltaHV*gammaDH*HV

    dVtot <- mu*S

    return(list(c(dS,dE,dI,dH,dFu,dR,dRH,dD,dItot,dDinst,dV,dEV,dIV,dHV,dFuV,dVtot)))

  })
}
