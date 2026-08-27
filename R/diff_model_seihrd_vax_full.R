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

    N <- S + E + I + H +
      V + EV + IV + HV +
      R + RH + D + DH

    ## unvaccinated susceptible

    dS <- -beta*S*I/N
    -betaV*S*IV/N
    -betaH*S*H/N
    -betaHV*S*HV/N
    -mu*S

    ## unvaccinated exposed

    dE <- beta*S*I/N
    +betaV*S*IV/N
    +betaH*S*H/N
    +betaHV*S*HV/N
    -sigma*E

    ## unvaccinated infectious

    dI <- sigma*E -
      (phi*gammaH +
         (1-phi)*delta*gammaDI +
         (1-phi)*(1-delta)*gammaRI)*I

    ## unvaccinated hospitalized

    dH <- phi*gammaH*I -
      (deltaH*gammaDH +
         (1-deltaH)*gammaRH)*H


    ####################################################
    ## vaccinated
    ####################################################

    dV <- mu*S
    -betaI*I*V/N
    -betaVV*IV*V/N
    -betaIH*H*V/N
    -betaHVV*HV*V/N

    dEV <- betaI*I*V/N
    +betaVV*IV*V/N
    +betaIH*H*V/N
    +betaHVV*HV*V/N
    -sigma*EV

    dIV <- sigma*EV -
      (phi*gammaH +
         (1-phi)*deltaV*gammaDI +
         (1-phi)*(1-deltaV)*gammaRI)*IV

    dHV <- phi*gammaH*IV -
      (deltaHV*gammaDH +
         (1-deltaHV)*gammaRH)*HV


    ####################################################
    ## outcomes
    ####################################################

    dR <- (1-phi)*(1-delta)*gammaRI*I +
      (1-phi)*(1-deltaV)*gammaRI*IV

    dRH <- (1-deltaH)*gammaRH*H +
      (1-deltaHV)*gammaRH*HV

    dD <- (1-phi)*delta*gammaDI*I +
      (1-phi)*deltaV*gammaDI*IV

    dDH <- deltaH*gammaDH*H +
      deltaHV*gammaDH*HV


    ####################################################
    ## cumulative
    ####################################################

    dItot <- sigma*E + sigma*EV

    dDtot <- dD + dDH

    return(
      list(
        c(
          dS,dE,dI,dH,
          dR,dRH,dD,dDH,
          dItot,dDtot,
          dV,dEV,dIV,dHV
        )
      )
    )

  })

}
