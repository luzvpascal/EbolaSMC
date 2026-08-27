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

    N <- S + E + I + H + Fu +
      V + EV + IV + HV + FuV +
      R + RH + D


    ################################################
    # UNVACCINATED
    ################################################

    dS <- -beta*S*I/N
    -betaV*S*IV/N
    -betaH*S*H/N
    -betaHV*S*HV/N
    -betaFu*S*Fu/N
    -betaFuV*S*FuV/N
    -mu*S


    dE <- beta*S*I/N
    +betaV*S*IV/N
    +betaH*S*H/N
    +betaHV*S*HV/N
    +betaFu*S*Fu/N
    +betaFuV*S*FuV/N
    -sigma*E


    dI <- sigma*E -
      (phi*gammaH +
         (1-phi)*delta*gammaDI +
         (1-phi)*(1-delta)*gammaRI)*I


    dH <- phi*gammaH*I -
      (deltaH*gammaDH +
         (1-deltaH)*gammaRH)*H


    dFu <- deltaH*gammaDH*H +
      (1-phi)*delta*gammaDI*I -
      gammaFu*Fu


    ################################################
    # VACCINATED
    ################################################

    dV <- mu*S
    -betaI*I*V/N
    -betaVV*IV*V/N
    -betaIH*H*V/N
    -betaHVV*HV*V/N
    -betaVFu*Fu*V/N
    -betaFuVV*FuV*V/N


    dEV <- betaI*I*V/N
    +betaVV*IV*V/N
    +betaIH*H*V/N
    +betaHVV*HV*V/N
    +betaVFu*Fu*V/N
    +betaFuVV*FuV*V/N
    -sigma*EV


    dIV <- sigma*EV -
      (phi*gammaH +
         (1-phi)*deltaV*gammaDI +
         (1-phi)*(1-deltaV)*gammaRI)*IV


    dHV <- phi*gammaH*IV -
      (deltaHV*gammaDH +
         (1-deltaHV)*gammaRH)*HV


    dFuV <- deltaHV*gammaDH*HV +
      (1-phi)*deltaV*gammaDI*IV -
      gammaFu*FuV


    ################################################
    # OUTCOMES
    ################################################

    dR <- (1-phi)*(1-delta)*gammaRI*I +
      (1-phi)*(1-deltaV)*gammaRI*IV


    dRH <- (1-deltaH)*gammaRH*H +
      (1-deltaHV)*gammaRH*HV


    dD <- gammaFu*Fu +
      gammaFu*FuV


    ################################################
    # CUMULATIVE OUTPUTS
    ################################################

    dItot <- sigma*E + sigma*EV


    dDinst <- deltaH*gammaDH*H +
      (1-phi)*delta*gammaDI*I +
      deltaHV*gammaDH*HV +
      (1-phi)*deltaV*gammaDI*IV


    return(
      list(
        c(
          dS,dE,dI,dH,dFu,
          dR,dRH,dD,
          dItot,dDinst,
          dV,dEV,dIV,dHV,dFuV
        )
      )
    )

  })

}
