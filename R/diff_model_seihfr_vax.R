#' @title Differential model for SEIHFR
#' @description
#' Differential model for SEIR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, Fu, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dFu, dR
#' @export

diff_model_seihfr_vax <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- S + E + I + H + Fu + R + V

    dS <- -beta*S*I/N - betaH*S*H/N - betaFu*S*Fu/N -mu*S
    #exposed from I,   from H,      from Fu

    dE <- beta*S*I/N + betaH*S*H/N + betaFu*S*Fu/N - sigma*E
    #exposed from I,   from H,      from Fu         #leaves E to I

    dI <- sigma*E - #from E
      (phi*gammaH + (1-phi)*delta*gammaD +  (1-phi)*(1-delta)*gammaI)*I
    #leave to H       leaves to F             leaves to R

    dH <- phi*gammaH*I - ((1-delta2)*gammaIH + delta2*gammaDH)*H
    #from I         leaves to R             leaves to F

    dFu <- (1-phi)*delta*gammaD*I + delta2*gammaDH*H - gammaFu*Fu
    #from I                 from H              leaves to R

    dR <- (1-phi)*(1-delta)*gammaI*I + (1-delta2)*gammaIH*H + gammaFu*Fu
    #from I                     from H              from F

    dV <- mu*S

    return(list(c(dS, dE, dI, dH, dFu, dR, dV)))
  })
}
