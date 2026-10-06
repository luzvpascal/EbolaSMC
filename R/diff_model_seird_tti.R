#' @title Differential model for SEIRD with contact tracing, testing and isolation
#' @description
#' Differential model for SEIHR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dR, dD
#' @export
#'
diff_model_seird_tti <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- Su + Eu + Iu + Ru + Du + Id + Rd + Dd


    # dSu <- - beta_hat*c*Su*(Iu+Id)/N
    dSu <- - beta*Su*(Iu+Id)/N

    # dEu <- beta_hat*c*Su*(Iu+Id)/N - sigma*Eu
    dEu <- beta*Su*(Iu+Id)/N - sigma*Eu

    dIu <- sigma*Eu - (1-theta)*gammaR*Iu - theta*eta*Iu

    dRu <- (1-theta)*(1-delta)*gammaR*Iu

    dDu <- (1-theta)*delta*gammaR*Iu

    #detectable
    dId <- theta*eta*Iu - gammaD*Id

    dRd <- (1-delta)*gammaD * Id

    dDd <- delta*gammaD*Id

    #total detectable
    dIdtot <- theta*eta*Iu
    dDtot <- dDu+dDd

    return(list(c(dSu, dEu, dIu, dRu, dDu, dId, dRd, dDd, dIdtot, dDtot)))

    #with contact tracing
    # N <- Su + Eu + Iu + Ru + Du + Sd + Ed + Id + Rd + Dd
    #
    #
    # dSu <- - beta_hat*c*Su*Iu/N - beta_hat*(1-epsilon)*c*Su*Id/N - eta*theta*chi*NCiSu + kappa*Sd
    #
    # dEu <- beta_hat*c*Su*Iu/N + beta_hat*(1-epsilon)*c*Su*Id/N + beta_hat*(1-epsilon)*c*Sd*Iu/N + beta_hat*(1-epsilon)**2*c*Sd*Id/N - sigma*Eu - eta*theta*chi*Eu
    #
    # dIu <- sigma*Eu - gamma*Iu - theta*Iu - eta*theta*chi*Iu
    #
    # dRu <- (1-delta)*gamma*Iu
    #
    # dDu <- delta*gamma*Iu
    #
    #
    # dSd <- eta*theta*chi*NCiSu - kappa*Sd - beta_hat*(1-epsilon)*c*Sd*Iu/N - beta_hat*(1-epsilon)**2*c*Sd*Id/N
    #
    # dEd <- eta*theta*chi*Eu - sigma*Ed
    #
    # dId <- sigma*Ed + theta*Iu + eta*theta*chi*Iu - gamma*Id
    #
    # dRd <- (1-delta)*gamma * Id
    #
    # dDd <- delta*gamma*Id
    #
    # #aux var
    # dNCiSu <- (1-beta_hat)*c*Iu*Su/N + (1-beta_hat)*(1-epsilon)*c*Id*Su/N - beta_hat*c*Iu*NCiSu/N - beta_hat*(1-epsilon)*c*Id*NCiSu/N - (gamma+eta*chi*theta)*NCiSu

    # return(list(c(dSu, dEu, dIu, dRu, dDu, dSd, dEd, dId, dRd, dDd,dNCiSu)))
  })
}
