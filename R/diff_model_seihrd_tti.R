#' @title Differential model for SEHIRD with contact tracing, testing and isolation
#' @description
#' Differential model for SEIHR
#' @param t time of prediction
#' @param state named vector of S, E, I, H, R
#' @param pars vector of parameters
#' @return a list object with dS, dE, dI, dR, dD
#' @export
#'
diff_model_seihrd_tti <- function(t, state, pars){
  with(as.list(c(state, pars)), {

    N <- Su + Eu + Iu + Ru + Du + Id + Rd + Dd


    # dSu <- - beta_hat*c*Su*(Iu+Id)/N
    dSu <- - beta*Su*(Iu+Id+Ic)/N - betaH*H*Su/N

    # dEu <- beta_hat*c*Su*(Iu+Id)/N - sigma*Eu
    dEu <- beta*Su*(Iu+Id)/N + betaH*H*Su/N - sigma*Eu

    dIu <- sigma*Eu - (1-theta)*gammaR*Iu - theta*eta*Iu

    dRu <- (1-theta)*(1-delta)*gammaR*Iu

    dDu <- (1-theta)*delta*gammaR*Iu #we have data about


    #detected
    dId <- theta*eta*Iu - phi*gammaH*Id - (1-phi)*delta*gammaDI*Id - (1-phi)*(1-delta)*gammaRI*Id

    #back in the community
    dRc <- (1-phi)*(1-delta)*gammaRI*Id #no data

    dDc <- (1-phi)*delta*gammaDI*Id #no data

    #hospital
    dH <- phi*gammaH*Id - (deltaH*gammaDH + (1-deltaH)*gammaRH)*H

    dRH <- (1-deltaH)*gammaRH*H #we have data

    dDH <- deltaH*gammaDH*H #we have data


    #total detected alive
    dIdtot <- theta*eta*Iu

    return(list(c(dSu,
                  dEu,
                  dIu,
                  dRu,
                  dDu,

                  dId,

                  dRc,
                  dDc,

                  dH,
                  dRH,
                  dDH,

                  dIdtot
                  )))

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
