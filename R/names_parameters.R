#' @title Get names of parameters in model
#' @description
#' Get names of parameters in model
#' @param model_type string describing model type
#' @return Names of parameters
#' @export
#'
names_parameters <- function(model_type){
  if (model_type=="SEIRD"){
    names <- c("beta",
               "sigma",
               "delta",
               "gammaR")
    names_vax <- c("mu",
                   "betaV",
                   "betaI",
                   "betaVV",
                   "deltaV")

    names_tti <- c("beta",
                   # "c",
                   "sigma",
                   "delta",
                   "gammaR",
                   "gammaD",
                   "theta",
                   "eta"
                   )
    # names_tti <- c("beta_hat",
    #                "c",
    #                "sigma",
    #                "delta",
    #                "gammaR",
    #                "gammaD",
    #                "theta",
    #                "eta"
    #                )

  return(list(names=names,
              names_vax=names_vax,
              names_tti=names_tti))
  } else if (model_type=="SEIFRD"){
    names <- c("beta",
               "betaFu",
               "sigma",
               "delta",
               "gammaR",
               "gammaFu")
    names_vax <- c("mu",
                   "betaV",
                   "betaI",
                   "betaVV",
                   "betaFuV",
                   "betaVFu",
                   "betaFuVV",
                   "deltaV")

  return(list(names=names,
              names_vax=names_vax))
  } else if (model_type=="SEIHRD"){
    names <- c("beta",
               "betaH",
               "sigma",
               "phi",
               "gammaH",
               "delta",
               "deltaH",
               "gammaR",
               "gammaRH")
    names_vax <- c("mu",
                   "betaV",
                   "betaI",
                   "betaVV",
                   "betaHV",
                   "betaIH",
                   "betaHVV",
                   "deltaV",
                   "deltaHV")

  return(list(names=names,
              names_vax=names_vax))
  } else if (model_type=="SEIHFRD"){
    names <- c("beta",
               "betaH",
               "betaFu",
               "sigma",
               "phi",
               "gammaH",
               "delta",
               "deltaH",
               "gammaR",
               "gammaRH",
               "gammaFu"
               )
    names_vax <- c("mu",
                   "betaV",
                   "betaI",
                   "betaVV",
                   "betaHV",
                   "betaIH",
                   "betaHVV",
                   "betaFuV",
                   "betaVFu",
                   "betaFuVV",
                   "deltaV",
                   "deltaHV")

  return(list(names=names,
              names_vax=names_vax))
  }
}
