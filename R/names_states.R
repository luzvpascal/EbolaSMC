#' @title Get names of states in model
#' @description
#' Get names of states in model
#' @param model_type string describing model type
#' @return list like object of $names (names of states to call model functions)  and $all_names: all the names of states returned by model function
#' @export
#'
names_states <- function(model_type){
  if (model_type=="SEIRD"){
    names <- c("S","E","I","R","D")
    all_names <- c("S","E","I","R","D","Itot")
    names_vax <- c("V","EV","IV")
  } else if (model_type=="SEIFRD"){
    names <- c("S","E","I","Fu","R","D")
    all_names <- c("S","E","I","Fu","R","D","Itot","Dinst")
    names_vax <- c("V", "EV", "IV", "FuV")
  } else if (model_type=="SEIHRD"){
    names <- c("S", "E", "I", "H", "R", "RH", "D", "DH")
    all_names <- c("S", "E", "I", "H", "R", "RH", "D", "DH", "Itot", "Dtot")
    names_vax <- c("V", "EV","IV", "HV")
  } else if (model_type=="SEIHFRD"){
    names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D")
    all_names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D", "Itot", "Dinst", "DinstH")
    names_vax <- c("V","EV","IV", "HV", "FuV" )
  }

  return(list(names=names,
              all_names=all_names,
              names_vax=names_vax))
}
