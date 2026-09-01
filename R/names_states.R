#' @title Get names of states in model
#' @description
#' Get names of states in model
#' @param model_type string describing model type
#' @return list like object of $names (names of states to call model functions)  and $all_names: all the names of states returned by model function
#' @export
#'
names_states <- function(model_type){
  if (model_type=="SEIRDItotNB"){
    names <- c("S","E","I","R","D")
    all_names <- c("S","E","I","R","D","Itot")
  } else if (model_type=="SEIFRDItotNB"){
    names <- c("S","E","I","Fu","R","D")
    all_names <- c("S","E","I","Fu","R","D","Itot","Dinst")
  } else if (model_type=="SEIHRDItotNB"){
    names <- c("S", "E", "I", "H", "R", "RH", "D", "DH")
    all_names <- c("S", "E", "I", "H", "R", "RH", "D", "DH", "Itot", "Dtot")
  } else if (model_type=="SEIHFRDItotNB"){
    names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D")
    all_names <- c("S", "E", "I", "H", "Fu", "R", "RH", "D", "Itot", "Dinst")
  }

  return(list(names=names,
              all_names=all_names))
}
