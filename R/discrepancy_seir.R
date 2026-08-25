#' @title Simulate constraints and calculate discrepancy function for SEIR model
#' @description
#' Simulate constraints and calculate discrepancy function. For now, there are no constraints
#' @param parameters a vector of parameters.
#' @param args a list of arguments as defined in the function \link[SMCfeatures]{SMCfeatures}
#' @return a list with
#' simulation list like object with t_const, y_const, K.
#' a summary statistic (discrepancy measure).
#' @export
#'
#'
discrepancy_seir <- function (parameters, args){
  simulation <- list()
  simulation$t_const <- 0
  simulation$y_const <- 0

  if (args$include_expert_constraints) {
    # maybe try to find how to change the parameters so that
    # max_population_in_5 <- simulation$y_const[1]
    # tolerance_in_5 <- 10
    discrepancy <- 0 #max(max_population_in_5 - tolerance_in_5, 0)

    return(list(simulation, discrepancy))
  }
  else {
    return(list(simulation, args$discrepancy_final))
  }
}
