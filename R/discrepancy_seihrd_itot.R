#' @title Simulate constraints and calculate discrepancy function for SEIR model
#' @description
#' Simulate constraints and calculate discrepancy function. For now, there are no constraints
#' @param parameters a vector of parameters.
#' @param args a list of arguments as defined in the function \link[SMCfeatures]{SMCfeatures}
#' @return a list with
#' simulation list like object with t_const, y_const, K.
#' a summary statistic (discrepancy measure).
#' @export

discrepancy_seihrd_itot <- function (parameters, args){
  simulation <- list()
  simulation$t_const <- args$input_data
  simulation$y_const <- SMCfeatures::model_seihrd_itot(parameters,simulation$t_const)

  if (args$include_expert_constraints) {
    #hospital model must be within data bounds

    discrepancy_low <- max(args$output_data$gaussian$H - simulation$y_const$H, 0) #data H must be lower than prediction

    discrepancy_high <- max(simulation$y_const$H - args$output_data$gaussian$total_H, 0) #prediction must be lower than total number of hospitalized

    discrepancy <- discrepancy_low+discrepancy_high

    return(list(simulation, discrepancy))
  }
  else {
    return(list(simulation, args$discrepancy_final))
  }
}
