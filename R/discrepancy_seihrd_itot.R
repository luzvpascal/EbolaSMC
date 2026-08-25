#' @title Simulate constraints and calculate discrepancy function for SEIHRD model
#' @description
#' Simulate constraints and calculate discrepancy function. For now, there are no constraints
#' @param parameters a vector of parameters.
#' @param args a list of arguments as defined in the function \link[EbolaSMC]{SMCfeatures}
#' @return a list with
#' simulation list like object with t_const, y_const, K.
#' a summary statistic (discrepancy measure).
#' @export

discrepancy_seihrd_itot <- function (parameters, simulation, args){

  simulation$t_const <- args$output_data$constraint$time #might not start at time 0
  t_index <- simulation$t_const+1
  simulation$y_const <- simulation$y_data[t_index,]

  if (args$include_expert_constraints) {
    #hospital model must be within data bounds
    discrepancy_low <- max(args$output_data$constraint$H - simulation$y_const$H, 0) #data H must be lower than prediction

    discrepancy_high <- max(simulation$y_const$H - args$output_data$constraint$total_H, 0) #prediction must be lower than total number of hospitalized

    discrepancy <- discrepancy_low+discrepancy_high

    return(list(simulation, discrepancy))
  }
  else {
    return(list(simulation, args$discrepancy_final))
  }
}
