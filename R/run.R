
#'@title Run the simulation
#'@description parameterises and runs the model. A table with the simulation
#'output for each timestep will be returned.
#'
#'Default parameters can be found in the inst/default/ directory
#'
#'@param model, A named list of parameters to override in the
#'"model_parameters.txt" file
#'@param farauti, A named list of parameters to override in the
#'"farauti_parameters.txt" file
#'@param punctulatus, A named list of parameters to override in the
#'"punctulatus_parameters.txt" file
#'@param koliensis, A named list of parameters to override in the
#'"koliensis_parameters.txt" file
#'@param interventions, A named list of intervention parameters
#'@examples
#'output <- run_simulation(model = list(end_time = 1991, bb = 2), farauti=list(mu_M = .5))
#'dim(output)
#'
run_simulation <- function(
  model = NULL,
  farauti = NULL,
  punctulatus = NULL,
  koliensis = NULL,
  interventions = list()
  ) {
  basedir <- system.file('defaults', package = 'vivax', mustWork = TRUE)
  model_param_specs <- list(
    list(name='model_parameters.txt', overrides=model),
    list(name='farauti_parameters.txt', overrides=farauti),
    list(name='punctulatus_parameters.txt', overrides=punctulatus),
    list(name='koliensis_parameters.txt', overrides=koliensis)
  )

  tables <- lapply(model_param_specs, function(spec) {
    fixup(
      read.table(file.path(basedir, spec$name)),
      spec$overrides
    )
  })

  model_param_paths <- vapply(tables, function(table) {
    path <- tempfile()
    write.table(table, path, row.names = FALSE, col.names = FALSE)
    path
  }, character(1))

  intervention_param_path <- tempfile()
  generate_intervention_file(interventions, intervention_param_path)

  out_path <- tempfile()
  run_simulation_from_path(
    model_param_paths[[1]],
    model_param_paths[[2]],
    model_param_paths[[3]],
    model_param_paths[[4]],
    intervention_param_path
  )
}

#'@title Fixup a data table
#'@description overrides parameter entries in a table with alternatives from a list
#'@param table the parameter table to override
#'@param overrides the replacement parameters and values
fixup <- function(table, overrides) {
  if (is.null(overrides)) {
    return(table)
  }
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    matches <- table$V1 == name
    if (!any(matches)) {
      stop(paste('unknown parameter', name))
    }
    table[matches, 'V2'] <- overrides[[name]]
  }

  table
}
