
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
#'@examples
#'output <- run_simulation(model = list(end_time = 1991, bb = 2), farauti=list(mu_M = .5))
#'dim(output)
#'
run_simulation <- function(
  model = NULL,
  farauti = NULL,
  punctulatus = NULL,
  koliensis = NULL
  ) {
  basedir <- system.file('defaults', package = 'vivax', mustWork = TRUE)
  param_specs <- list(
    list(name='model_parameters.txt', overrides=model),
    list(name='farauti_parameters.txt', overrides=farauti),
    list(name='punctulatus_parameters.txt', overrides=punctulatus),
    list(name='koliensis_parameters.txt', overrides=koliensis),
    list(name='intervention_coverage.txt', overrides=NULL)
  )

  tables <- lapply(param_specs, function(spec) {
    fixup(
      read.table(file.path(basedir, spec$name)),
      spec$overrides
    )
  })

  paths <- vapply(tables, function(table) {
    path <- tempfile()
    write.table(table, path, row.names = FALSE, col.names = FALSE)
    path
  }, character(1))

  out_path <- tempfile()
  run_simulation_from_path(
    paths[[1]],
    paths[[2]],
    paths[[3]],
    paths[[4]],
    paths[[5]],
    out_path
  )
  present_output(read.table(out_path))
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
