
#'@title Run the simulation
#'@description parameterises and runs the model. The output will be returned
#'@param model parameters to override
#'@param farauti parameters to override
#'@param punctulatus parameters to override
#'@param koliensis parameters to override
run_simulation <- function(
  model = NULL,
  farauti = NULL,
  punctulatus = NULL,
  koliensis = NULL
  ) {
  basedir <- system.file('defaults', package = 'vivax')
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
  read.table(out_path)
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
