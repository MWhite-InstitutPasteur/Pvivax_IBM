
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
#'@param prev_min_ages, A vector of minimum ages (in years) to disaggregate prevalence
#'statistics by.
#'@param prev_max_ages, A vector of maximum ages (in years) to disaggregate prevalence
#'statistics by. Individuals will be summarised by the interval:
#'[prev_min_ages[i], prev_max_ages[i])
#'@param incidence_min_ages, A vector of minimum ages (in years) to disaggregate incidence 
#'statistics by
#'@param incidence_max_ages, A vector of maximum ages (in years) to disaggregate incidence
#'statistics by. Individuals will be summarised by the interval:
#'[incidence_min_ages[i], incidence_max_ages[i])
#'@param a0 If set, the seasonality will be modelled with the fourier method and
#'a0 will be the first parameter
#'@param a_seasonality the a parameters for the fourier model of seasonality
#'@param b_seasonality the b parameters for the fourier model of seasonality
#'@examples
#'output <- run_simulation(model = list(end_time = 1991, bb = 2), farauti=list(mu_M = .5))
#'dim(output)
#'
run_simulation <- function(
  model = NULL,
  farauti = NULL,
  punctulatus = NULL,
  koliensis = NULL,
  interventions = list(),
  prev_min_ages = 2,
  prev_max_ages = 10,
  incidence_min_ages = c(0, 5,  15),
  incidence_max_ages = c(5, 15, 80),
  a0 = NULL,
  a_seasonality = NULL,
  b_seasonality = NULL
  ) {
  basedir <- system.file('defaults', package = 'vivax', mustWork = TRUE)
  model_param_specs <- list(
    list(name='model_parameters.txt', overrides=model),
    list(name='farauti_parameters.txt', overrides=farauti),
    list(name='punctulatus_parameters.txt', overrides=punctulatus),
    list(name='koliensis_parameters.txt', overrides=koliensis)
  )

  if (!all(vapply(c(
      prev_min_ages,
      prev_max_ages,
      incidence_min_ages,
      incidence_max_ages
    ), is.numeric, logical(1)))) {
    stop('summary ages must be numeric')
  }

  if (!(length(prev_min_ages) == length(prev_max_ages))) {
    stop('prevalence ages should match up')
  }

  if (!all(prev_min_ages < prev_max_ages)) {
    stop('prevalence min ages should be less than the max')
  }

  if (!(length(incidence_min_ages) == length(incidence_max_ages))) {
    stop('incidence ages should match up')
  }

  if (!all(incidence_min_ages < incidence_max_ages)) {
    stop('incidence min ages should be less than the max')
  }

  use_fourier <- !is.null(a0)
  if (use_fourier) {
    if (length(a_seasonality) == 0) {
      stop('seasonality vectors should be non-empty')
    }
    if (!(length(a_seasonality) == length(b_seasonality))) {
      stop('seasonality vectors should match up')
    }
  } else {
    a0 <- 0
    a_seasonality <- c(0)
    b_seasonality <- c(0)
  }

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
    intervention_param_path,
    prev_max_ages,
    prev_min_ages,
    incidence_max_ages,
    incidence_min_ages,
    use_fourier,
    a0,
    a_seasonality,
    b_seasonality
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
