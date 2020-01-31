test_that("simulation runs with default inputs", {
  basedir <- system.file('defaults', package = 'vivax')
  outfile <- tempfile()

  run_simulation_from_path(
    file.path(basedir, 'model_parameters.txt'),
    file.path(basedir, 'farauti_parameters.txt'),
    file.path(basedir, 'punctualtus_parameters.txt'),
    file.path(basedir, 'koliensis_parameters.txt'),
    file.path(basedir, 'intervention_coverage.txt'),
    file.path(outfile)
  )

  output <- read.table(outfile)
  expect_equal(dim(output), c(14600, 63))
})
