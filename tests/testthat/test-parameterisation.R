test_that("default parameters set up correctly", {
  defaultdir <- system.file('defaults', package = 'vivax', mustWork = TRUE)
  m <- mockery::mock()
  with_mock(
    `vivax::run_simulation_from_path` = function(...) {
      write.table(table(seq(10)), list(...)[[6]])
      m(...)
    },
    run_simulation()
  )
  args <- mockery::mock_args(m)
  expect_mapequal(
    read.table(file.path(defaultdir, 'model_parameters.txt')),
    read.table(args[[1]][[1]])
  )
  expect_mapequal(
    read.table(file.path(defaultdir, 'farauti_parameters.txt')),
    read.table(args[[1]][[2]])
  )
  expect_mapequal(
    read.table(file.path(defaultdir, 'punctulatus_parameters.txt')),
    read.table(args[[1]][[3]])
  )
  expect_mapequal(
    read.table(file.path(defaultdir, 'koliensis_parameters.txt')),
    read.table(args[[1]][[4]])
  )
  expect_mapequal(
    as.data.frame(matrix(-1, nrow=1, ncol=24)),
    read.table(args[[1]][[5]])
  )
})

test_that("you can override parameters", {
  defaultdir <- system.file('defaults', package = 'vivax', mustWork = TRUE)
  m <- mockery::mock()
  with_mock(
    `vivax::run_simulation_from_path` = function(...) {
      write.table(table(seq(10)), list(...)[[6]])
      m(...)
    },
    run_simulation(model = list(bb = 2), farauti=list(mu_M = .5))
  )
  args <- mockery::mock_args(m)
  default_model_params <- read.table(file.path(defaultdir, 'model_parameters.txt'))
  default_fara_params <- read.table(file.path(defaultdir, 'farauti_parameters.txt'))
  model_params <- read.table(args[[1]][[1]])
  fara_params <- read.table(args[[1]][[2]])
  expect_equal(
    model_params[model_params$V1 == 'bb', 'V2'],
    2
  )
  expect_equal(
    fara_params[fara_params$V1 == 'mu_M', 'V2'],
    .5
  )
  expect_mapequal(
    default_model_params[default_model_params$V1 != 'bb',],
    model_params[model_params$V1 != 'bb',]
  )
  expect_mapequal(
    default_fara_params[default_fara_params$V1 != 'mu_M',],
    fara_params[fara_params$V1 != 'mu_M',]
  )
})

test_that("overriding kindly lets you know if you've made a typo", {
  m <- mockery::mock()
  expect_error(
    with_mock(
      `vivax::run_simulation_from_path` = function(...) {
        write.table(table(seq(10)), list(...)[[6]])
        m(...)
      },
      run_simulation(model = list(bbb = 2), farauti=list(mu_M = .5))
    ),
    'unknown parameter bbb'
  )
})

test_that("you can model bed nets", {
  path <- tempfile()
  generate_intervention_file(
    list(
      LLIN_years = c(1995, 1996),
      LLIN_cover = c(.6, .9)
    ),
    path
  )
  expected <- as.data.frame(matrix(c(
    1995, .6, rep(-1, 22),
    1996, .9, rep(-1, 22),
    rep(-1, 24)
  ), nrow=3, ncol=24, byrow=TRUE))
  expect_mapequal(expected, read.table(path))
})

test_that("intervention file generator tells you if your vectors don't match up", {
  path <- tempfile()
  expect_error(
    generate_intervention_file(
      list(
        LLIN_years = c(1995, 1996),
        LLIN_cover = c(.6)
      ),
      path
    ),
    'LLIN vectors are not the same size'
  )
})
