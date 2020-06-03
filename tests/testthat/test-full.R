test_that("simulation runs with default inputs", {
  output <- run_simulation(model=list(start_time=1997, end_time=1998))
  expect_type(output, "list")
})
