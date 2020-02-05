test_that("simulation runs for a year", {
  output <- run_simulation(model=list(end_time=1991))
  expect_equal(dim(output), c(365, 63))
})
