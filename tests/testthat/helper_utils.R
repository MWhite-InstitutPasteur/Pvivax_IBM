
expect_tableequal <- function(a, b) {
  row.names(a) <- NULL
  row.names(b) <- NULL
  expect_mapequal(a, b)
}
