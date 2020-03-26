fake_model_output <- function(outpath) {
  write.table(as.data.frame(matrix(1:63, nrow = 1, ncol = 63)), outpath, row.names=FALSE, col.names=FALSE)
}
