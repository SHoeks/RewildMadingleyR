writeCohortOutputBins = function(temp_dir, vector = c(0.0001,0.0010,0.0100,0.1000,1.0000,10.0000,100.0000,1000.0000,10000.0000)) {

  input_dir = paste0(temp_dir, "/input/")
  sink(paste0(input_dir, "CohortCohortOutputBins.csv"))
  for (i in vector) {
    cat(paste0(i))
    cat("\n")
  }
  sink()
}
