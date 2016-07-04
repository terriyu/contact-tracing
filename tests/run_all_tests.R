library(testthat)

tryCatch(test_dir('tests/testthat'), error = function(e) {
  e$message <- paste0(e$message, "\nTry running this script in the root directory.")
    stop(e)
  })
