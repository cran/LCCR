print.confLCCR = function(x, ...){

#---- print output of confLCCR ---
  cat("\nConfidence interval for the populations size based on latent class models")
  cat("\nfor capture-recapture data\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\nAvailable objects:\n")
  print(names(x))
  cat("\nLevl:\n")
  print(x$level)
  cat("\nInterval limits:\n")
  print(x$conf)

}
