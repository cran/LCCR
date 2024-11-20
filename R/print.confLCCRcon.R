print.confLCCRcon = function(x, ...){

#---- print output of confLCCRcon ---
  cat("\nConfidence interval for the populations size based on latent class models and CML method")
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
