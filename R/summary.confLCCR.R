summary.confLCCR <-
function(object, ...){

  cat("\nConfidence interval for the populations size based on latent class models")
  cat("\nfor capture recapture data\n")
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\nAvailable objects:\n")
  print(names(object))
  cat("\nLevl:\n")
  print(object$level)
  cat("\nInterval limits:\n")
  print(object$conf)

}
